#!/usr/bin/env python3

"""
This is an alternative to the current (as of 18/10/2024) approach in use for communicating phenotype data in cohorts

Talos has so far been built on a bespoke pedigree format, featuring arbitrarily extended columns to store
phenotypic and externalID information. This is a starting point, but to make support easier, we should move to a
standard format, such as the phenopackets format (https://phenopacket-schema.readthedocs.io/en/latest/index.html).

The phenopackets data interchange format is a GA4GH standard, and is being developed by the GA4GH community.

Here we represent our cohort as a Cohort object, consisting of multiple Members. Each Member is a Phenopacket,
encapsulating their phenotypic data and relevant ontological details.
"""

import re
from argparse import ArgumentParser
from collections import defaultdict

import phenopackets.schema.v2 as pps2
from google.protobuf.json_format import MessageToJson
from pyhpo import Ontology

from metamist.graphql import gql, query
from talos.static_values import get_logger

HPO_KEY = 'HPO Terms (present)'
HPO_RE = re.compile(r'HP:\d+')
PARTICIPANT_QUERY = gql(
    """
query MyQuery($project: String!, $sequencing_type: String!, $technology: String!) {
  project(name: $project) {
    pedigree
    sequencingGroups(technology: {eq: $technology}, type:  {eq: $sequencing_type}) {
      id
      sample {
        participant {
          externalId
          phenotypes
          reportedSex
          families {
            externalId
          }
        }
      }
    }
  }
}""",
)

# create an ontology object, once
_ = Ontology()

# map the integer reported sex values to the enum
reported_sex_map = {1: pps2.Sex.MALE, 2: pps2.Sex.FEMALE}


def find_hpo_labels(metamist_data: dict) -> dict[str, list[dict[str, str]]]:
    """
    match HPO terms to their plaintext names

    NB. we are making a decision here to strip out any participant HPO terms which are descendants of the term
    'Mode of Inheritance' (HP...5). Some clinicians may use this term to describe the suspected inheritance patterns
    for a familial disease, which can enable better variant curation. Downstream, we can discover associations between
    disease genes and their MOI terms directly, which messes with what we are trying to do:
        associate patient _phenotypes_ with variant _genes_.

    As we're now using hpo3/pyHPO which ships with a built-in ontology, we can use assume presence of a valid Ontology,
    simplifying the code

    Args:
        metamist_data ():

    Returns:
        dict, participant IDs to HPO:labels
    """
    per_sg_hpos: dict[str, list[dict[str, str]]] = defaultdict(list)

    for sg in metamist_data['project']['sequencingGroups']:
        # select all HPO terms, so long as they are not a child of 'Mode of Inheritance' (HP:0000005)
        hpos = {
            hpo_term
            for hpo_term in HPO_RE.findall(sg['sample']['participant']['phenotypes'].get(HPO_KEY, ''))
            if 'HP:0000005' not in Ontology.get_hpo_object(hpo_term).all_parents
        }

        # allow for HPO terms to be missing from this edition of the ontology
        for hpo in hpos:
            try:
                per_sg_hpos[sg['id']].append({'id': hpo, 'label': Ontology.get_hpo_object(hpo).name})
            except ValueError:
                get_logger(__file__).error(f'HPO term was absent from the tree: {hpo}')
                per_sg_hpos[sg['id']].append({'id': hpo, 'label': 'Unknown'})

    return per_sg_hpos


def assemble_phenopackets(dataset: str, metamist_data: dict, hpo_lookup: dict[str, list[dict[str, str]]]):
    """
    Assemble a cohort phenopacket from the metamist data

    Args:
        dataset (str): the dataset to query for
        metamist_data (dict): the data from metamist
        hpo_lookup (dict): a lookup of HPO terms to plaintext labels

    Returns:
        a cohort phenopacket
    """

    # create a cohort shell
    cohort = pps2.Cohort(
        id=dataset,
        description=f'Phenotypic data for {dataset}',
        members=[],
        meta_data=pps2.MetaData(
            created_by='Talos',
            # created=pps2.(datetime.today().isoformat()),
            resources=[
                pps2.Resource(
                    id='hp',
                    name='Human Phenotype Ontology',
                    url='http://www.human-phenotype-ontology.org',
                    version='2024-08-13',
                    namespace_prefix='HP',
                    iri_prefix='http://purl.obolibrary.org/obo/HP_',
                ),
            ],
        ),
    )

    # iterate over all SG Entities
    for sg in metamist_data['project']['sequencingGroups']:
        ext_id = sg['sample']['participant']['externalId']

        # allow for wild situation where we have no families
        if not sg['sample']['participant']['families']:
            family_id = ext_id
        else:
            family_id = sg['sample']['participant']['families'][0]['externalId']

        # create a Phenopacket for this individual
        ppack = pps2.Phenopacket(
            # primary ID is the internal CPG ID
            id=sg['id'],
            subject=pps2.Individual(
                # ID here is the Family ID; the Cohort model doesn't explicitly accommodate Family ID, but we need it
                id=family_id,
                # alternate_ids captures the external ID, and is optional
                alternate_ids=[ext_id],
                date_of_birth=sg['sample']['participant']['phenotypes'].get('Date of birth', None),
                sex=reported_sex_map.get(sg['sample']['participant']['reportedSex'], pps2.Sex.UNKNOWN_SEX),
            ),
            phenotypic_features=[
                pps2.PhenotypicFeature(type=pps2.OntologyClass(**each_hpo)) for each_hpo in hpo_lookup[sg['id']]
            ],
        )

        cohort.members.append(ppack)

    return cohort


def main(output: str, dataset: str, seq_type: str, tech: str = 'short-read'):
    """
    Assemble a cohort phenopacket from the metamist data
    Args:
        output (str): where to write the phenopacket
        dataset (str): the dataset to query for
        seq_type (str): exome/genome
        tech (str): type of sequence data to query for
    """

    # pull all the relevant data from metamist
    metamist_data = query(
        PARTICIPANT_QUERY,
        variables={'project': dataset, 'sequencing_type': seq_type, 'technology': tech},
    )

    # match names to HPO terms
    labelled_hpos = find_hpo_labels(metamist_data=metamist_data)

    # build the cohort
    cohort = assemble_phenopackets(dataset=dataset, metamist_data=metamist_data, hpo_lookup=labelled_hpos)

    with open(output, 'w', encoding='utf-8') as handle:
        json = MessageToJson(cohort)
        handle.write(json)


def cli_main():
    parser = ArgumentParser(description='Generate a PED file for Talos')
    parser.add_argument('--dataset', help='The dataset to query for')
    parser.add_argument('--output', help='The output file')
    parser.add_argument('--type', help='Sequencing type (exome or genome)')
    parser.add_argument('--tech', help='Sequencing technology', default='short-read')
    args = parser.parse_args()

    main(dataset=args.dataset, output=args.output, seq_type=args.type, tech=args.tech)


if __name__ == '__main__':
    cli_main()
