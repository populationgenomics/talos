"""
script for converting the 'extended PED format' to a regular PED file and Phenopackets file
"""

from argparse import ArgumentParser
from peds import open_ped

import phenopackets.schema.v2 as pps2
from google.protobuf.json_format import MessageToJson


# map the integer reported sex values to the enum
reported_sex_map = {1: pps2.Sex.MALE, 2: pps2.Sex.FEMALE}


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('ped_file', help='path to the extended PED file')
    parser.add_argument('output', help='stem path to write phenopacket and new PED file')
    args = parser.parse_args()

    extended_ped = open_ped(args.ped_file)

    # create the cohort shell
    cohort = pps2.Cohort(
        id='unknown',
        description=f'Phenotypic data from {args.ped_file}',
        members=[],
        meta_data=pps2.MetaData(
            created_by='Talos',
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

    for participant in extended_ped.members:
        # todo save stripped back pedigree entry

        # create a Phenopacket for this individual, and add it to the cohort
        cohort.members.append(
            pps2.Phenopacket(
                # primary ID is the internal CPG ID
                id=participant.id,
                subject=pps2.Individual(
                    # ID is the Family ID; Cohort model doesn't explicitly accommodate Family ID, but we need it
                    id=participant.family,
                    # alternate_ids captures the external ID, and is optional
                    alternate_ids=[participant.data[0]] if participant.data else None,
                    sex=reported_sex_map.get(participant.sex, pps2.Sex.UNKNOWN_SEX),
                ),
                # if HPO terms were provided, add them
                phenotypic_features=[
                    pps2.PhenotypicFeature(type=pps2.OntologyClass(id=hpo)) for hpo in participant.data[1:]
                ],
            )
        )

    # todo save both things


if __name__ == '__main__':
    cli_main()
