"""
Methods for taking the final output and generating static report content
"""

# pylint: disable=too-many-instance-attributes

import logging
import sys
from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import jinja2
import pandas as pd
from peddy.peddy import Ped

from cpg_utils import to_path
from cpg_utils.config import get_config
from reanalysis.utils import read_json_from_path


CATEGORY_ORDERING = ['any', '1', '2', '3', 'de_novo', '5', 'support']
JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent / 'templates'


@dataclass
class DataTable:
    """
    Representation of a DataTables table that the Jinja2 templating system renders.
    """

    id: str
    columns: list[str]
    rows: list[Any]
    heading: str = ''
    description: str = ''


def variant_in_forbidden_gene(variant_dict, forbidden_genes):
    """
    Check if gene id or gene symbol is on forbidden gene list
    """
    for gene_id in variant_dict['gene'].split(','):
        if gene_id in forbidden_genes:
            return True

    # Allow for exclusion by Symbol too
    for tx_con in variant_dict['var_data']['transcript_consequences']:
        if tx_con['symbol'] in forbidden_genes:
            return True

    return False


class HTMLBuilder:
    """
    Takes the input, makes the output
    """

    def __init__(self, results: str, panelapp: str, pedigree: Ped):
        """
        Args:
            results ():
            panelapp ():
            pedigree ():
        """
        self.panelapp = read_json_from_path(panelapp)
        self.pedigree = Ped(pedigree)

        # If it exists, read the forbidden genes as a set
        self.forbidden_genes = read_json_from_path(
            get_config()['dataset_specific'].get('forbidden', 'missing'), set()
        )

        logging.warning(f'There are {len(self.forbidden_genes)} forbidden genes')

        # Use config to find CPG-to-Seqr ID JSON; allow to fail
        seqr_path = get_config()['dataset_specific'].get('seqr_lookup')
        self.seqr = {}

        if seqr_path and to_path(seqr_path).exists():
            self.seqr = read_json_from_path(seqr_path)

            # Force user to correct config file if seqr URL/project are missing
            for seqr_key in ['seqr_instance', 'seqr_project']:
                assert get_config()['dataset_specific'].get(
                    seqr_key
                ), f'Seqr-related key required but not present: {seqr_key}'

        # Read results file
        results_dict = read_json_from_path(results)
        self.metadata = results_dict['metadata']
        self.panel_names = {panel['name'] for panel in self.metadata['panels']}

        # Process samples and variants
        self.samples = []
        for sample, content in results_dict['results'].items():
            self.samples.append(
                Sample(
                    name=sample,
                    metadata=content['metadata'],
                    variants=content['variants'],
                    html_builder=self,
                )
            )
        self.samples.sort(key=lambda x: x.ext_id)

    def get_summary_stats(
        self,
    ) -> tuple[pd.DataFrame, list[str]]:
        """
        Run the numbers across all variant categories
        :return:
        """

        category_count = {key: [] for key in CATEGORY_ORDERING}
        unique_variants = {key: set() for key in CATEGORY_ORDERING}

        samples_with_no_variants = []

        for sample in self.samples:

            if len(sample.variants) == 0:
                samples_with_no_variants.append(sample.ext_id)

            sample_variants = {key: set() for key in CATEGORY_ORDERING}

            # iterate over the list of variants
            for variant in sample.variants:
                var_string = str(variant)
                unique_variants['any'].add(var_string)
                sample_variants['any'].add(var_string)

                # find all categories associated with this variant
                # for each category, add to corresponding list and set
                for category_value in variant.var_data.get('categories'):
                    sample_variants[category_value].add(var_string)
                    unique_variants[category_value].add(var_string)

            category_count['any'].append(len(sample_variants['any']))

            # update the global lists with per-sample counts
            for key, key_list in category_count.items():
                key_list.append(len(sample_variants[key]))

        summary_dicts = [
            {
                'Category': key,
                'Total': sum(category_count[key]),
                'Unique': len(unique_variants[key]),
                'Peak #/sample': max(category_count[key]),
                'Mean/sample': sum(category_count[key]) / len(category_count[key]),
            }
            for key in CATEGORY_ORDERING
            if category_count[key]
        ]

        df = pd.DataFrame(summary_dicts)
        df['Mean/sample'] = df['Mean/sample'].round(3)

        # the table re-sorts when parsed into the DataTable
        # so this forced ordering doesn't work
        df.Category = df.Category.astype('category')
        df.Category = df.Category.cat.set_categories(CATEGORY_ORDERING)
        df = df.sort_values(by='Category')

        return df, samples_with_no_variants

    def read_metadata(self) -> dict[str, pd.DataFrame]:
        """
        parses into a general table and a panel table
        """

        tables = {
            'Panels': pd.DataFrame(self.metadata['panels']),
            'Meta': pd.DataFrame(
                {'Data': key.capitalize(), 'Value': self.metadata[key]}
                for key in ['cohort', 'input_file', 'run_datetime', 'container']
            ),
            'Families': pd.DataFrame(
                [
                    {'family_size': fam_type, 'tally': fam_count}
                    for fam_type, fam_count in sorted(
                        self.metadata['family_breakdown'].items()
                    )
                ]
            ),
        }
        return tables

    def write_html(self, output_path: str):
        """
        Uses the results to create the HTML tables
        writes all content to the output path
        """
        summary_table, zero_categorised_samples = self.get_summary_stats()

        template_context = {
            'metadata': self.metadata,
            'samples': self.samples,
            'seqr_url': get_config()['dataset_specific'].get('seqr_instance'),
            'seqr_project': get_config()['dataset_specific'].get('seqr_project'),
            'meta_tables': {},
            'forbidden_genes': [],
            'zero_categorised_samples': [],
            'summary_table': None,
        }

        for title, meta_table in self.read_metadata().items():
            template_context['meta_tables'][title] = DataTable(
                id=f'{title.lower()}-table',
                heading=title,
                description='',
                columns=list(meta_table.columns),
                rows=list(meta_table.to_records(index=False)),
            )

        if self.forbidden_genes:
            template_context['forbidden_genes'] = sorted(self.forbidden_genes)

        if len(zero_categorised_samples) > 0:
            template_context['zero_categorised_samples'] = zero_categorised_samples

        template_context['summary_table'] = DataTable(
            id='summary-table',
            heading='Per-Category Summary',
            description='',
            columns=list(summary_table.columns),
            rows=list(summary_table.to_records(index=False)),
        )

        # write all HTML content to the output file in one go
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR),
        )
        template = env.get_template('index.html.jinja')
        content = template.render(**template_context)
        to_path(output_path).write_text(
            '\n'.join(line for line in content.split('\n') if line.strip())
        )


class Sample:
    """
    Sample related logic
    """

    def __init__(
        self,
        name: str,
        metadata: dict,
        variants: list[dict[str, Any]],
        html_builder: HTMLBuilder,
    ):
        self.name = name
        self.family_id = metadata['family_id']
        self.family_members = metadata['members']
        self.phenotypes = metadata['phenotypes']
        self.ext_id = metadata['ext_id']
        self.panel_ids = metadata['panel_ids']
        self.panel_names = metadata['panel_names']
        self.seqr_id = html_builder.seqr.get(name, name)
        self.html_builder = html_builder

        # Ingest variants excluding any on the forbidden gene list
        self.variants = [
            Variant(variant_dict, self, html_builder.panelapp['genes'])
            for variant_dict in variants
            if not variant_in_forbidden_gene(variant_dict, html_builder.forbidden_genes)
        ]

    def __str__(self):
        return self.name


class Variant:
    """
    Handle as much of per variant logic as we can here. Hopefully, this is just simple
    data munging and mapping operations.

    Try not to put presentation logic here - keep it in the jinja templates
    """

    def __init__(
        self,
        variant_dict: dict,
        sample: Sample,
        gene_map: dict[str, Any],
    ):
        self.chrom = variant_dict['var_data']['coords']['chrom']
        self.pos = variant_dict['var_data']['coords']['pos']
        self.ref = variant_dict['var_data']['coords']['ref']
        self.alt = variant_dict['var_data']['coords']['alt']
        self.first_seen: str = variant_dict['first_seen']
        self.var_data = variant_dict['var_data']
        self.supported = variant_dict['supported']
        self.support_vars = variant_dict['support_vars']
        self.flags = variant_dict['flags']
        self.reasons = variant_dict['reasons']
        self.genotypes = variant_dict['genotypes']
        self.sample = sample

        # List of (gene_id, symbol)
        self.genes: list[tuple[str, str]] = []
        for gene_id in variant_dict['gene'].split(','):
            symbol = gene_map.get(gene_id, {'symbol': gene_id})['symbol']
            self.genes.append((gene_id, symbol))

        # Separate phenotype match flags from waring flags
        # TODO: should we keep these separate in the report?
        self.warning_flags = []
        self.panel_flags = []
        for flag in variant_dict['flags']:
            if flag in sample.html_builder.panel_names:
                self.panel_flags.append(flag)
            else:
                self.warning_flags.append(flag)

        # Summaries CSQ strings
        (
            self.mane_consequences,
            self.non_mane_consequences,
            self.mane_hgvsps,
        ) = self.parse_csq()

    def __str__(self) -> str:
        return f'{self.chrom}-{self.pos}-{self.ref}-{self.alt}'

    def parse_csq(self):
        """
        Parse CSQ variant string returning:
            - set of "consequences" from MANE transcripts
            - set of "consequences" from non-MANE transcripts
            - Set of variant effects in p. nomenclature (or c. if no p. is available)
        """
        mane_consequences = set()
        non_mane_consequences = set()
        mane_hgvsps = set()

        for csq in self.var_data.get('transcript_consequences', []):
            if 'consequence' not in csq:
                continue

            # if csq['mane_select'] or csq['mane_plus_clinical']:
            if csq['mane_select']:
                mane_consequences.update(csq['consequence'].split('&'))
                if csq['hgvsp']:
                    mane_hgvsps.add(csq['hgvsp'].split(':')[1])
                elif csq['hgvsc']:
                    mane_hgvsps.add(csq['hgvsc'].split(':')[1])
            else:
                non_mane_consequences.add(csq['consequence'])

        return mane_consequences, non_mane_consequences, mane_hgvsps


if __name__ == '__main__':

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    logging.info(f'Using templates in {JINJA_TEMPLATE_DIR}')

    parser = ArgumentParser()
    parser.add_argument('--results', help='Path to analysis results', required=True)
    parser.add_argument('--pedigree', help='PED file', required=True)
    parser.add_argument('--panelapp', help='PanelApp data', required=True)
    parser.add_argument('--out_path', help='results path', required=True)
    args = parser.parse_args()

    html = HTMLBuilder(
        results=args.results, panelapp=args.panelapp, pedigree=args.pedigree
    )
    html.write_html(output_path=args.out_path)
