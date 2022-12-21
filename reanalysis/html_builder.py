"""
Methods for taking the final output and generating static report content
"""

# pylint: disable=too-many-instance-attributes

import logging
import sys
from collections import defaultdict
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


JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent / 'templates'

GNOMAD_TEMPLATE = (
    '<a href="https://gnomad.broadinstitute.org/variant/'
    '{variant}?dataset=gnomad_r3" target="_blank">{value:.5f}</a>'
)
PANELAPP_TEMPLATE = (
    '<a href="https://panelapp.agha.umccr.org/panels/entities/{symbol}" '
    'target="_blank">{symbol}</a>'
)
SEQR_TEMPLATE = (
    '<a href="{seqr}/variant_search/variant/{variant}/family/{family}" '
    'target="_blank">{variant}</a>'
)
FAMILY_TEMPLATE = (
    '<a href="{seqr}/project/{project}/family_page/{family}" '
    'target="_blank">{sample}</a>'
)
TOOLTIP_TEMPLATE = (
    '<a href="https://panelapp.agha.umccr.org/panels/{panelid}" '
    'data-toggle="tooltip" title="{panelname}" target="_blank">{display}</a>'
)

STRONG_STRING = '<strong>{content}</strong>'
COLOR_STRING = '<span style="color: {color}"><strong>{content}</strong></span>'
COLORS = {
    '1': '#FF0000',
    '2': '#FF9B00',
    '3': '#1B00FF',
    'de_novo': '#FF0000',
    '5': '#006e4e',
    'support': '#00FF08',
}
CATEGORY_ORDERING = ['any', '1', '2', '3', 'de_novo', '5']


@dataclass
class DataTable:
    """
    Representation of a DataTables table that the Jinja2 templating system renders.
    """

    id: str
    columns: list[str]
    rows: list[Any]
    heading: str = None
    description: str = None


def make_coord_string(var_coord: dict[str, str]) -> str:
    """
    make a quick string representation from vardata
    """
    return (
        f'{var_coord["chrom"]}-{var_coord["pos"]}-{var_coord["ref"]}-{var_coord["alt"]}'
    )


def color_csq(all_csq: set[str], mane_csq: set[str]) -> str:
    """
    takes the collection of all consequences, and MANE csqs
    if a CSQ occurs on MANE, write in bold,
    if non-MANE, write in red
    return the concatenated string

    NOTE: I really hate how I've implemented this
    :param all_csq:
    :param mane_csq:
    :return: the string filling the consequence box in the HTML
    """
    csq_strings = []
    for csq in all_csq:
        # bold, in Black
        if csq in mane_csq:
            csq_strings.append(STRONG_STRING.format(content=csq))

        # bold, and red
        else:
            csq_strings.append(COLOR_STRING.format(color=COLORS['1'], content=csq))

    return ', '.join(csq_strings)


def get_csq_details(variant: dict[str, Any]) -> tuple[str, str]:
    """
    populates a single string of all relevant consequences
    UPDATE - take MANE into account
    """

    csq_set = set()
    mane_transcript = set()
    mane_csq = set()

    # iterate over all consequences, special care for MANE
    for each_csq in variant['var_data'].get('transcript_consequences', []):

        # allow for variants with no transcript CSQs
        if 'consequence' not in each_csq:
            continue

        row_csq = set(each_csq['consequence'].split('&'))

        # record the transcript ID(s), and CSQ(s)
        # we only expect 1, but set operations are useful
        mane_trans = each_csq.get('mane_select')
        if mane_trans is not None:
            mane_csq.update(row_csq)
            mane_transcript.add(mane_trans)
        csq_set.update(row_csq)

    # we expect one... but this would make the addition of MANE plus_clinical easy
    mane_string = STRONG_STRING.format(content=', '.join(mane_transcript))

    return color_csq(csq_set, mane_csq), mane_string


class HTMLBuilder:
    """
    takes the input, makes the output
    """

    def __init__(self, results: str, panelapp: str, pedigree: Ped):
        """
        pass
        Args:
            results ():
            panelapp ():
            pedigree ():
        """
        self.panelapp = read_json_from_path(panelapp)
        self.pedigree = Ped(pedigree)

        # if it exists, read the forbidden genes as a dict
        self.forbidden_genes = (
            {
                ensg: self.panelapp['genes'].get(ensg, {}).get('symbol', ensg)
                for ensg in set(
                    read_json_from_path(
                        get_config().get('dataset_specific', {})['forbidden']
                    )
                )
            }
            if get_config().get('dataset_specific', {}).get('forbidden')
            else {}
        )

        logging.warning(f'There are {len(self.forbidden_genes)} forbidden genes')

        results_dict = read_json_from_path(results)
        self.metadata = results_dict['metadata']

        # pre-filter the results to remove forbidden genes
        self.variants = self.remove_forbidden_genes(results_dict['results'])

        self.panel_tooltips = self.generate_tooltips()

        # map of internal:external IDs for translation in results (optional)
        ext_lookup = get_config().get('dataset_specific', {}).get('external_lookup')
        self.external_map = read_json_from_path(ext_lookup) if ext_lookup else {}

        # use config to find CPG-to-Seqr ID JSON; allow to fail
        seqr_path = get_config().get('dataset_specific', {}).get('seqr_lookup')
        self.seqr = {}

        if seqr_path:
            self.seqr = read_json_from_path(seqr_path)

            # force user to correct config file if seqr URL/project are missing
            for seqr_key in ['seqr_instance', 'seqr_project']:
                assert (
                    get_config().get('dataset_specific', {}).get(seqr_key)
                ), f'Seqr-related key required but not present: {seqr_key}'

    def generate_tooltips(self) -> dict:
        """
        generate the per-panel tooltip lookups from self.metadata

        Returns:
            a dictionary of each panel and the corresponding tooltip
        """

        panel_dict = {}

        # iterate over the list of panels
        for panel in self.metadata['panels']:
            panel_dict[panel['id']] = {
                'star': TOOLTIP_TEMPLATE.format(
                    panelid=panel['id'],
                    panelname=f'{panel["name"]} - {panel["version"]}',
                    display='*',
                ),
                'name': TOOLTIP_TEMPLATE.format(
                    panelid=panel['id'],
                    panelname=f'{panel["name"]} - {panel["version"]}',
                    display='*',
                ),
            }

        return panel_dict

    def remove_forbidden_genes(
        self, variant_dictionary: dict[str, Any]
    ) -> dict[str, Any]:
        """
        takes the results from the analysis and purges forbidden-gene variants
        """
        clean_results = defaultdict(list)
        for sample, content in variant_dictionary.items():
            sample_vars = []
            for variant in content:
                skip_variant = False
                for gene_id in variant['gene'].split(','):
                    if gene_id in self.forbidden_genes.keys():
                        skip_variant = True
                        continue
                if not skip_variant:
                    sample_vars.append(variant)

            # add any retained variants to the new per-sample list
            clean_results[sample] = sample_vars

        return clean_results

    def get_summary_stats(
        self,
    ) -> tuple[pd.DataFrame, list[str]]:
        """
        run the numbers across all variant categories
        :return:
        """

        category_count = {key: [] for key in CATEGORY_ORDERING}
        unique_variants = {key: set() for key in CATEGORY_ORDERING}

        samples_with_no_variants = []

        for sample, variants in self.variants.items():

            if len(variants) == 0:
                samples_with_no_variants.append(self.external_map.get(sample, sample))

            sample_variants = {key: set() for key in CATEGORY_ORDERING}

            # iterate over the list of variants
            for variant in variants:
                var_string = make_coord_string(variant['var_data']['coords'])
                unique_variants['any'].add(var_string)
                sample_variants['any'].add(var_string)

                # find all categories associated with this variant
                # for each category, add to corresponding list and set
                for category_value in variant['var_data'].get('categories'):
                    if category_value == 'support':
                        continue
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

        # swap those panels out for hyperlinks

        for panel in self.metadata['panels']:
            panel['name'] = self.panel_tooltips[panel['id']]['name']

        tables = {
            'Panels': pd.DataFrame(self.metadata['panels']),
            'Meta': pd.DataFrame(
                {'Data': key.capitalize(), 'Value': self.metadata[key]}
                for key in ['cohort', 'input_file', 'run_datetime', 'commit_id']
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
        uses the results to create the HTML tables
        writes all content to the output path
        """
        summary_table, zero_categorised_samples = self.get_summary_stats()
        html_tables = self.create_html_tables()

        template_context = {
            'meta_tables': [],
            'forbidden_genes': [],
            'zero_categorised_samples': [],
            'summary_table': None,
            'sample_tables': [],
        }

        for title, meta_table in self.read_metadata().items():
            template_context['meta_tables'].append(
                DataTable(
                    id=f'{title.lower()}-table',
                    heading=title,
                    description='',
                    columns=list(meta_table.columns),
                    rows=list(meta_table.to_records(index=False)),
                )
            )

        if self.forbidden_genes:
            # this should be sorted/arranged better
            forbidden_list = [
                f'{ensg} ({self.forbidden_genes[ensg]})'
                for ensg in self.forbidden_genes.keys()
            ]
            template_context['forbidden_genes'] = forbidden_list

        if len(zero_categorised_samples) > 0:
            template_context['zero_categorised_samples'] = zero_categorised_samples

        template_context['summary_table'] = DataTable(
            id='summary-table',
            heading='Per-Category Summary',
            description='',
            columns=list(summary_table.columns),
            rows=list(summary_table.to_records(index=False)),
        )

        for sample, table in html_tables.items():
            if sample in self.external_map and sample in self.seqr:
                sample_string = FAMILY_TEMPLATE.format(
                    seqr=get_config().get('dataset_specific', {}).get('seqr_instance'),
                    project=get_config()
                    .get('dataset_specific', {})
                    .get('seqr_project'),
                    family=self.seqr[sample],
                    sample=self.external_map[sample],
                )
            else:
                sample_string = self.external_map.get(sample, sample)

            tid = f'{sample_string}-variants-table'
            table = DataTable(
                id=tid,
                heading=f'Sample: {sample_string}',
                description='',
                columns=list(table.columns),
                rows=list(table.to_records(index=False)),
            )
            template_context['sample_tables'].append(table)

        # write all HTML content to the output file in one go
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR),
        )
        template = env.get_template('index.html.jinja')
        content = template.render(**template_context)
        to_path(output_path).write_text(
            '\n'.join(line for line in content.split('\n') if line.strip())
        )

    def create_html_tables(self) -> dict[str, pd.DataFrame]:
        """
        make a table describing the outputs
        :return:
        """
        candidate_dictionaries = {}
        sample_tables = {}

        for sample, variants in self.variants.items():

            for variant in variants:

                # pull out the string representation
                var_string = make_coord_string(variant['var_data']['coords'])

                csq_string, mane_string = get_csq_details(variant)
                candidate_dictionaries.setdefault(variant['sample'], []).append(
                    {
                        'variant': self.make_seqr_link(
                            var_string=var_string, sample=sample
                        ),
                        'flags': ' '.join(
                            [
                                self.panel_tooltips[flag]['star']
                                for flag in variant['flags']
                                if flag in self.panel_tooltips
                            ]
                        ),
                        'categories': ', '.join(
                            list(
                                map(
                                    lambda x: COLOR_STRING.format(
                                        color=COLORS[x], content=x
                                    ),
                                    variant['var_data']['categories'],
                                )
                            )
                        ),
                        # allow for multiple symbols on the same row
                        'symbol': ','.join(
                            [
                                PANELAPP_TEMPLATE.format(
                                    symbol=self.panelapp['genes'][symbol]['symbol']
                                )
                                for symbol in variant['gene'].split(',')
                            ]
                        ),
                        'CSQ': csq_string,
                        'MOIs': ', '.join(variant['reasons']),
                        'support': ', '.join(
                            [
                                self.make_seqr_link(
                                    var_string=partner,
                                    sample=sample,
                                )
                                for partner in variant['support_vars']
                            ]
                        )
                        if variant['supported']
                        else 'N/A',
                        'mane_select': mane_string,
                        'gnomad': GNOMAD_TEMPLATE.format(
                            variant=var_string,
                            value=float(variant['var_data']['info']['gnomad_af']),
                        ),
                        'gnomad_AC': variant['var_data']['info']['gnomad_ac'],
                        'exomes_hom': variant['var_data']['info']['gnomad_ex_hom'],
                        'genomes_hom': variant['var_data']['info']['gnomad_hom'],
                        'gnomad_hemi': (
                            variant['var_data']['info']['gnomad_hemi']
                            if 'x' in var_string.lower()
                            else 'N/A'
                        ),
                    }
                )

        for sample, variants in candidate_dictionaries.items():
            sample_tables[sample] = pd.DataFrame(variants)

        return sample_tables

    def make_seqr_link(self, var_string: str, sample: str) -> str:
        """
        either return just the variant as a string, or a seqr link if possible

        Args:
            var_string ():
            sample ():

        Returns:
            a formatting string if possible, else just the variant
        """
        if sample not in self.seqr:
            return var_string
        return SEQR_TEMPLATE.format(
            seqr=get_config().get('dataset_specific', {}).get('seqr_instance'),
            variant=var_string,
            family=self.seqr.get(sample),
        )


if __name__ == '__main__':

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

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
