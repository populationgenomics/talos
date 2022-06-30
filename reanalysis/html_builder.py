"""
Methods for taking the final output and generating static report content
"""

import logging
from collections import defaultdict
from argparse import ArgumentParser
from typing import Any

from cloudpathlib import AnyPath
import pandas as pd
from peddy.peddy import Ped

from reanalysis.utils import read_json_from_path, good_string


GNOMAD_TEMPLATE = (
    '<a href="https://gnomad.broadinstitute.org/variant/'
    '{variant}?dataset=gnomad_r3" target="_blank">{value:.5f}</a>'
)
PANELAPP_TEMPLATE = (
    '<a href="https://panelapp.agha.umccr.org/panels/137/gene/{symbol}/" '
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

STRONG_STRING = '<strong>{content}</strong>'
COLOR_STRING = '<span style="color: {color}"><strong>{content}</strong></span>'
COLORS = {
    '1': '#FF0000',
    '2': '#FF9B00',
    '3': '1B00FF',
    'de_novo': '#FF0000',
    '5': '#006e4e',
    'support': '#00FF08',
}
CATEGORY_ORDERING = ['any', '1', '2', '3', 'de_novo', '5']


def category_strings(var_data: dict[str, Any], sample: str) -> list[str]:
    """
    get a list of strings representing the categories present on this variant
    :param var_data:
    :param sample:
    :return:
    """
    strings = [
        category
        for category in ['1', '2', '3', 'support', '5']
        if var_data.get(f'category_{category}', False)
    ]
    if sample in var_data.get('category_4', []):
        strings.append('de_novo')

    return strings


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
        mane_trans = each_csq.get('mane_select', '')
        if mane_trans != '':
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

    def __init__(
        self,
        results_dict: str,
        panelapp_data: str,
        config: str,
        pedigree: Ped,
    ):
        """

        :param results_dict:
        :param panelapp_data:
        :param config:
        :param pedigree:
        """
        self.results = read_json_from_path(results_dict)

        self.config = read_json_from_path(config)['output']

        # map of internal:external IDs for translation in results
        self.external_map = read_json_from_path(self.config['external_lookup']) or {}

        # use config to find CPG-to-Seqr ID JSON; allow to fail
        self.seqr = {}
        try:
            # update the seqr instance location
            if self.config.get('seqr_lookup'):
                self.seqr = read_json_from_path(self.config.get('seqr_lookup'))
                for seqr_key in ['seqr_instance', 'seqr_project']:
                    assert seqr_key in self.config and good_string(
                        self.config[seqr_key]
                    ), f'Seqr-related key required but not present: {seqr_key}'

        except AttributeError:
            logging.error(
                f'Failure parsing Seqr lookup from {self.config.get("seqr_lookup")}'
            )

        self.panelapp = read_json_from_path(panelapp_data)

        # read the ped file with Peddy
        self.pedigree = Ped(pedigree)

    def get_summary_stats(
        self,
    ) -> tuple[str, list[str]]:
        """
        run the numbers across all variant categories
        :return:
        """

        category_count = {key: [] for key in CATEGORY_ORDERING}
        unique_variants = defaultdict(set)

        samples_with_no_variants = []

        for sample, variants in self.results.items():

            if len(variants) == 0:
                samples_with_no_variants.append(self.external_map.get(sample, sample))

                # update all indices; 0 variants for this sample
                for category_list in category_count.values():
                    category_list.append(0)
                continue

            # how many variants were attached to this sample?
            # this set is for chr-pos-ref-alt
            # i.e. don't double-count if the variant is dominant and compound-het
            category_count['any'].append(
                len({key.split('__')[0] for key in variants.keys()})
            )

            # create a per-sample object to track variants for each category
            sample_count = defaultdict(int)

            # iterate over the variants
            for var_key, variant in variants.items():

                var_string = var_key.split('__')[0]

                # find all categories associated with this variant
                # for each category, add to corresponding list and set
                for category_value in category_strings(
                    variant['var_data'], sample=sample
                ):
                    if category_value == 'support':
                        continue
                    sample_count[category_value] += 1
                    unique_variants[category_value].add(var_string)

                # update the set of all unique variants
                unique_variants['any'].add(var_string)

            # update the global lists with per-sample counts
            for key, key_list in category_count.items():
                key_list.append(sample_count[key])

        summary_dicts = [
            {
                'Category': key,
                'Total': sum(category_count[key]),
                'Unique': len(unique_variants[key]),
                'Peak #/sample': max(category_count[key]),
                'Mean/sample': sum(category_count[key]) / len(category_count[key]),
            }
            for key in CATEGORY_ORDERING
        ]

        return (
            pd.DataFrame(summary_dicts).to_html(index=False, escape=False),
            samples_with_no_variants,
        )

    def write_html(self, output_path: str):
        """
        uses the results to create the HTML tables
        writes all content to the output path
        """

        summary_table, zero_categorised_samples = self.get_summary_stats()
        html_tables, category_2_genes = self.create_html_tables()
        category_2_table = self.category_2_table(category_2_genes)

        html_lines = ['<head>\n</head>\n<body>\n']

        if category_2_table:
            html_lines.extend(
                [
                    '<h3>MOI changes used for Cat.2</h3>',
                    category_2_table,
                    '<br/>',
                    f'<h3>Samples without Categorised Variants '
                    f'({len(zero_categorised_samples)})</h3>',
                ]
            )
        else:
            html_lines.append('<h3>No Cat.2 variants found</h3>')

        if len(zero_categorised_samples) > 0:
            html_lines.append(f'<h5>{", ".join(zero_categorised_samples)}</h3>')
        html_lines.append('<br/>')

        html_lines.append('<h3>Per-Category summary</h3>')
        html_lines.append(summary_table)
        html_lines.append('<br/>')

        html_lines.append('<h1>Per Sample Results</h1>')
        html_lines.append(
            '<br/>Note: "csq" shows all unique csq from all protein_coding txs'
        )
        html_lines.append('<br/>Any black "csq" appear on a MANE transcript')
        html_lines.append('<br/>Any red "csq" don\'t appear on a MANE transcript<br/>')

        for sample, table in html_tables.items():
            if sample in self.external_map and sample in self.seqr:
                sample_string = FAMILY_TEMPLATE.format(
                    seqr=self.config.get('seqr_instance'),
                    project=self.config.get('seqr_project'),
                    family=self.seqr[sample],
                    sample=self.external_map[sample],
                )
            else:
                sample_string = sample
            html_lines.append(fr'<h3>Sample: {sample_string}</h3>')
            html_lines.append(table)
        html_lines.append('\n</body>')

        # write all HTML content to the output file in one go
        # fix formatting later
        AnyPath(output_path).write_text(''.join(html_lines))

    def create_html_tables(self):
        """
        make a table describing the outputs
        :return:
        """
        candidate_dictionaries = {}
        sample_tables = {}

        category_2_genes = set()

        for sample, variants in self.results.items():
            for var_key, variant in variants.items():

                # pull out the string representation
                var_string = var_key.split('__')[0]

                # find list of all categories assigned
                variant_categories = category_strings(
                    variant['var_data'], sample=sample
                )

                if '2' in variant_categories:
                    category_2_genes.update(set(variant['gene'].split(',')))

                csq_string, mane_string = get_csq_details(variant)
                candidate_dictionaries.setdefault(variant['sample'], []).append(
                    {
                        'variant': self.make_seqr_link(
                            var_string=var_string, sample=sample
                        ),
                        'categories': ', '.join(
                            list(
                                map(
                                    lambda x: COLOR_STRING.format(
                                        color=COLORS[x], content=x
                                    ),
                                    variant_categories,
                                )
                            )
                        ),
                        # allow for multiple symbols on the same row
                        'symbol': ','.join(
                            [
                                PANELAPP_TEMPLATE.format(
                                    symbol=self.panelapp[symbol]['symbol']
                                )
                                for symbol in variant['gene'].split(',')
                            ]
                        ),
                        'csq': csq_string,
                        'mane_select': mane_string,
                        'gnomad': GNOMAD_TEMPLATE.format(
                            variant=var_string,
                            value=float(variant['var_data']['info']['gnomad_af']),
                        ),
                        'gnomad_AC': variant['var_data']['info']['gnomad_ac'],
                        'exac:hom': variant['var_data']['info']['exac_ac_hom'],
                        'gnomad_hom': variant['var_data']['info']['gnomad_hom'],
                        'gnomad_hemi': (
                            variant['var_data']['info']['gnomad_hemi']
                            if 'chrx' in var_string.lower()
                            else 'N/A'
                        ),
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
                    }
                )

        for sample, variants in candidate_dictionaries.items():
            sample_tables[sample] = pd.DataFrame(variants).to_html(
                index=False, render_links=True, escape=False
            )

        return sample_tables, category_2_genes

    def category_2_table(self, category_2_variants: set[str]) -> str:
        """
        takes all Cat. 2 variants, and documents relevant genes
        cat. 2 is now 'new genes', not 'new, or altered MOI'
        table altered to account for this changed purpose

        :param category_2_variants:
        :return:
        """

        if len(category_2_variants) == 0:
            return ''

        current_key = (
            f'MOI in v{self.panelapp["panel_metadata"].get("current_version")}'
        )

        gene_dicts = []
        for gene in category_2_variants:
            gene_data = self.panelapp.get(gene)
            gene_dicts.append(
                {
                    'gene': gene,
                    'symbol': PANELAPP_TEMPLATE.format(symbol=gene_data.get('symbol')),
                    current_key: gene_data.get('moi'),
                }
            )
        return pd.DataFrame(gene_dicts).to_html(
            index=False, render_links=True, escape=False
        )

    def make_seqr_link(self, var_string: str, sample: str) -> str:
        """
        either return just the variant as a string, or a seqr link if possible
        :param sample:
        :param var_string:
        :return:
        """
        if sample not in self.seqr:
            return var_string
        return SEQR_TEMPLATE.format(
            seqr=self.config.get('seqr_instance'),
            variant=var_string,
            family=self.seqr.get(sample),
        )


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument(
        '--results', help='Path to JSON containing analysis results', required=True
    )
    parser.add_argument(
        '--config_path', help='path to the runtime JSON config', required=True
    )
    parser.add_argument('--pedigree', help='Path to joint-call PED file', required=True)
    parser.add_argument(
        '--panelapp', help='Path to JSON file of PanelApp data', required=True
    )
    parser.add_argument(
        '--out_path', help='Path to write JSON results to', required=True
    )
    args = parser.parse_args()

    html_generator = HTMLBuilder(
        results_dict=args.results,
        panelapp_data=args.panelapp,
        pedigree=args.pedigree,
        config=args.config_path,
    )
    html_generator.write_html(output_path=args.out_path)
