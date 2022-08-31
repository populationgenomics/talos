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

from reanalysis.utils import read_json_from_path


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


def make_coord_string(var_coord: dict[str, str]) -> str:
    """
    make a quick string representation from vardata
    """
    return (
        f'{var_coord["chrom"]}-{var_coord["pos"]}-{var_coord["ref"]}-{var_coord["alt"]}'
    )


def category_strings(var_data: dict[str, Any], sample: str) -> list[str]:
    """
    get a list of strings representing the categories present on this variant
    :param var_data:
    :param sample:
    :return:
    """
    strings = [
        cat.replace('categoryboolean', '')
        for cat in var_data['info'].keys()
        if cat.startswith('categoryboolean') and var_data['info'][cat]
    ]
    if var_data['info'].get('categorysupport'):
        strings.append('support')

    if sample in var_data['info'].get('categorysample4', []):
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
        before parsing data, purge any forbidden genes

        :param results_dict:
        :param panelapp_data:
        :param config:
        :param pedigree:
        """
        self.config = read_json_from_path(config)
        self.panelapp = read_json_from_path(panelapp_data)
        self.pedigree = Ped(pedigree)

        # if it exists, read the forbidden genes as a dict
        self.forbidden_genes = (
            {
                ensg: self.panelapp.get(ensg, {}).get('symbol', ensg)
                for ensg in set(
                    read_json_from_path(self.config['output'].get('forbidden'))
                )
            }
            if self.config['output'].get('forbidden') is not None
            else {}
        )

        logging.warning(f'There are {len(self.forbidden_genes)} forbidden genes')

        # pre-filter the results to remove forbidden genes
        self.results = self.remove_forbidden_genes(read_json_from_path(results_dict))

        # map of internal:external IDs for translation in results (optional)
        ext_lookup = self.config['output'].get('external_lookup')
        self.external_map = read_json_from_path(ext_lookup) if ext_lookup else {}

        # use config to find CPG-to-Seqr ID JSON; allow to fail
        seqr_path = self.config['output'].get('seqr_lookup')
        self.seqr = {}

        if seqr_path:
            self.seqr = read_json_from_path(seqr_path)

            # force user to correct config file if seqr URL/project are missing
            for seqr_key in ['seqr_instance', 'seqr_project']:
                assert self.config['output'].get(
                    seqr_key
                ), f'Seqr-related key required but not present: {seqr_key}'

    def remove_forbidden_genes(
        self, variant_dictionary: dict[str, Any]
    ) -> dict[str, Any]:
        """
        takes the results from the analysis and purges forbidden-gene variants
        """
        clean_results = defaultdict(list[dict[str, Any]])
        for sample, content in variant_dictionary.items():
            if sample == 'metadata':
                clean_results['metadata'] = content
                continue
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
    ) -> tuple[str, list[str]]:
        """
        run the numbers across all variant categories
        :return:
        """

        category_count = {key: [] for key in CATEGORY_ORDERING}
        unique_variants = {key: set() for key in CATEGORY_ORDERING}

        samples_with_no_variants = []

        for sample, variants in self.results.items():

            if sample == 'metadata':
                continue

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
                for category_value in category_strings(
                    variant['var_data'], sample=sample
                ):
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
        ]

        return (
            pd.DataFrame(summary_dicts).to_html(index=False, escape=False),
            samples_with_no_variants,
        )

    def read_metadata(self) -> dict[str, str]:
        """
        reads self.config[metadata]
        parses into a general table and a panel table
        """
        tables = {
            'Panels': pd.DataFrame(self.results['metadata']['panels']).to_html(
                index=False, escape=False
            ),
            'Meta': pd.DataFrame(
                {key.capitalize(): self.results['metadata'][key]}
                for key in ['cohort', 'run_datetime', 'input_file']
            ).to_html(index=False, escape=False),
            'Families': pd.DataFrame(
                [
                    {'family_size': fam_type, 'tally': fam_count}
                    for fam_type, fam_count in sorted(
                        self.results['metadata']['family_breakdown'].items()
                    )
                ]
            ).to_html(index=False, escape=False),
        }
        return tables

    def write_html(self, output_path: str):
        """
        uses the results to create the HTML tables
        writes all content to the output path
        """
        summary_table, zero_categorised_samples = self.get_summary_stats()
        html_tables = self.create_html_tables()

        html_lines = ['<head></head>\n<body>\n']

        for title, meta_table in self.read_metadata().items():
            html_lines.append(f'<h3>{title}</h3>')
            html_lines.append(meta_table)
            html_lines.append('<br/>')

        if self.forbidden_genes:
            # this should be sorted/arranged better
            forbidden_list = [
                f'{ensg} ({self.forbidden_genes[ensg]})'
                for ensg in self.forbidden_genes.keys()
            ]
            html_lines.append('<h3>Forbidden Gene IDs:</h3>')
            html_lines.append(f'<h4>{", ".join(forbidden_list)}</h4>')
        else:
            html_lines.append('<h3>No Forbidden Genes</h3>')
        html_lines.append('<br/>')

        if len(zero_categorised_samples) > 0:
            html_lines.append('<h3>Samples with no Reportable Variants</h3>')
            html_lines.append(f'<h4>{", ".join(zero_categorised_samples)}</h3>')
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
                    seqr=self.config['output'].get('seqr_instance'),
                    project=self.config['output'].get('seqr_project'),
                    family=self.seqr[sample],
                    sample=self.external_map[sample],
                )
            else:
                sample_string = self.external_map.get(sample, sample)

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

        for sample, variants in self.results.items():

            if sample == 'metadata':
                continue

            for variant in variants:

                # pull out the string representation
                var_string = make_coord_string(variant['var_data']['coords'])

                # find list of all categories assigned
                variant_categories = category_strings(
                    variant['var_data'], sample=sample
                )

                csq_string, mane_string = get_csq_details(variant)
                candidate_dictionaries.setdefault(variant['sample'], []).append(
                    {
                        'variant': self.make_seqr_link(
                            var_string=var_string, sample=sample
                        ),
                        'flags': ', '.join(variant['flags']),
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
                        'exomes_hom': variant['var_data']['info']['gnomad_ex_hom'],
                        'genomes_hom': variant['var_data']['info']['gnomad_hom'],
                        'gnomad_hemi': (
                            variant['var_data']['info']['gnomad_hemi']
                            if 'x' in var_string.lower()
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

        return sample_tables

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
            seqr=self.config['output'].get('seqr_instance'),
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
