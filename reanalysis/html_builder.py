"""
Methods for taking the final output and generating static report content
"""

from typing import Any, Dict, List, Set, Tuple
from argparse import ArgumentParser
from cloudpathlib import AnyPath
import pandas as pd
from peddy.peddy import Ped

from reanalysis.utils import read_json_dict_from_path


SEQR_TEMPLATE = (
    '<a href="https://seqr.populationgenomics.org.au/variant_search/'
    'variant/{variant}/family/{family}" target="_blank">{variant}</a>'
)

GNOMAD_TEMPLATE = (
    '<a href="https://gnomad.broadinstitute.org/variant/'
    '{variant}?dataset=gnomad_r3" target="_blank">{value:.5f}</a>'
)
PANELAPP_TEMPLATE = (
    '<a href="https://panelapp.agha.umccr.org/panels/137/gene/{symbol}/"'
    ' target="_blank">{symbol}</a>'
)

STRONG_STRING = '<strong>{content}</strong>'
COLOR_STRING = '<span style="color: {color}"><strong>{content}</strong></span>'


def numerical_categories(var_data: Dict[str, Any], sample: str) -> List[str]:
    """
    get a list of strings representing the categories present on this variant
    :param var_data:
    :param sample:
    :return:
    """
    strings = [
        category
        for category in ['1', '2', '3', 'support']
        if var_data[f'category_{category}']
    ]
    if sample in var_data['category_4']:
        strings.append('de_novo')

    return strings


def string_format_coords(coords: Dict[str, Any]):
    """
    forms a string representation
    chr-pos-ref-alt
    """
    return f'{coords["chrom"]}-{coords["pos"]}-{coords["ref"]}-{coords["alt"]}'


def set_up_colors():
    """
    generates the dictionary of color strings
    colours are hard coded, rather than parsed from config
    """

    return {
        '1': '#FF0000',
        '2': '#FF9B00',
        '3': '1B00FF',
        'de_novo': '#FF0000',
        'support': '#00FF08',
    }


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
        self.results = read_json_dict_from_path(results_dict)
        self.config = read_json_dict_from_path(config)['output']

        # map of internal:external IDs for translation in results
        self.external_map = (
            read_json_dict_from_path(self.config['external_lookup']) or {}
        )

        # use config to find CPG-to-Seqr ID JSON; allow to fail
        try:
            self.seqr = read_json_dict_from_path(self.config.get('seqr_lookup'))
        except AttributeError:
            self.seqr = {}

        self.panelapp = read_json_dict_from_path(panelapp_data)

        # peddy can't read cloudpaths, so create a local copy
        with open('i_am_a_temporary.ped', 'w', encoding='utf-8') as handle:
            handle.write(AnyPath(pedigree).read_text())
        self.pedigree = Ped('i_am_a_temporary.ped')
        self.colors = set_up_colors()

    def get_summary_stats(
        self,
    ) -> Tuple[str, List[str]]:  # pylint: disable=too-many-locals
        """
        run the numbers across all variant categories
        :return:
        """

        category_count = {
            '1': [],
            '2': [],
            '3': [],
            'de_novo': [],
            'all': [],
        }
        category_strings = {
            '1': set(),
            '2': set(),
            '3': set(),
            'de_novo': set(),
            'support': set(),
            'all': set(),
        }

        samples_with_no_variants = []

        for sample, variants in self.results.items():

            # skip over any sample entries that
            if not self.pedigree[sample].affected:
                assert (
                    len(variants) == 0
                ), f'Sample {sample} was unaffected, but had categorised variants'
                continue

            if len(variants) == 0:
                samples_with_no_variants.append(self.external_map.get(sample, sample))

                # update all indices; 0 variants for this sample
                for category_list in category_count.values():
                    category_list.append(0)
                continue

            # how many variants were attached to this sample?
            # this set is for chr-pos-ref-alt
            # i.e. don't double-count if the variant is dominant and compound-het
            category_count['all'].append(
                len({key.split('__')[0] for key in variants.keys()})
            )

            # create a per-sample object to track variants for each category
            sample_count = {'1': 0, '2': 0, '3': 0, 'de_novo': 0}

            # iterate over the variants
            for var_key, variant in variants.items():

                var_string = var_key.split('__')[0]

                # find all categories associated with this variant
                # for each category, add to corresponding list and set
                for category_value in numerical_categories(
                    variant['var_data'], sample=sample
                ):
                    if category_value == 'support':
                        continue
                    sample_count[category_value] += 1
                    category_strings[category_value].add(var_string)

                # update the set of all unique variants
                category_strings['all'].add(var_string)

            # update the global lists with per-sample counts
            for key, value in sample_count.items():
                category_count[key].append(value)

        summary_dicts = [
            {
                'Category': title,
                'Total': sum(category_count[key]),
                'Unique': len(category_strings[key]),
                'Peak #/sample': max(category_count[key]),
                'Mean/sample': sum(category_count[key]) / len(category_count[key]),
            }
            for title, key in [
                ('Total', 'all'),
                ('Cat1', '1'),
                ('Cat2', '2'),
                ('Cat3', '3'),
                ('de_novo', 'de_novo'),
            ]
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
            html_lines.append(
                fr'<h3>Sample: {self.external_map.get(sample, sample)}</h3>'
            )
            html_lines.append(table)
        html_lines.append('\n</body>')

        # write all HTML content to the output file in one go
        # fix formatting later
        AnyPath(output_path).write_text(''.join(html_lines))

    def color_csq(self, all_csq: Set[str], mane_csq: Set[str]) -> str:
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
                csq_strings.append(
                    COLOR_STRING.format(color=self.colors['1'], content=csq)
                )

        return ', '.join(csq_strings)

    def get_csq_details(self, variant: Dict[str, Any]) -> Tuple[str, str]:
        """
        populates a single string of all relevant consequences
        UPDATE - take MANE into account
        """

        csq_set = set()
        mane_transcript = set()
        mane_csq = set()

        # iterate over all consequences, special care for MANE
        for each_csq in variant['var_data']['transcript_consequences']:
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

        return self.color_csq(csq_set, mane_csq), mane_string

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
                variant_categories = numerical_categories(
                    variant['var_data'], sample=sample
                )

                if '2' in variant_categories:
                    category_2_genes.add(variant['gene'])

                csq_string, mane_string = self.get_csq_details(variant)
                candidate_dictionaries.setdefault(variant['sample'], []).append(
                    {
                        'variant': self.make_seqr_link(
                            var_string=var_string, sample=sample
                        ),
                        'categories': ', '.join(
                            list(
                                map(
                                    lambda x: COLOR_STRING.format(
                                        color=self.colors[x], content=x
                                    ),
                                    variant_categories,
                                )
                            )
                        ),
                        'symbol': PANELAPP_TEMPLATE.format(
                            symbol=self.panelapp[variant['gene']]['symbol']
                        ),
                        'csq': csq_string,
                        'mane_select': mane_string,
                        'gnomad': GNOMAD_TEMPLATE.format(
                            variant=var_string,
                            value=float(variant['var_data']['info']['gnomad_af']),
                        ),
                        'gnomad_AC': variant['var_data']['info']['gnomad_ac'],
                        'exac:hom': variant['var_data']['info']['exac_ac_hom'],
                        'g_exome:hom': variant['var_data']['info']['gnomad_ex_hom'],
                        'g_genome:hom': variant['var_data']['info']['gnomad_hom'],
                        'MOIs': ', '.join(variant['reasons']),
                        'support': self.make_seqr_link(
                            var_string=var_string,
                            sample=sample,
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

    def category_2_table(self, category_2_variants: Set[str]) -> str:
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
