"""
Methods for taking the final output and generating static report content
"""

import re
import sys
from argparse import ArgumentParser
from collections import defaultdict
from dataclasses import dataclass
from itertools import chain
from pathlib import Path
from typing import Any

import jinja2
import pandas as pd
from cloudpathlib.anypath import to_anypath

from talos.config import config_retrieve
from talos.models import PanelApp, PanelDetail, ReportVariant, ResultData, SmallVariant, StructuralVariant
from talos.utils import get_logger, read_json_from_path

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent / 'templates'

# above this length we trim the actual bases to just an int
MAX_INDEL_LEN: int = 10

# regex pattern - number, number, not-number
KNOWN_YEAR_PREFIX = re.compile(r'\d{2}\D')
CDNA_SQUASH = re.compile(r'(?P<type>ins|del)(?P<bases>[ACGT]+)$')


def cli_main():
    get_logger(__file__).info('Running HTML builder')
    parser = ArgumentParser()
    parser.add_argument('--results', help='Path to analysis results', required=True)
    parser.add_argument('--panelapp', help='PanelApp data', required=True)
    parser.add_argument('--output', help='Final HTML filename', required=True)
    parser.add_argument('--latest', help='Optional second report, latest variants only')
    parser.add_argument('--split_samples', help='divides samples into sub-reports', type=int)
    args = parser.parse_args()
    main(args.results, args.panelapp, args.output, args.latest, args.split_samples)


def main(results: str, panelapp: str, output: str, latest: str | None = None, split_samples: int | None = None):
    """

    Args:
        results (str): path to the MOI-tested results file
        panelapp (str): path to the panelapp data
        output (str): where to write the HTML file
        latest (str, optional): where to write a latest-results only file
        split_samples (int, optional): if this cohort should be subdivided into multiple reports
    """

    html = HTMLBuilder(results=results, panelapp_path=panelapp)
    # if this fails with a NoVariantsFoundException, there were no variants to present in the whole cohort
    # catch this, but fail gracefully so that the process overall is a success
    try:
        get_logger().info('Finding whole-cohort categorised variants')
        html.write_html(output_filepath=output)
    except NoVariantsFoundError:
        get_logger().warning('No Categorised variants found in this whole cohort')
        sys.exit(0)

    # If the latest arg is used, filter the results
    # write the HTML if any results remain
    if latest and (date_filtered_object := check_date_filter(results=results)):
        # build the HTML for latest reports only
        latest_html = HTMLBuilder(results=date_filtered_object, panelapp_path=panelapp)
        # this can fail if there are no latest-in-this-run variants, but we continue to splitting
        try:
            latest_html.write_html(output_filepath=latest, latest=True)
        except NoVariantsFoundError:
            get_logger().info('No latest-only variants found, but continuing on to subset splitting')

    # if no splitting, just exit here
    if not split_samples:
        get_logger().info('No splitting required in this run, exiting')
        sys.exit(0)

    # do something to split the output into separate datasets
    # either look for an ID convention, or go with a random split
    # originally this used Path(X).parent, but that translates gs:// to gs:/
    # gs:/ as a schema is not recognised as a GCP path, leading to write errors
    default_report_name = Path(output).name
    html_base = output.rstrip(default_report_name)

    for data, report, latest in split_data_into_sub_reports(results, split_samples):
        html = HTMLBuilder(results=data, panelapp_path=panelapp)
        try:
            get_logger().info(f'Attempting to create {report}')
            html.write_html(output_filepath=f'{html_base}{report}')

        except NoVariantsFoundError:
            get_logger().info('No variants in that report, skipping')

        # If the latest arg is used, filter the results
        # write the HTML if any results remain
        if latest and (date_filtered_object := check_date_filter(results=data)):
            # build the HTML for latest reports only
            latest_html = HTMLBuilder(results=date_filtered_object, panelapp_path=panelapp)
            try:
                get_logger().info(f'Attempting to create {latest_html}')
                latest_html.write_html(output_filepath=f'{html_base}{latest}', latest=True)
            except NoVariantsFoundError:
                get_logger().info('No variants in that latest report, skipping')


class NoVariantsFoundError(Exception):
    """raise if a report subset contains no data"""


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


def variant_in_forbidden_gene(variant_obj: ReportVariant, forbidden_genes):
    """
    Check if gene id or gene symbol is on forbidden gene list
    """
    for gene_id in variant_obj.gene.split(','):
        if gene_id in forbidden_genes:
            return True

    if isinstance(variant_obj.var_data, StructuralVariant):
        return False

    # Allow for exclusion by Symbol too
    return any(tx_con['symbol'] in forbidden_genes for tx_con in variant_obj.var_data.transcript_consequences)


class HTMLBuilder:
    """
    Takes the input, makes the output
    """

    def __init__(self, results: str | ResultData, panelapp_path: str):
        """
        Args:
            results (str | ResultData): path to the results JSON, or the results object
            panelapp_path (str): where to read panelapp data from
        """
        self.panelapp: PanelApp = read_json_from_path(panelapp_path, return_model=PanelApp)

        # If it exists, read the forbidden genes as a list
        self.forbidden_genes = config_retrieve(['GeneratePanelData', 'forbidden_genes'], [])
        assert isinstance(self.forbidden_genes, list)
        get_logger().warning(f'There are {len(self.forbidden_genes)} forbidden genes')

        # Use config to find CPG-to-Seqr ID JSON; allow to fail
        self.seqr: dict[str, str] = {}

        if seqr_path := config_retrieve(['CreateTalosHTML', 'seqr_lookup'], ''):
            self.seqr = read_json_from_path(seqr_path, default=self.seqr)
            assert isinstance(self.seqr, dict)

            # Force user to correct config file if seqr URL/project are missing
            for seqr_key in ['seqr_instance', 'seqr_project']:
                assert config_retrieve(['CreateTalosHTML', seqr_key], False), f'Seqr key absent: {seqr_key}'

        # Optionally read in the labels file
        # This file should be a nested dictionary of sample IDs and variant identifiers
        # with a list of corresponding label values, e.g.:
        # {
        #     "sample1": {
        #         "1-123456-A-T": ["label1", "label2"],
        #         "1-123457-A-T": ["label1"]
        #     },
        # }
        ext_labels = config_retrieve(['CreateTalosHTML', 'external_labels'], {})
        assert isinstance(ext_labels, dict)
        self.ext_labels: dict[str, dict] = ext_labels

        # Read results file, or take it directly
        results_dict = read_json_from_path(results, return_model=ResultData) if isinstance(results, str) else results

        assert isinstance(results_dict, ResultData)

        self.metadata = results_dict.metadata
        self.panel_names = {panel.name for panel in self.metadata.panels}

        # pull out forced panel matches
        cohort_panels = config_retrieve(['GeneratePanelData', 'forced_panels'], [])
        self.forced_panels: list = [panel for panel in self.metadata.panels if panel.id in cohort_panels]
        self.forced_panel_names = {panel.name for panel in self.metadata.panels if panel.id in cohort_panels}

        # Process samples and variants
        self.samples: list[Sample] = []
        self.solved: list[str] = []
        for sample, content in results_dict.results.items():
            if content.metadata.solved:
                self.solved.append(sample)
                continue
            self.samples.append(
                Sample(
                    name=sample,
                    metadata=content.metadata,
                    variants=content.variants,
                    ext_labels=self.ext_labels.get(sample, {}),
                    html_builder=self,
                ),
            )
        self.samples.sort(key=lambda x: x.ext_id)

    def get_summary_stats(self) -> tuple[pd.DataFrame, list[str], list[dict]]:
        """
        Run the numbers across all variant categories
        Treat each primary-secondary comp-het pairing as one event
        i.e. the thing being counted here is the number of events
        which passed through the MOI process, not the absolute number
        of variants in the report
        """
        ordered_categories = ['any', *list(config_retrieve('categories', {}).keys())]
        category_count: dict[str, list[int]] = {key: [] for key in ordered_categories}
        unique_variants: dict[str, set[str]] = {key: set() for key in ordered_categories}

        samples_with_no_variants: list[str] = []
        ext_label_map: dict = self.ext_labels.copy() if self.ext_labels else {}

        for sample in self.samples:
            if len(sample.variants) == 0:
                samples_with_no_variants.append(sample.ext_id)

            sample_variants: dict[str, set[str]] = {key: set() for key in ordered_categories}

            # iterate over the list of variants
            for variant in sample.variants:
                var_string = variant.var_data.coordinates.string_format
                unique_variants['any'].add(var_string)
                sample_variants['any'].add(var_string)

                # find all categories associated with this variant
                # for each category, add to corresponding list and set
                for category_value in variant.categories:
                    sample_variants[category_value].add(var_string)
                    unique_variants[category_value].add(var_string)

                # remove any external labels associated with this sample/variant.
                if sample.name in ext_label_map:
                    ext_label_map[sample.name].pop(var_string, None)

            # update the global lists with per-sample counts
            for key, key_list in category_count.items():
                key_list.append(len(sample_variants[key]))

        # Extract the list of unused ext labels
        unused_ext_labels = [
            {
                'sample': sample_id,
                'sample_ext': self.seqr.get(sample_id, sample_id),
                'variant': var_id,
                'labels': labels,
            }
            for sample_id, var_dict in ext_label_map.items()
            for var_id, labels in var_dict.items()
        ]

        summary_dicts = [
            {
                'Category': key,
                'Total': sum(category_count[key]),
                'Unique': len(unique_variants[key]),
                'Peak #/sample': max(category_count[key]),
                'Mean/sample': sum(category_count[key]) / len(category_count[key]),
            }
            for key in ordered_categories
            if category_count[key]
        ]

        # this can fail if there are no categorised variants... at all
        if not summary_dicts:
            raise NoVariantsFoundError('No categorised variants found')

        my_df: pd.DataFrame = pd.DataFrame(summary_dicts)
        my_df['Mean/sample'] = my_df['Mean/sample'].round(3)

        # the table re-sorts when parsed into the DataTable
        # so this forced ordering doesn't work
        my_df.Category = my_df.Category.astype('category')
        my_df.Category = my_df.Category.cat.set_categories(ordered_categories)
        my_df = my_df.sort_values(by='Category')

        return my_df, samples_with_no_variants, unused_ext_labels

    def read_metadata(self) -> dict[str, pd.DataFrame]:
        """
        parses into a general table and a panel table
        """

        tables = {
            'Panels': pd.DataFrame(
                {'ID': panel.id, 'Version': panel.version, 'Name': panel.name} for panel in self.metadata.panels
            ),
            'Meta': pd.DataFrame(
                {
                    'Data': key.capitalize(),
                    'Value': self.metadata.__getattribute__(key),
                }
                for key in ['run_datetime', 'version']
            ),
            'Families': pd.DataFrame(
                [
                    {'family_size': fam_type, 'tally': fam_count}
                    for fam_type, fam_count in sorted(self.metadata.family_breakdown.items())
                ],
            ),
        }

        if self.forced_panels:
            tables['Cohort Matched Panels'] = pd.DataFrame(self.forced_panels)

        return tables

    def write_html(self, output_filepath: str, latest: bool = False):
        """
        Uses the results to create the HTML tables
        writes all content to the output path

        Args:
            output_filepath (str): where to write the results to
            latest (bool):
        """

        # if no variants were found, this can fail with a NoVariantsFoundException error
        # we ignore that here, and catch it in the outer scope
        (summary_table, zero_cat_samples, unused_ext_labels) = self.get_summary_stats()

        report_title = 'Talos Report (Latest Variants Only)' if latest else 'Talos Report'

        template_context = {
            'metadata': self.metadata,
            'samples': self.samples,
            'seqr_url': config_retrieve(['CreateTalosHTML', 'seqr_instance'], ''),
            'seqr_project': config_retrieve(['CreateTalosHTML', 'seqr_project'], ''),
            'meta_tables': {},
            'forbidden_genes': sorted(self.forbidden_genes),
            'zero_categorised_samples': zero_cat_samples,
            'unused_ext_labels': unused_ext_labels,
            'summary_table': None,
            'report_title': report_title,
            'solved': self.solved,
        }

        for title, meta_table in self.read_metadata().items():
            template_context['meta_tables'][title] = DataTable(
                id=f'{title.lower()}-table',
                heading=title,
                description='',
                columns=list(meta_table.columns),
                rows=list(meta_table.to_records(index=False)),
            )

        template_context['summary_table'] = DataTable(
            id='summary-table',
            heading='Per-Category Summary',
            description='',
            columns=list(summary_table.columns),
            rows=list(summary_table.to_records(index=False)),
        )

        # write all HTML content to the output file in one go
        env = jinja2.Environment(loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR), autoescape=True)
        template = env.get_template('index.html.jinja')
        content = template.render(**template_context)
        to_anypath(output_filepath).open('w').writelines(
            '\n'.join(line for line in content.split('\n') if line.strip()),
        )


class Sample:
    """
    Sample related logic
    """

    def __init__(
        self,
        name: str,
        metadata,
        variants: list[ReportVariant],
        ext_labels: dict[str, list[str]],
        html_builder: HTMLBuilder,
    ):
        self.name = name
        self.family_id = metadata.family_id
        self.family_members = metadata.members
        self.phenotypes = metadata.phenotypes
        self.ext_id = metadata.ext_id
        self.panel_ids = metadata.panel_ids
        self.panel_names = metadata.panel_names
        self.seqr_id = html_builder.seqr.get(name, None)
        self.ext_labels = ext_labels
        self.html_builder = html_builder

        # Ingest variants excluding any on the forbidden gene list
        self.variants = [
            Variant(
                report_variant,
                self,
                ext_labels.get(report_variant.var_data.coordinates.string_format, []),
                html_builder.panelapp.genes,
            )
            for report_variant in variants
            if not variant_in_forbidden_gene(report_variant, html_builder.forbidden_genes)
        ]

    def __str__(self):
        return self.name


class Variant:
    """
    Handle as much of per variant logic as we can here. Hopefully, this is just simple
    data munging and mapping operations.

    Try not to put presentation logic here - keep it in the jinja templates
    """

    def get_var_change(self) -> str:
        """
        Find the variant change for the variant
        - we want to truncate huge small variant InDels (ballooning column width)
           - e.g. LOLOLOLOLOLOLOLOLOLOLOLOLOLOLOLOLO->A -> del 34bp
        - SVs always need to be presented differently
           - e.g INS 4079bp
        """
        if isinstance(self.var_data, SmallVariant):
            if len(self.ref) > MAX_INDEL_LEN or len(self.alt) > MAX_INDEL_LEN:
                ref_len = len(self.ref)
                alt_len = len(self.alt)
                if ref_len > alt_len:
                    return f'del {ref_len - alt_len}bp'
                if ref_len == alt_len:
                    return f'complex delins {ref_len}bp'
                return f'ins {alt_len - ref_len}bp'

            return f'{self.ref}->{self.alt}'
        if isinstance(self.var_data, StructuralVariant):
            return f"{self.var_data.info['svtype']} {self.var_data.info['svlen']}bp"

        raise ValueError(f'Unknown variant type: {self.var_data.__class__.__name__}')

    def __init__(
        self,
        report_variant: ReportVariant,
        sample: Sample,
        ext_labels: list,
        gene_map: dict[str, PanelDetail],
    ):
        self.var_data = report_variant.var_data
        self.var_type = report_variant.var_data.__class__.__name__
        self.chrom = report_variant.var_data.coordinates.chrom
        self.pos = report_variant.var_data.coordinates.pos
        self.ref = report_variant.var_data.coordinates.ref
        self.alt = report_variant.var_data.coordinates.alt
        self.change = self.get_var_change()
        self.categories = report_variant.categories
        self.first_tagged: str = report_variant.first_tagged
        self.support_vars = report_variant.support_vars
        self.warning_flags = report_variant.flags
        self.panel_flags = report_variant.panels.matched
        self.forced_matches = report_variant.panels.forced
        self.reasons = report_variant.reasons
        self.genotypes = report_variant.genotypes
        self.sample = sample
        self.ext_labels = ext_labels

        # List of (gene_id, symbol)
        self.genes: list[tuple[str, str]] = []
        for gene_id in report_variant.gene.split(','):
            symbol = gene_map.get(gene_id, PanelDetail(symbol=gene_id)).symbol
            self.genes.append((gene_id, symbol))

        # Summaries CSQ strings
        if isinstance(self.var_data, SmallVariant):
            (self.mane_csq, self.non_mane_csq, self.mane_hgvsps) = self.parse_csq()

        # pull up the highest AlphaMissense score, if present
        am_scores = (
            [
                float(csq['am_pathogenicity'])
                for csq in self.var_data.transcript_consequences
                if csq.get('am_pathogenicity')
            ]
            if isinstance(self.var_data, SmallVariant)
            else []
        )

        self.var_data.info['alpha_missense_max'] = max(am_scores) if am_scores else 'missing'

        # this is the weird gnomad callset ID
        if (
            isinstance(self.var_data, StructuralVariant)
            and 'gnomad_v2.1_sv_svid' in self.var_data.info
            and isinstance(self.var_data.info['gnomad_v2.1_sv_svid'], str)
        ):
            self.var_data.info['gnomad_key'] = self.var_data.info['gnomad_v2.1_sv_svid'].split('v2.1_')[-1]

    def __str__(self) -> str:
        return f'{self.chrom}-{self.pos}-{self.ref}-{self.alt}'

    def parse_csq(self):
        """
        Parse CSQ variant string returning:
            - set of "consequences" from MANE transcripts
            - set of "consequences" from non-MANE transcripts
            - Set of variant effects in p. nomenclature (or c. if no p. is available)

        condense massive cdna annotations, e.g.
        c.4978-2_4978-1insAGGTAAGCTTAGAAATGAGAAAAGACATGCACTTTTCATGTTAATGAAGTGATCTGGCTTCTCTTTCTA
        """
        mane_consequences = set()
        non_mane_consequences = set()
        mane_hgvsps = set()

        for csq in self.var_data.transcript_consequences:
            if 'consequence' not in csq:
                continue

            # if csq['mane_select'] or csq['mane_plus_clinical']:
            if csq['mane_select']:
                mane_consequences.update(csq['consequence'].split('&'))
                if csq['hgvsp']:
                    mane_hgvsps.add(csq['hgvsp'].split(':')[1])
                elif csq['hgvsc']:
                    hgvsc = csq['hgvsc'].split(':')[1]

                    # if massive indel base stretches are included, replace with a numerical length
                    if match := CDNA_SQUASH.search(hgvsc):
                        hgvsc.replace(match.group('bases'), str(len(match.group('bases'))))

                    mane_hgvsps.add(hgvsc)
            else:
                non_mane_consequences.add(csq['consequence'])

        return mane_consequences, non_mane_consequences, mane_hgvsps


def check_date_filter(results: str | ResultData, filter_date: str | None = None) -> ResultData | None:
    """
    Check if there's a date filter in the config
    if there is, load the results JSON and filter out variants

    Extra consideration - if one part of a comp-het variant pair is new,
    retain both sides in the report

    Args:
        results (str): path to the results file
        filter_date (str | None): path to the results file
    """

    # take both types
    if isinstance(results, str):
        # Load the results JSON
        results_dict: ResultData = read_json_from_path(results, return_model=ResultData)
    else:
        results_dict = results

    # pick up the current date from datetime or config
    if filter_date is None:
        filter_date = results_dict.metadata.run_datetime

    # Filter out variants based on date
    for content in results_dict.results.values():
        # keep only this run's new variants, or partners thereof
        vars_to_keep = [variant for variant in content.variants if variant.first_tagged == filter_date]

        pairs_to_keep = set(chain.from_iterable(var.support_vars for var in vars_to_keep))
        content.variants = [
            variant
            for variant in content.variants
            if (variant.first_tagged == filter_date or variant.var_data.coordinates.string_format in pairs_to_keep)
        ]

    # pop off all the samples with no variants
    for sample_id in list(results_dict.results.keys()):
        if not results_dict.results[sample_id].variants:
            results_dict.results.pop(sample_id)

    # check if there's anything to return
    if results_dict.results:
        get_logger().info(f'Filtered results obtained for {filter_date}')
        return results_dict

    get_logger().info(f'No filtered results obtained for {filter_date}')
    return None


def known_date_prefix_check(all_results: ResultData) -> list[str]:
    """
    Check for known date prefixes in the results

    Args:
        all_results (): the whole summary dataset

    Returns:
        a list of all found prefixes, or empty list
    """

    known_prefixes: dict[str, int] = defaultdict(int)
    for content in all_results.results.values():
        if match := KNOWN_YEAR_PREFIX.match(content.metadata.ext_id):
            known_prefixes[match.group()[0:2]] += 1
        else:
            get_logger().info(f'At least one sample lacks a consistent prefix: {content.metadata.ext_id}')
            return []

    get_logger().info(f'Sample distribution by prefix: {dict(known_prefixes)}')
    return sorted(known_prefixes.keys())


def split_data_into_sub_reports(data_path: str, split_samples: int) -> list[tuple[ResultData, str, str]]:
    """
    Split the data into sub-reports
    """
    all_results = read_json_from_path(data_path, return_model=ResultData)
    assert isinstance(all_results, ResultData)
    return_results: list[tuple[ResultData, str, str]] = []

    if prefixes := known_date_prefix_check(all_results):
        for prefix in prefixes:
            this_rd = ResultData(
                metadata=all_results.metadata,
                results={
                    sample: content
                    for sample, content in all_results.results.items()
                    if content.metadata.ext_id.startswith(prefix)
                },
                version=all_results.version,
            )
            get_logger().info(f'Found {len(this_rd.results)} with prefix {prefix}')
            return_results.append((this_rd, f'subset_{prefix}.html', f'subset_{prefix}_latest.html'))
        return return_results

    # calculate the number of samples per sub-report
    samples_per_report = len(all_results.results.keys()) // split_samples

    # split the data into sub-reports
    sub_reports = []
    for i in range(split_samples):
        start = i * samples_per_report
        end = (i + 1) * samples_per_report
        sub_report = ResultData(
            metadata=all_results.metadata,
            results=dict(list(all_results.results.items())[start:end]),
        )
        sub_reports.append((sub_report, f'subset_{i}.html', f'subset_{i}_latest.html'))

    return sub_reports


if __name__ == '__main__':
    cli_main()
