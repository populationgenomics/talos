"""
Methods for taking the final output and generating static report content
"""

from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import jinja2
import pandas as pd
from peddy.peddy import Ped

from cpg_utils import to_path

from reanalysis.models import (
    PanelApp,
    PanelDetail,
    ParticipantMeta,
    ReportVariant,
    ResultData,
    SmallVariant,
    StructuralVariant,
)
from reanalysis.utils import (
    get_cohort_config,
    get_cohort_seq_type_conf,
    get_config,
    get_logger,
    read_json_from_path,
)

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent / 'templates'
DATASET_CONFIG = None  # type: ignore
DATASET_SEQ_CONFIG = None  # type: ignore


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
    for tx_con in variant_obj.var_data.transcript_consequences:
        if tx_con['symbol'] in forbidden_genes:
            return True

    return False


class HTMLBuilder:
    """
    Takes the input, makes the output
    """

    def __init__(self, results: str | ResultData, panelapp_path: str, pedigree: Ped):
        """
        Args:
            results (str | ResultData): path to the results JSON, or the results object
            panelapp_path (str): where to read panelapp data from
            pedigree (str): path to the PED file
        """
        self.panelapp: PanelApp = read_json_from_path(panelapp_path, return_model=PanelApp)  # type: ignore

        self.pedigree = Ped(pedigree)

        # If it exists, read the forbidden genes as a set
        self.forbidden_genes = read_json_from_path(
            DATASET_CONFIG.get('forbidden', 'missing'), default=set()  # type: ignore
        )

        get_logger().warning(f'There are {len(self.forbidden_genes)} forbidden genes')

        # Use config to find CPG-to-Seqr ID JSON; allow to fail
        seqr_path = DATASET_SEQ_CONFIG.get('seqr_lookup')  # type: ignore
        self.seqr: dict[str, str] = {}

        if seqr_path:
            self.seqr = read_json_from_path(seqr_path, default=self.seqr)  # type: ignore
            assert isinstance(self.seqr, dict)

            # Force user to correct config file if seqr URL/project are missing
            for seqr_key in ['seqr_instance', 'seqr_project']:
                assert DATASET_SEQ_CONFIG.get(  # type: ignore
                    seqr_key
                ), f'Seqr-related key required but not present: {seqr_key}'

        # Optionally read in the labels file
        # This file should be a nested dictionary of sample IDs and variant identifiers
        # with a list of corresponding label values, e.g.:
        # {
        #     "sample1": {
        #         "1-123456-A-T": ["label1", "label2"],
        #         "1-123457-A-T": ["label1"]
        #     },
        # }
        self.ext_labels: dict[str, dict] = read_json_from_path(  # type: ignore
            DATASET_SEQ_CONFIG.get('external_labels'), {}  # type: ignore
        )

        # Read results file, or take it directly
        if isinstance(results, str):
            results_dict = read_json_from_path(results, return_model=ResultData)  # type: ignore
        else:
            results_dict = results

        assert isinstance(results_dict, ResultData)

        self.metadata = results_dict.metadata
        self.panel_names = {panel.name for panel in self.metadata.panels}

        # pull out forced panel matches
        cohort_panels = DATASET_CONFIG.get('cohort_panels', [])  # type: ignore
        self.forced_panels = [
            panel for panel in self.metadata.panels if panel.id in cohort_panels
        ]
        self.forced_panel_names = {
            panel.name for panel in self.metadata.panels if panel.id in cohort_panels
        }

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
                )
            )
        self.samples.sort(key=lambda x: x.ext_id)

    def get_summary_stats(
        self,
    ) -> tuple[pd.DataFrame, list[str], list[dict]]:
        """
        Run the numbers across all variant categories
        Treat each primary-secondary comp-het pairing as one event
        i.e. the thing being counted here is the number of events
        which passed through the MOI process, not the absolute number
        of variants in the report
        """
        ordered_categories = ['any'] + list(get_config()['categories'].keys())
        category_count: dict[str, list[int]] = {key: [] for key in ordered_categories}
        unique_variants: dict[str, set[str]] = {
            key: set() for key in ordered_categories
        }

        samples_with_no_variants: list[str] = []
        ext_label_map: dict = self.ext_labels.copy() if self.ext_labels else {}

        for sample in self.samples:
            if len(sample.variants) == 0:
                samples_with_no_variants.append(sample.ext_id)

            sample_variants: dict[str, set[str]] = {
                key: set() for key in ordered_categories
            }

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

        df: pd.DataFrame = pd.DataFrame(summary_dicts)
        df['Mean/sample'] = df['Mean/sample'].round(3)

        # the table re-sorts when parsed into the DataTable
        # so this forced ordering doesn't work
        df.Category = df.Category.astype('category')
        df.Category = df.Category.cat.set_categories(ordered_categories)
        df = df.sort_values(by='Category')

        return df, samples_with_no_variants, unused_ext_labels

    def read_metadata(self) -> dict[str, pd.DataFrame]:
        """
        parses into a general table and a panel table
        """

        tables = {
            'Panels': pd.DataFrame(
                {'ID': panel.id, 'Version': panel.version, 'Name': panel.name}
                for panel in self.metadata.panels
            ),
            'Meta': pd.DataFrame(
                {
                    'Data': key.capitalize(),
                    'Value': self.metadata.__getattribute__(key),
                }
                for key in ['cohort', 'input_file', 'run_datetime', 'container']
            ),
            'Families': pd.DataFrame(
                [
                    {'family_size': fam_type, 'tally': fam_count}
                    for fam_type, fam_count in sorted(
                        self.metadata.family_breakdown.items()
                    )
                ]
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
            output_filepath ():
            latest (bool):
        """

        (summary_table, zero_cat_samples, unused_ext_labels) = self.get_summary_stats()

        report_title = 'AIP Report (Latest Variants Only)' if latest else 'AIP Report'

        template_context = {
            'metadata': self.metadata,
            'samples': self.samples,
            'seqr_url': DATASET_SEQ_CONFIG.get('seqr_instance', ''),  # type: ignore
            'seqr_project': DATASET_SEQ_CONFIG.get('seqr_project', ''),  # type: ignore
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
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR),
        )
        template = env.get_template('index.html.jinja')
        content = template.render(**template_context)
        to_path(output_filepath).write_text(
            '\n'.join(line for line in content.split('\n') if line.strip())
        )


class Sample:
    """
    Sample related logic
    """

    def __init__(
        self,
        name: str,
        metadata: ParticipantMeta,
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
            if not variant_in_forbidden_gene(
                report_variant, html_builder.forbidden_genes
            )
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
        report_variant: ReportVariant,
        sample: Sample,
        ext_labels: list,
        gene_map: dict[str, PanelDetail],
    ):
        self.var_type = report_variant.var_data.__class__.__name__
        self.chrom = report_variant.var_data.coordinates.chrom
        self.pos = report_variant.var_data.coordinates.pos
        self.ref = report_variant.var_data.coordinates.ref
        self.alt = report_variant.var_data.coordinates.alt
        self.categories = report_variant.categories
        self.first_seen: str = report_variant.first_seen
        self.var_data = report_variant.var_data
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

        self.var_data.info['alpha_missense_max'] = (
            max(am_scores) if am_scores else 'missing'
        )

        # this is the weird gnomad callset ID
        if isinstance(self.var_data, StructuralVariant):
            if 'gnomad_v2.1_sv_svid' in self.var_data.info:
                self.var_data.info['gnomad_key'] = self.var_data.info[  # type: ignore
                    'gnomad_v2.1_sv_svid'
                ].split('v2.1_')[-1]

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

        for csq in self.var_data.transcript_consequences:
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


def check_date_filter(
    results: str, filter_date: str | None = None
) -> ResultData | None:
    """
    Check if there's a date filter in the config
    if there is, load the results JSON and filter out variants

    Args:
        results (str): path to the results file
        filter_date (str | None): path to the results file
    """

    # Load the results JSON
    results_dict: ResultData = read_json_from_path(results, return_model=ResultData)  # type: ignore

    # pick up the current date from datetime or config
    if filter_date is None:
        filter_date = results_dict.metadata.run_datetime

    # Filter out variants based on date
    for content in results_dict.results.values():
        # keep only this run's new variants
        content.variants = [
            variant for variant in content.variants if variant.first_seen == filter_date
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


if __name__ == '__main__':
    get_logger(__file__).info('Running HTML builder')

    parser = ArgumentParser()
    parser.add_argument('--results', help='Path to analysis results', required=True)
    parser.add_argument('--pedigree', help='PED file', required=True)
    parser.add_argument('--panelapp', help='PanelApp data', required=True)
    parser.add_argument('--output', help='Final HTML filename', required=True)
    parser.add_argument('--latest', help='Optional second report, latest variants only')
    parser.add_argument('--dataset', help='Optional, dataset to use', default=None)
    args = parser.parse_args()

    DATASET_CONFIG = get_cohort_config(args.dataset)
    DATASET_SEQ_CONFIG = get_cohort_seq_type_conf(args.dataset)

    # build the HTML using all results
    html = HTMLBuilder(
        results=args.results, panelapp_path=args.panelapp, pedigree=args.pedigree
    )
    html.write_html(output_filepath=args.output)

    # If the latest arg is used, filter the results
    # write the HTML if any results remain
    if args.latest and (
        date_filtered_object := check_date_filter(results=args.results)
    ):
        # build the HTML for latest reports only
        latest_html = HTMLBuilder(
            results=date_filtered_object,
            panelapp_path=args.panelapp,
            pedigree=args.pedigree,
        )
        latest_html.write_html(output_filepath=args.latest, latest=True)
