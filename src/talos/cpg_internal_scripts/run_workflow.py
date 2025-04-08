"""
point of entry script to run this pipeline using cpg-flow
"""

from argparse import ArgumentParser

from talos.cpg_internal_scripts.talos_stages import MinimiseOutputForSeqr, UploadTalosHtml, CreateTalosHtml
from cpg_flow.workflow import run_workflow


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--dry-run', action='store_true', help='Print the commands that would be run')
    parser.add_argument('--nextflow', action='store_true', help='Run the nextflow workflow Stage')
    args = parser.parse_args()
    main(dry_run=args.dry_run)


def main(dry_run: bool = False):
    run_workflow(
        stages=[
            MinimiseOutputForSeqr,
            UploadTalosHtml,
            CreateTalosHtml,
        ],
        dry_run=dry_run,
    )


if __name__ == '__main__':
    cli_main()
