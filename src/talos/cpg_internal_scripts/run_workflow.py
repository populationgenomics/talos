"""
point of entry script to run this pipeline using cpg-flow
"""

from argparse import ArgumentParser

from cpg_flow import workflow

from talos.cpg_internal_scripts import talos_stages


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--dry-run', action='store_true', help='Print the commands that would be run')
    args = parser.parse_args()
    main(dry_run=args.dry_run)


def main(dry_run: bool = False):
    workflow.run_workflow(
        stages=[talos_stages.MinimiseOutputForSeqr, talos_stages.UpdateIndexFile],
        dry_run=dry_run,
    )


if __name__ == '__main__':
    cli_main()
