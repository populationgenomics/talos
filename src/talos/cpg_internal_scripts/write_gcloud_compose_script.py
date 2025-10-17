""" """

import argparse
import os

import loguru

from cpg_utils import config

from talos import utils

COMPOSE_COMMAND = 'gcloud storage objects compose'


def make_compose_string(fragment_list: list[str], output: str) -> str:
    """
    Generate a gcloud compose command string to concatenate a list of VCF fragments
    Args:
        fragment_list ():
        output ():

    Returns:
        the compose instructions to generate the output from the fragments
    """

    return f'{COMPOSE_COMMAND} {" ".join(fragment_list)} {output}'


def main(input_manifest: str, vcf_dir: str, output_vcf: str, output_script: str, tempdir: str) -> None:
    """
    Generate a script to run gcloud compose commands based on the input manifest.

    Args:
        input_manifest (str): Path to the input manifest file containing VCF fragments.
        vcf_dir (str): directory containing VCF fragments
        output_vcf (str): Path to the final composition product.
        output_script (str): Path to the output script file to be generated.
        tempdir (str): Temporary GCP directory for intermediate files.
    """

    # please the linter
    vcf_fragments = []
    with open(input_manifest) as read_handle:
        vcf_fragments = [os.path.join(vcf_dir, line.strip()) for line in read_handle if line.strip()]

    if not vcf_fragments:
        raise ValueError('Input manifest is empty or contains no valid VCF fragments.')

    merge_round = 1
    condense_strings = []
    chunk_size = config.config_retrieve(['gcloud_condense', 'chunk_size'], 32)
    while len(vcf_fragments) > 1:
        new_fragments = []
        loguru.logger.info(f'Processing round {merge_round} with {len(vcf_fragments)} fragments.')

        # If we have more fragments than the chunk size, we need to merge them in chunks
        for merge_chunk, fragment_list in enumerate(utils.chunks(vcf_fragments, chunk_size)):
            output = f'{tempdir}/{merge_round}/temp_chunk_{merge_chunk}.vcf.gz'
            condense_strings.append(make_compose_string(fragment_list=fragment_list, output=output))
            new_fragments.append(output)

        vcf_fragments = new_fragments
        merge_round += 1

    # one left, _move_ it to the non-tmp bucket
    # inter-bucket condense operations aren't valid, so we can't 'compose' from tmp to main
    condense_strings.append(f'gcloud storage mv {vcf_fragments[0]} {output_vcf}')
    # Write the final script to the output file
    with open(output_script, 'w') as script_file:
        script_file.write('#!/bin/bash\n\n')
        for condense_string in condense_strings:
            script_file.write(f'{condense_string}\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate a script to run gcloud compose commands.')
    parser.add_argument('--input', required=True, help='Path to fragments and fragment manifest file.')
    parser.add_argument('--vcf_dir', required=True, help='Directory containing manifest/fragments.')
    parser.add_argument('--output', required=True, help='Output path for the final VCF.')
    parser.add_argument('--script', required=True, help='Output filename for the generated script.')
    parser.add_argument('--tmp', required=True, help='tmp directory for intermediate files.')
    args = parser.parse_args()
    main(
        input_manifest=args.input,
        vcf_dir=args.vcf_dir,
        output_vcf=args.output,
        output_script=args.script,
        tempdir=args.tmp,
    )
