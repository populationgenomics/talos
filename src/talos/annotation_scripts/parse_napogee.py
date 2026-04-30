"""
takes the MitoTip tsv as input, outputs a new VCF file.

input columns: LOADS. THERE'S LOADS OF COLUMNS
relevant columns:
chr start ref alt nAPOGEE_score APOGEE_oob_score nAPOGEE_posterior_probability pathogenicity_assessment
"""

import gzip
import io
import zipfile
from argparse import ArgumentParser
from csv import DictReader
from importlib import resources


def main(input_napogee: str, output: str):
    with zipfile.ZipFile(input_napogee, 'r') as napogee_opened:
        filename = napogee_opened.namelist()[0]
        with (
            napogee_opened.open(filename, 'r') as handle,
            gzip.open(output, 'wt') as out,
            (resources.files('talos') / 'vcf_headers' / 'napogee_header.txt').open() as head_in,
        ):
            # transcribe the required header
            for header_line in head_in:
                out.write(header_line)

            # some real dancing around to get here:
            # a zipfile can contain multiple files, so we open the outer zip, then the inner file by name
            # that delivers bytes, so we need to pass that to a textWrapper to get at the contents
            # that in turn can be eaten as a DictReader to get to the true contents
            handle_data = io.TextIOWrapper(handle)
            reader = DictReader(handle_data, delimiter='\t')
            for line in reader:
                chrom = line['chr']
                pos = line['start']
                ref = line['ref']
                alt = line['alt']
                prob = line['nAPOGEE_posterior_probability']
                score = line['nAPOGEE_score']
                napogee = line['pathogenicity_assessment']
                clingen = line['ClinGen_classification_jan2025']

                out.write(
                    f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\tclingen={clingen};napogee_score={score};napogee={napogee};napogee_prob={prob}\n',
                )


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='Input zipped napogee text file')
    parser.add_argument('--output', help='output VCF')
    args = parser.parse_args()
    main(input_napogee=args.input, output=args.output)
