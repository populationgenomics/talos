"""
takes the MitoTip tsv as input, outputs a new VCF file.

input columns: LOADS. THERE'S LOADS OF COLUMNS
relevant columns:
Chr Start Ref Alt APOGEE2 APOGEE2_probability APOGEE2_score
"""

import gzip
import io
import zipfile
from argparse import ArgumentParser
from csv import DictReader
from importlib import resources


def main(input_mitimpact: str, output: str):
    with zipfile.ZipFile(input_mitimpact, 'r') as mitimpact_opened:
        filename = mitimpact_opened.namelist()[0]
        with (
            mitimpact_opened.open(filename, 'r') as handle,
            gzip.open(output, 'wt') as out,
            (resources.files('talos') / 'vcf_headers' / 'mitimpact_header.txt').open() as head_in,
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
                chrom = line['Chr']
                pos = line['Start']
                ref = line['Ref']
                alt = line['Alt']
                apogee2_prob = line['APOGEE2_probability']
                apogee2_score = line['APOGEE2_score']
                apogee2 = line['APOGEE2']

                out.write(
                    f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\tapogee2={apogee2};apogee2_score={apogee2_score};apogee2_probability={apogee2_prob}\n',
                )


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='Input zipped mitimpact text file')
    parser.add_argument('--output', help='output VCF')
    args = parser.parse_args()
    main(input_mitimpact=args.input, output=args.output)
