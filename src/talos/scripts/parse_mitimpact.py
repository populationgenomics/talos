"""
Takes the MitoTip ZIP file as input, outputs a new VCF file. Strictly uses MitImpact v4.


input columns: LOADS. THERE'S LOADS OF COLUMNS. Relevant columns:
Chr Start Ref Alt APOGEE2 APOGEE2_probability APOGEE2_score nAPOGEE_score nAPOGEE nAPOGEE_probability
"""

import gzip
import io
import zipfile
from argparse import ArgumentParser
from csv import DictReader
from importlib import resources


def _clean(value: str | None) -> str:
    value = (value or '').strip()
    return '' if value in ('', '.', 'NA', 'nan', 'NaN') else value


def main(input_mitimpact: str, output: str):
    # nested filehandles, so can't be consolidated to a single 'with'
    with zipfile.ZipFile(input_mitimpact, 'r') as mitimpact_opened:  # noqa: SIM117
        with (
            mitimpact_opened.open('MitImpact4_V059.csv', 'r') as handle,
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
            reader = DictReader(handle_data, delimiter=',')
            for line in reader:
                chrom = line['Chr']
                pos = line['Start']
                ref = line['Ref']
                alt = line['Alt']

                # pull details from the CSV
                content = {
                    'apogee2': _clean(line['APOGEE2']),
                    'apogee2_prob': _clean(line['APOGEE2_probability']),
                    'apogee2_score': _clean(line['APOGEE2_score']),
                    'napogee': _clean(line['nAPOGEE']),
                    'napogee_prob': _clean(line['nAPOGEE_probability']),
                    'napogee_score': _clean(line['nAPOGEE_score']),
                }

                # write that content as a String
                content_string = ';'.join('{k}={v}' for k, v in content.items())

                out.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t{content_string}\n')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', required=True, help='Input zipped mitimpact text file')
    parser.add_argument('--output', required=True, help='output VCF')
    args = parser.parse_args()
    main(input_mitimpact=args.input, output=args.output)
