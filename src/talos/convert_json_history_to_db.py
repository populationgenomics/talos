import json
from argparse import ArgumentParser

from cloudpathlib.anypath import to_anypath

from talos.utils import (
    INSERT_CAT_QUERY,
    INSERT_PARTNER_QUERY,
    INSERT_VAR_QUERY,
    create_or_open_db,
    translate_category,
)

# one query for all 4 insertions
INSERT_VARSTARS = 'INSERT INTO var_stars (var_id, sample_id, clinvar_stars, first_pheno_match) VALUES (?, ?, ?, ?);'


def main(input_file: str, output_file: str):
    local_file = 'local.db'
    connection = create_or_open_db(local_file)
    cursor = connection.cursor()

    # empty db
    all_variants: dict[str, int] = {}

    with to_anypath(input_file).open() as handle:
        history = json.load(handle)

    # do the things here
    for sample_id, var_dict in history['results'].items():
        for var_key, content in var_dict.items():
            if var_key in all_variants:
                var_id = all_variants[var_key]

            else:
                chrom, pos, ref, alt = var_key.split('-')
                pos = int(pos)
                cursor.execute(
                    INSERT_VAR_QUERY,
                    (chrom, pos, ref, alt),
                )
                maybe_var_id = cursor.lastrowid
                if isinstance(maybe_var_id, int):
                    var_id = maybe_var_id
                    all_variants[var_key] = var_id
                else:
                    connection.close()
                    raise ValueError(f'Failed to retrieve lastrowid after inserting new variant: {var_key}')

            if content['support_vars']:
                cursor.execute(INSERT_PARTNER_QUERY, (var_id, sample_id, ','.join(sorted(content['support_vars']))))

            clinvar_stars = content.get('clinvar_stars', None)
            first_pheno_match = content.get('first_pheno_match', None)

            connection.execute(INSERT_VARSTARS, (var_id, sample_id, clinvar_stars, first_pheno_match))

            for category, date_string in content['categories'].items():
                if category.lower() == 'support':
                    continue
                cursor.execute(INSERT_CAT_QUERY, (var_id, sample_id, translate_category(category), date_string))
        connection.commit()
    connection.close()

    with open(local_file, 'rb') as db_in, to_anypath(output_file).open('wb') as db_out:
        db_out.write(db_in.read())


if __name__ == '__main__':
    parser = ArgumentParser(description='Convert a JSON history file to a SQLite database.')
    parser.add_argument('json_file', type=str, help='Path to the JSON history file.')
    parser.add_argument('db_file', type=str, help='Path to the output SQLite database file.')
    args = parser.parse_args()

    main(args.json_file, args.db_file)


"""
def find_latest_file(results_folder: str, ext: str = 'json') -> str | None:

    logger.info(f'Using results from {results_folder}')

    # this is currently a CloudPath to access globbing for files in cloud or local settings
    date_files = {}
    for filename in to_anypath(results_folder).glob(f'*.{ext}'):
        # if the filename is a valid date, add it to the dict, otherwise catch the error and skip it
        try:
            date_files[date_from_string(filename.name)] = filename
        except ValueError:
            logger.info(f'File {filename} did not have a valid date, skipping')

    if not date_files:
        logger.warning(f'The folder {results_folder} was provided, but did not contain any valid files')
        return None

    return str(date_files[max(date_files.keys())].absolute())

"""
