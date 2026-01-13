"""
One script to entertain two different behaviours in the ClinvArbitration/Zenodo process

1. given a zenodo ID, we want to locate the latest version of the record
2. given the latest version of the record, we want to get the download URL for the record's file

This is split into two steps as (AFAIK) Nextflow doesn't handle dynamic output names particularly well, and once we find
the latest name for a file, we want to check if we already have that downloaded. That makes the process flow:

- start with a zenodo ID
- use the script to get the latest version of the same record
- do we already have a file called "clinvarbitration_<that ID>.tar.gz"
  - if yes, use it as an input channel
  - if no, rerun this step to get the URL for that file/content, download, and return as a channel
"""

from argparse import ArgumentParser

import httpx


def get_zenodo_record_details(record_id: int, download: bool = False) -> str:
    """
    Resolve information for a Zenodo record.

    Take a zenodo record ID, find the latest version of that record, return latest ID.
    Optionally if the behaviour switch `--download` is used, print out the file/content link.

    args:
        record_id (int): The numeric Zenodo record identifier to resolve.
        download (bool): If ``False``, resolve and return the latest record ID;
            if ``True``, return the download URL for the first file of the given record instead.
    return:
        A string containing either the latest Zenodo record ID or the file download URL
    """

    manual_retries = 3

    # use this record ID to find the content if we're downloading, or the latest version if we're searching
    record_url = f'https://zenodo.org/api/records/{record_id}'
    if not download:
        record_url += '/versions/latest'

    data = None
    while manual_retries > 0:
        manual_retries -= 1
        try:
            response = httpx.get(record_url, headers={'Accept': 'application/json'}, timeout=60, follow_redirects=True)
            data = response.json()
            if not download:
                return str(data['recid'])
        except httpx.HTTPError:
            # Failed to get latest zenodo record ID due to an HTTP/network error, retrying
            continue

    if data is None:
        raise httpx.DecodingError(f'Failed to get latest zenodo record ID: {record_id}, exhausted retries')

    # navigate to the files, and find the content URL we can use to download
    return str(data['files'][0]['links']['self'])


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('zenodo_id', type=int)
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    print(get_zenodo_record_details(args.zenodo_id, download=args.download))
