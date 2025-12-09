import json
import logging
import os
import re
import subprocess
import sys


def extract_version_from_file(file_path: str) -> str | None:
    """
    Extract the version from a Dockerfile by searching for a line like:
    ENV VERSION=1.0.0
    """
    with open(file_path) as f:
        content = f.read()
    pattern = re.compile(r'^\s*ENV\s+VERSION\s*=\s*([^\s]+)', re.MULTILINE)
    match = pattern.search(content)
    return match.group(1) if match else None


def get_next_version_tag(folder: str, version: str) -> str:
    """
    Query GCP to list tags for the given image and determine the next available
    version suffix for the extracted version.
    """
    base_image_path_prod = os.environ.get(
        'GCP_BASE_IMAGE',
        'australia-southeast1-docker.pkg.dev/cpg-common/images',
    )
    base_image_path_archive = os.environ.get(
        'GCP_BASE_ARCHIVE_IMAGE',
        'australia-southeast1-docker.pkg.dev/cpg-common/images-archive',
    )
    full_image_name_prod = f'{base_image_path_prod}/{folder}'
    full_image_name_archive = f'{base_image_path_archive}/{folder}'

    tags_list = []

    logging.basicConfig(level=logging.ERROR)

    for full_image_name in [full_image_name_prod, full_image_name_archive]:
        cmd = [
            'gcloud',
            'container',
            'images',
            'list-tags',
            full_image_name,
            '--format=json',
        ]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)  # noqa: S603
        except subprocess.CalledProcessError:
            logging.error(f'Failed to list tags for {full_image_name}')
            continue

        # If existing tags are found, proceed to determine the next version.
        tags_list += json.loads(result.stdout)

    max_suffix = 0
    pattern = re.compile(rf'^{re.escape(version)}-(\d+)$')

    # If no tags are found, it returns the next version as -1.
    for entry in tags_list:
        tags = entry.get('tags', [])
        for tag in tags:
            match = pattern.match(tag)
            if match:
                num = int(match.group(1))
                max_suffix = max(max_suffix, num)
    new_suffix = max_suffix + 1
    return f'{version}-{new_suffix}'


def main():
    dockerfile_name = 'Dockerfile'
    container_name = 'talos'
    current_version = extract_version_from_file(dockerfile_name)
    if current_version is None:
        # Throw an error here
        raise NotImplementedError('The Dockerfile needs to contain a version string in the format "ENV VERSION=x.x.x"')

    # Determine the next available tag based on current_version.
    new_tag = get_next_version_tag(container_name, current_version)

    include_entries = [{'name': container_name, 'tag': new_tag}]

    # Build the final matrix structure.
    matrix = {'include': include_entries}
    print(json.dumps(matrix, separators=(',', ':')), file=sys.stderr)
    print(json.dumps(matrix, separators=(',', ':')), end='')


main()
