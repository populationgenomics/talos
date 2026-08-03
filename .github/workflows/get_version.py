import json
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
    cmd = [
        'gcloud',
        'container',
        'images',
        'list-tags',
        f'{base_image_path_prod}/{folder}',
        '--format=json',
    ]

    max_suffix = 0
    pattern = re.compile(rf'^{re.escape(version)}-(\d+)$')

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)  # noqa: S603
        for block in json.loads(result.stdout):
            if not block.get('tags', []):
                continue
            match = pattern.match(block['tags'][0])
            if match:
                num = int(match.group(1))
                max_suffix = max(max_suffix, num)
    except subprocess.CalledProcessError as err:
        raise RuntimeError('Failed to list tags for the given image') from err

    return f'{version}-{max_suffix + 1}'


def main():
    dockerfile_name = 'docker/Dockerfile'
    container_name = 'talos'
    current_version = extract_version_from_file(dockerfile_name)
    if current_version is None:
        # Throw an error here
        raise NotImplementedError('The Dockerfile needs to contain a version string in the format "ENV VERSION=x.x.x"')

    # Determine the next available tag based on current_version.
    new_tag = get_next_version_tag(container_name, current_version)

    print(new_tag, file=sys.stderr)
    print(new_tag, end='')


main()
