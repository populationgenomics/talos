#!/usr/bin/env python3

"""
Script for generating a properly-formatted basic workflow TOML
file suitable for use with cpg-utils convenience functions 
(dataset_path(), remote_tmpdir(), etc.)

Requires the CPG_DEPLOY_CONFIG environment variable to be set to an
appropriate Json deployment config, or a pointer to a Json file
with the appropriate deployment config. Requires any referenced 
dataset to have a project mapping either in the server-config secret
of a full live CPG analysis-runner deployment or in the passed-in
server_config json argument.

python3 generate_workflow_config.py \
    --dataset severalgenomes \
    --access_level test \
    --driver_image image \
    --output_prefix gregsmi \
    --extra_datasets severalgenomes rgp \
    --deploy_config ~/sources/cpg/cpg-deploy/azure/deploy-config.prod.json \
    --server_config ~/sources/cpg/cpg-deploy/aip/terraform.tfvars.json \
    --image_base azcpg001acr.azurecr.io/cpg-common/images \
    --reference_base hail-az://azcpg001sa/reference \
    --extra_configs reanalysis_global.toml reanalysis_cohort.toml \
    --print_only
"""

import importlib.util as imp_util
import json
import os
from argparse import ArgumentParser
from typing import Any, Dict

import requests
import toml
from cpg_utils import to_path
from cpg_utils.config import read_configs, update_dict
from cpg_utils.deploy_config import (DeployConfig, get_deploy_config,
                                     set_deploy_config, set_server_config)
from cpg_utils.storage import get_dataset_bucket_config, get_dataset_bucket_url

IMAGES_FILE = "https://raw.githubusercontent.com/populationgenomics/images/main/images.toml"
REFERENCES_FILE = "https://raw.githubusercontent.com/populationgenomics/references/main/references.py"


def generate_image_section(image_base: str) -> Dict[str, str]:
    """
    Copied from https://github.com/populationgenomics/images/blob/main/.github/workflows/prep_config.py.
    Generates a set of absolute image paths from a registry prefix, based on the set in images.toml.
    """
    r = requests.get(IMAGES_FILE)
    d = toml.loads(r.text)
    # Manually adding cpg_workflows image.
    d["cpg_workflows"] = "cpg_workflows:latest"
    d = {k: f'{image_base}/{k}:{v}' for k, v in d.items()}
    return d


def generate_reference_section(reference_base: str) -> Dict[str, Any]:
    """
    Partially copied from https://github.com/populationgenomics/references/blob/main/.github/workflows/prep_config.py.
    Generates a set of absolute reference paths from a storage prefix, based on the set in references.py.
    """
    r = requests.get(REFERENCES_FILE)
    spec = imp_util.spec_from_loader('references', loader=None, origin=REFERENCES_FILE)
    references = imp_util.module_from_spec(spec)
    os.environ["PROJECT"] = "dummy-project"
    exec(r.content, references.__dict__)

    d = {'genome_build': references.GENOME_BUILD}

    for source in references.SOURCES:
        dst_path = f"{reference_base}/{source.dst}"
        if not source.files:
            d[source.name] = dst_path
        else:
            d[source.name] = { k: f"{dst_path}/{suffix}" for k, suffix in source.files.items() }
    return d


parser = ArgumentParser()
parser.add_argument("--dataset", help="name of primary dataset to target", required=True)
parser.add_argument("--access_level", help="whether to target test or main buckets", required=True)
parser.add_argument("--driver_image", help="fully-qualified Hail driver image name to target", required=True)
parser.add_argument("--output_prefix", help="subpath on dataset's access-level bucket to use for output", required=True)
parser.add_argument("--extra_datasets", help="additional datasets to add storage url blocks for", nargs="*")
parser.add_argument("--deploy_config", help="deployment configuration file to use during generation")
parser.add_argument("--server_config", help="deployment dataset configuration file to use during generation")
parser.add_argument("--image_base", help="base container registry path to use in creating an [images] section")
parser.add_argument("--reference_base", help="base storage path to use in creating a [references] section")
parser.add_argument("--extra_configs", help="additional workflow-specific toml files to merge in", nargs="*")
parser.add_argument("--print_only", help="print to screen, don't write to file", action="store_true")
parser.add_argument("-o", help="output file path and name (cloudpathlib.AnyPath-compatible paths are supported)")
args = parser.parse_args()

if args.deploy_config:
    with open(args.deploy_config, "r") as config_file:
        set_deploy_config(DeployConfig.from_dict(json.load(config_file)))

if args.server_config:
    with open(args.server_config, "r") as config_file:
        set_server_config(json.load(config_file)["datasets"])

assert args.access_level == "test" or args.access_level == "main"
print("generating workflow toml using deploy_config:")
print(json.dumps(get_deploy_config().to_dict(include_datasets=True), indent=2))

workflow_config = {
    'access_level': args.access_level,
    'dataset': args.dataset,
    'driver_image': args.driver_image,
    'output_prefix': args.output_prefix
}

hail_config = {
    "bucket": get_dataset_bucket_url(args.dataset, "hail"),
    "billing_project": args.dataset
}

storage_config = { 
    ds : get_dataset_bucket_config(ds, args.access_level) 
    for ds in [args.dataset] + (args.extra_datasets or [])
}

config = {
    "hail": hail_config,
    "workflow": workflow_config,
    "storage": storage_config,
    "CPG_DEPLOY_CONFIG": get_deploy_config().to_dict()
}

if args.image_base:
    config["images"] = generate_image_section(args.image_base)

if args.reference_base:
    config["references"] = generate_reference_section(args.reference_base)

if args.extra_configs:
    extra_config = read_configs(args.extra_configs)
    update_dict(config, extra_config)

if args.print_only:
    print(toml.dumps(config))
else:
    output_name = args.o or f"{args.dataset}_{args.access_level}.toml"
    print(f"writing config toml to {output_name}")
    config_path = to_path(output_name)
    with config_path.open('w') as f:
        toml.dump(config, f)
