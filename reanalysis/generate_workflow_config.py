#!/usr/bin/env python3

"""
Script for generating a properly-formatted basic workflow TOML
file suitable for use with cpg-utils convenience functions 
(dataset_path(), remote_tmpdir(), etc.)

Requires CPG_DEPLOY_CONFIG environment variable to be set to an
appropriate Json deployment config, or a pointer to a Json file
with the appropriate deployment config.

python3 generate_workflow_config.py \
    --dataset severalgenomes \
    --access_level test \
    --driver_image image \
    --output_prefix gregsmi \
    --extra_datasets severalgenomes rgp \
    --deploy_config ~/sources/cpg/analysis-runner/deploy-config.prod.json \
    --extra_configs reanalysis_global.toml reanalysis_cohort.toml
"""

from argparse import ArgumentParser
import json

import toml
from cpg_utils.config import read_configs, update_dict
from cpg_utils.deploy_config import DeployConfig, get_deploy_config, get_workflow_config, set_deploy_config
from cpg_utils.storage import get_dataset_bucket_config, get_dataset_bucket_url
from cpg_utils import to_path


parser = ArgumentParser()
parser.add_argument("--dataset", help="name of primary dataset to target", required=True)
parser.add_argument("--access_level", help="whether to target test or main buckets", required=True)
parser.add_argument("--driver_image", help="fully-qualified Hail driver image name to target", required=True)
parser.add_argument("--output_prefix", help="subpath on dataset's access-level bucket to use for output", required=True)
parser.add_argument("--extra_datasets", help="additional datasets to add storage url blocks for", nargs="*")
parser.add_argument("--extra_configs", help="additional workflow-specific toml files to merge in", nargs="*")
parser.add_argument("--deploy_config", help="deployment configuration file to use during generation")
parser.add_argument("--print_only", help="print to screen, don't write to file", action="store_true")
parser.add_argument("-o", help="output file path and name (cloudpathlib.AnyPath-compatible paths are supported)")
args = parser.parse_args()

if args.deploy_config:
    with open(args.deploy_config, "r") as config_file:
        set_deploy_config(DeployConfig.from_dict(json.load(config_file)))

assert args.access_level == "test" or args.access_level == "main"
print("generating workflow toml... using deploy_config:")
print(json.dumps(get_deploy_config().to_dict(), indent=2))

workflow_config = get_workflow_config(args.dataset, args.access_level, args.driver_image, args.output_prefix)

hail_config = {
    "bucket": get_dataset_bucket_url(args.dataset, "hail"),
    "billing_project": args.dataset
}

storage_config = { 
    ds : get_dataset_bucket_config(ds, args.access_level) 
    for ds in [args.dataset] + args.extra_datasets
}

config = {
    "hail": hail_config,
    "workflow": workflow_config,
    "storage": storage_config,
    "CPG_DEPLOY_CONFIG": get_deploy_config().to_dict()
}

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
