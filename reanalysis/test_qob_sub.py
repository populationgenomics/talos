#import os
import logging
import sys
import cpg_utils.hail_batch
#import hail
#import asyncio

#from hail.utils.java import Env

def main():

    cpg_utils.hail_batch.init_batch(
        jar_url="hail-az://hailms02batch/query/jars/1078abac8b8e1c14fe7743aa58bc25118b4108de.jar",
        driver_cores=8, 
        driver_memory='highmem')

    # if Env._hc:
    #     print("Hail Context already initialized")
    # else:
    #     logging.info("Initializing QoB")

    #     batch = asyncio.get_event_loop().run_until_complete(
    #         hail.init_batch(
    #             billing_project="severalgenomes", 
    #             remote_tmpdir="hail-az://sevgen002sa/cpg-severalgenomes-hail",
    #             jar_url="hail-az://hailms02batch/query/jars/1078abac8b8e1c14fe7743aa58bc25118b4108de.jar",
    #             driver_memory="highmem",
    #             driver_cores=8,
    #             default_reference="GRCh38"
    #         )
    #     )




if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    main()
