#!/usr/bin/env python3

import subprocess
import sys
import traceback
import time
import requests

import build


def main():
    t = time.time()
    print("-------------------------------------------------------------------------------------")
    subprocess.check_call("python filter_offtarget.py".split())
    print("-------------------------------------------------------------------------------------")
    subprocess.check_call("python filter_targets.py".split())
    print("-------------------------------------------------------------------------------------")
    subprocess.check_call("python make_gene_index.py".split())
    print("-------------------------------------------------------------------------------------")
    print("Complete rebuild took {:3.1f} seconds.".format(time.time() - t))
    return 0


if __name__ == "__main__":
    print("Builder of FLASH.  For usage, see README.")
    try:
        print("Poking offtarget server.  Timeout 300 seconds.")
        build.fetch_with_retries(["ACGT" * 5], 5, 9, 18, timeout=300)
        print("Offtarget server is alive.")
    except:
        traceback.print_exc()
        print("*********************************************************************************")
        print("***   Did you forget to start the offtarget server?  Did it finish loading?   ***")
        print("***   Please follow the instructions in README.TXT and try again.             ***")
        print("*********************************************************************************")
        sys.exit(-1)
    retcode = main()
    sys.exit(retcode)
