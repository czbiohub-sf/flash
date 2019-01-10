#!/usr/bin/env python3

import subprocess
import sys
import traceback
import time
import requests


# For given radius c5_c10_c20 around a target we can determine if there are
# any offtargets in that radius.  There are two radii of interest.
#
offtarget_proximity = {
    "far":  "5_9_18",
    "near":  "5_9_19"
}
assert offtarget_proximity["far"] < offtarget_proximity["near"]
#
# Def: A radius c5_c10_c20 specifies maximum Hamming distances on nested
# suffixes of 5, 10, 20 bases.  A 20-mer Y is within 5_9_18 radius of 20-mer X
# if and only if the 5-char suffixes of X and Y match exactly, the 10-char
# suffixes have at most 1 positional difference, and the 20-char suffixes
# (i.e. the entire X and Y) have at most 2 positional differences.

gvc_top = "generated_files/under_version_control"

# Contains .txt files with various indexes of target => target properties
target_index_dir = "{}/target_index".format(gvc_top)

# Produced by make_genes_and_identify_all_targets.py, formerly split_fasta
genes_dir = "{}/genes".format(gvc_top)
genes_temp_dir = genes_dir + "_temp"
all_targets_path = "{}/all_targets.txt".format(target_index_dir)
ambiguous_targets_path = "{}/ambiguous_targets.txt".format(target_index_dir)
antibiotics_by_gene_path = "{}/antibiotics_by_gene.txt".format(target_index_dir)
genes_by_antibiotic_path = "{}/genes_by_antibiotic.txt".format(target_index_dir)
antibiotics_path = "{}/antibiotics.txt".format(target_index_dir)

# Produced by filter_offtargets.py with input from all_targets_path
# and using the GO server for offtarget filtering.
off_targets_path = "{}/off_targets.txt".format(target_index_dir)

# Produced by filter_targets.py with input from all_targets_path
# and off_targets_path.
filtered_targets_path = "{}/filtered_targets.txt".format(target_index_dir)

# Produced by make_gene_index, formerly known as crispr_sites.py.
gene_index_dir = "{}/gene_index".format(gvc_top)
gene_index_temp_dir = "{}/gene_index_temp".format(gvc_top)

padding_input_path = "inputs/additional/padding.yaml"

special_token = "rebuild_version_controlled_genes"


# It is convenient to track changes with git for the smaller and more human
# readable generated files.
#
# When files are (re)generated, the build process presents a git status
# report showing which, if any, of the generated files have changed.
#
# The changes are automatically added to the git index, i.e., to the
# user's pending commit.


def git_reset_and_remove_generated_folder(output_dir, return_values={}):
    try:
        assert ' ' not in output_dir
        # t = time.time()
        print("Deleting every file in {}.".format(output_dir))
        subprocess.check_call("rm -rf {}".format(output_dir).split())
        subprocess.check_call("git rm -r --cached --force --ignore-unmatch --quiet {}".format(output_dir).split())
        # print("GIT houskeeping took {:3.1f} seconds.".format(time.time() - t))
        return_values['status'] = 'Success'
    except:
        traceback.print_exc()
        return_values['status'] = 'Failure'


def git_add_back_generated_folder(output_dir, return_values={}):
    try:
        assert output_dir
        assert ' ' not in output_dir
        print("Adding every file in {output_dir} to git.".format(output_dir=output_dir))
        subprocess.check_call("git add -A {}".format(output_dir).split())
        subprocess.check_call("git status --short {}".format(output_dir).split())
        return_values['status'] = 'Success'
    except:
        traceback.print_exc()
        return_values['status'] = 'Failure'


def git_remove_generated_file(output_file):
    assert ' ' not in output_file
    subprocess.check_call("rm -f {}".format(output_file).split())
    subprocess.check_call("git rm --cached --force --ignore-unmatch --quiet {}".format(output_file).split())


def git_add_back_generated_file(output_file, skip_status=False):
    assert output_file
    assert ' ' not in output_file
    print("Adding {} back to git.".format(output_file))
    subprocess.check_call("git add -A {}".format(output_file).split())
    if not skip_status:
        subprocess.check_call("git status --short {}".format(output_file).split())


def fetch_with_retries(targets, c5, c10, c20, timeout=300, max_response_time=600):
    sleep_time = 2
    failures = 0
    while True:
        try:
            url = "http://127.0.0.1:8080/search?targets=%s&limits=%s"
            url = url % (",".join(map(str, targets)), ",".join(map(str, [c5, c10, c20])))
            return requests.get(url, timeout=max_response_time)
        except requests.exceptions.ConnectionError:
            failures += 1
            if (failures * sleep_time) > timeout:
                raise
            time.sleep(sleep_time)
