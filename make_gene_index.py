#!/usr/bin/env python3

from Bio import SeqIO
import os, glob, sys, time
import threading, traceback, subprocess
from multiprocessing import Pool, cpu_count
from collections import defaultdict
import flash, target_index, build
from multiprocessing import Pool, cpu_count


def goodness_level(target, tags,
                   offtarget_distance_far=target_index.offtag(build.offtarget_proximity["far"]),
                   offtarget_distance_near=target_index.offtag(build.offtarget_proximity["near"])):
    # 3 -- only the best targets (distance to offtarget exceeds 5_9_18)
    # 2 -- close to an offtarget but good structure (distance to offtarget exceeds 5_9_19)
    # 1 -- could even be an offtarget itself but at least has good structure
    # 0 -- anything goes
    # Let D = distance(target, offtarget).  There are 3 cases for D.
    if tags == ["ok"]:
        # target has good structure and D > offtarget_distance_far
        # example target entry "AAAAACTCCGGTGGCCGCCT ok"
        return 3
    if tags == [offtarget_distance_far]:
        # target has good structure and offtarget_distance_near < D <= offtarget_distance_far
        # example target entry "AAAAACTGAATATCAAGGCA off_5_9_18"
        return 2
    if tags == [offtarget_distance_far, offtarget_distance_near]:
        # target has good structure but D <= offtarget_distance_near
        # example target entry "AAAAACTGAGTATCAATGCA off_5_9_18 off_5_9_19"
        return 1
    # target has poor structure, e.g.,
    # "AAAAACTAGAACTAAATATG off_5_9_18 gc_frequency"
    # "AAAAATATATATAAGGCGCT dinucleotide_repeats>3"
    return 0


def index_gene(output_temp_dir, input_file, filtered_targets, ambiguous_targets):
    gene_filename = input_file.split("/")[-1].rsplit(".", 1)[0]
    output_gene_dir = output_temp_dir + "/" + gene_filename
    os.makedirs(output_gene_dir)  # should not exist
    records = list(SeqIO.parse(input_file, "fasta"))
    assert len(records) == 1
    str_desc = str(records[0].description)
    str_seq = str(records[0].seq)
    targets = []
    for i in flash.kmers_range(str_seq, 23):
        if 'G' == str_seq[i+21] == str_seq[i+22]:
            targets.append((i, 'F', str_seq[i:i+20]))
        # Forward and reverse kmer's frequently overlap.
        if 'C' == str_seq[i] == str_seq[i+1]:
            targets.append((i, 'R', flash.reverse_complement(str_seq[i+3:i+23])))
    targets = sorted(targets)
    # Going down the files include fewer of the targets:
    filename = [
        "all.txt",                   # any structure, any distance to offtarget
        "dna_good_structure.txt",    # good structure, any distance to offtarget
        "dna_good_5_9_19.txt",       # good structure and distance to offtarget > 5_9_19
        "dna_good_5_9_18.txt",       # good structure and distance to offtarget > 5_9_18
    ]
    filename = [(output_gene_dir + "/" + fn) for fn in filename]
    with open(filename[0], "w") as f0, \
         open(filename[1], "w") as f1, \
         open(filename[2], "w") as f2, \
         open(filename[3], "w") as f3:
        file = [f0, f1, f2, f3]
        for pos, dir, kmer in targets:
            if kmer in ambiguous_targets:
                continue
            tags = filtered_targets[kmer]
            for j in range(0, goodness_level(kmer, tags) + 1):
                file[j].write("{}\t{}\n".format(pos, dir))


def index_all(output_temp_dir, gene_files, filtered_targets, ambiguous_targets):
    for input_file in gene_files:
        index_gene(output_temp_dir, input_file, filtered_targets, ambiguous_targets)


def rebuild_gene_index():
    t = time.time()
    gene_files = glob.glob('{}/*.fasta'.format(build.genes_dir))
    if not gene_files or not os.path.exists(build.filtered_targets_path):
        print("Could not find input files.")
        return -2
    output_dir = build.gene_index_dir
    output_temp_dir = build.gene_index_temp_dir
    # print("Deleting every file in {}.".format(output_temp_dir))
    subprocess.check_call("rm -rf {}".format(output_temp_dir).split())
    subprocess.check_call("rm -rf {}".format(output_dir).split())    
    os.makedirs(output_temp_dir)
    filtered_targets = target_index.read_tagged_targets(build.filtered_targets_path)
    ambiguous_targets = target_index.read_tagged_targets(build.ambiguous_targets_path)
    # index_all(output_temp_dir, gene_files, filtered_targets, ambiguous_targets)
    num_subprocs = min(cpu_count(), 8)  # there is probably a limit to what the filesystem can handle
    with Pool(num_subprocs) as p:
        size = int((len(gene_files) + num_subprocs - 1) / num_subprocs)
        # Pool.starmap requires Python3.3+
        p.starmap(index_all,
            [(output_temp_dir, gene_files[index : index + size], filtered_targets, ambiguous_targets)
             for index in range(0, len(gene_files), size)])
    print("Moving {} to {}.".format(output_temp_dir, output_dir))
    subprocess.check_call(["/bin/mv", output_temp_dir, output_dir])
    print("Completed rebuild_gene_index in {:3.1f} seconds.".format(time.time() - t))
    return 0


if __name__ == "__main__":
    retcode = rebuild_gene_index()
    sys.exit(retcode)
