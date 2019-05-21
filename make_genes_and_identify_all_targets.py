#!/usr/bin/env python3

import os
import os.path
import string
import glob
import subprocess
import sys
import time
import threading
import argparse
import math
from collections import defaultdict
import json
from Bio import SeqIO

import build
from common import MutationIndex
import flash
from padding import (
    get_gene_to_padding,
    set_to_none_if_padding_not_provided
)


# It splits large fasta files into individual gene files, naming each output
# fasta file according to its gene's locus and resistance.
# Checks that genes do not accidentally clobber each other's files.
# Dedups genes that occur in multiple source files based on DNA sequences alone.
#
# Output the (target, gene) incidence table in genes/all_targets.txt, like so:
#
# CCCCGACAATCTGCGGGCGC
#     embB__AE000516__ARO_3003465
#     embB__AL123456__ARO_3003326
#     embB__BX248333__ARO_3003325
#
# CCCCGACATCCCCCGGACCA
#     sul2__AF497970__sulphonamide
#     sul2__AM183225__sulphonamide
#     sul2__AY055428__ARO_3000412
#     sul2__AY232670__sulphonamide
#     sul2__AY524415__sulphonamide
#     sul2__FJ197818__sulphonamide
#
# ....
#
# The gene keys are listed only to help human reviewers;  downstream code
# only cares about the 20mers.
#
# A gene is identified by its LOCUS and ORIGIN.
# It is common for two gene keys to map to the same DNA, for example
#
#    vatA__L07778__ARO_3002840
#    vatA__L07778__macrolide
#


class GeneRecord(object):

    def __init__(self, str_seq, str_desc, input_path):
        self.name = None
        self.locus = None
        self.resistances = None
        self.disambiguation = None
        self.key = None
        self.origin = None

        # str_seq and str_desc are from fasta
        self.str_seq = str_seq
        self.str_desc = str_desc

        self.parse(str_desc)

    def parse(self, descr):
        self.origin = 'generic'
        dl = descr.split("|")
        self.key = None
        for part in dl:
            if part.startswith("flash_key"):
                self.key = part.split(":", 1)[1]
        if self.key is None:
            self.key = sanitize(dl[0])
        if sanitize(self.key) != self.key:
            raise RuntimeError("flash_key '{}' contains disallowed characters".format(self.key))

    def filename(self):
        return self.key


VALID_CHARS = "-_.%s%s" % (string.ascii_letters, string.digits)


def sanitize(filename):
    f = filename.replace(' ', '_') \
                .replace("'", "p") \
                .replace("[", "__") \
                .replace("]", "__") \
                .replace("___", "__") \
                .replace("___", "__") \
                .replace("___", "__") \
                .replace("/", "__") \
                .replace(":","__")
    if f.endswith("__"):
        f = f[:-2]

    return ''.join(c for c in f if c in VALID_CHARS)


def read_file(sequences, input):
    for record in SeqIO.parse(input, "fasta"):
        str_desc = str(record.description)
        str_seq = str(record.seq).upper()

        # Replace all ambiguous bases with N
        for c in ['U', 'Y', 'R', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V']:
            str_seq = str_seq.replace(c, 'N')

        # Remove '-', which indicates a deletion in an alignment
        str_seq = str_seq.replace('-', '')
        gr = GeneRecord(str_seq, str_desc, input)
        sequences[str_seq].append(gr)


def remove_tags(str_desc, *tags_list):
    """Pull out the values from the (possibly multiple appearances of) each tag,
    and remove from str_desc."""
    separator = "|"
    colon = ":"
    multi_values = defaultdict(list)
    for tag in tags_list:
        if not tag.endswith(colon):
            tag += colon
        while tag in str_desc:
            parts = str_desc.split(tag, 1)
            multi_values[tag[:-1]].append(parts[1].split(separator, 1)[0])
            if separator in parts[1]:
                parts[1] = parts[1].split(separator, 1)[-1]
            else:
                if parts[0].endswith(separator):
                    parts[0] = parts[0][:-1]
                parts[1] = ""
            str_desc = parts[0] + parts[1]
    return str_desc, multi_values


def output_unique_sequence(all_files_lc, all_targets, str_seq, gene_records, padding_seq, output_dir, mutation_index):
    canonical_gene = gene_records[0]
    canonical_key = gene_records[0].filename()
    inferred_resistances = set()
    aliases = []
    for gr in gene_records:
        gr.ofn = gr.filename()
        aliases.append(gr.ofn)
        if gr.resistances:
            inferred_resistances |= set(gr.resistances)
    inferred_resistances = sorted(inferred_resistances, key=str.lower)
    str_desc = "|".join(gr.str_desc.strip(">|")
        for gr in sorted(gene_records, key=lambda g: g.ofn.lower()))
    str_seq = canonical_gene.str_seq
    canonical_key = canonical_gene.ofn
    if canonical_key.lower() in all_files_lc:
        canonical_key = gr.filename()
    if canonical_key.lower() in all_files_lc:
        print("Two different DNA sequences want filename '{}'. Please revise the heuristics.".format(canonical_key))
        raise RuntimeError("Gene filenames clash.")
    # we add it after conversion to lowercase because MacOS X filesystem
    # is case-insensitive (but case-preserving, tricky that way)
    all_files_lc.add(canonical_key.lower())
    out_file = output_dir + "/" + canonical_key + ".fasta"
    # Remove pre-existing flash tags, we will re-add some back here with merged values.
    str_desc, multi_values = remove_tags(str_desc,
        "flash_resistance",
        "flash_key",
        "flash_padding",
        "flash_dna_key",
        "flash_aliases",
        "flash_mutation_ranges")
    if not inferred_resistances:
        str_desc += "|flash_resistance:flash_todo"
    else:
        str_desc += "|flash_resistance:{}".format(",".join(str(r) for r in inferred_resistances))
    str_desc +="|flash_key:{}".format(canonical_key)
    if padding_seq.get(canonical_key):
        # common case
        padding = padding_seq[canonical_key]
        del padding_seq[canonical_key]
    else:
        # probably will never happen
        padding = None
        for gk in aliases:
            if padding_seq.get(gk):
                padding = padding_seq[gk]
                del padding_seq[gk]
                break
    if padding != None:
        prefix, suffix = padding.prefix, padding.suffix
        str_seq = prefix + str_seq + suffix
        str_desc += "|flash_padding:{}_{}".format(len(prefix), len(suffix))
    if multi_values.get("flash_mutation_ranges"):
        str_desc += "|flash_mutation_ranges:{}".format(
            ",".join(multi_values["flash_mutation_ranges"]))
    if len(gene_records) > 1:
        str_desc += "|flash_aliases:{}".format(",".join(sorted(aliases, key=str.lower)))
    with open(out_file, "w") as output_handle:
        # SeqIO.write(record, output_handle, "fasta")
        output_handle.write(">" + str_desc + "\n")
        line_length=60
        output_handle.write("\n".join(str(line) for line in (str_seq[n:n+line_length] for n in range(0, len(str_seq), line_length))) + "\n")
    for i in flash.kmers_range(str_seq, 23):
        if 'G' == str_seq[i+21] == str_seq[i+22]:
            d = all_targets[str_seq[i:i+20]]
            d[canonical_key] = min(flash.cut_location((i, 'F')), d.get(canonical_key, math.inf))
        if 'C' == str_seq[i] == str_seq[i+1]:
            d = all_targets[flash.reverse_complement(str_seq[i+3:i+23])]
            d[canonical_key] = min(flash.cut_location((i, 'R')), d.get(canonical_key, math.inf))


def output_targets(output_path, targets_dict):
    with open(output_path, "w") as output_handle:
        sorted_targets = sorted(targets_dict.keys())
        for target in sorted_targets:
            output_handle.write(target + "\n    " +
                "\n    ".join((canonical_key + " " + str(targets_dict[target][canonical_key])) for canonical_key in sorted(targets_dict[target].keys(), key=str.lower)) + "\n\n")


def split_all(input_files, output_dir, all_targets_index_path, ambiguous_targets_index_path, padding_input_path):
    gene_to_padding = get_gene_to_padding(padding_input_path)
    mutation_index = MutationIndex()
    sequences = defaultdict(list)
    for input in input_files:
        read_file(sequences, input)
    all_targets = defaultdict(dict)
    all_files_lc = set()
    for str_seq, gene_records in sorted(sequences.items()):
        output_unique_sequence(all_files_lc, all_targets, str_seq, gene_records, gene_to_padding, output_dir, mutation_index)
    if gene_to_padding:
        print()
        print("WARNING:  UNUSED PADDING SEQUENCES:  CHECK KEYS, FILENAMES MAY HAVE CHANGED:")
        print()
        print(json.dumps(gene_to_padding))
    ambiguous_targets = {}
    ambiguous_genes = set()
    for target, gene_pos in all_targets.items():
        if target.strip("ACGT"):
            genes = sorted(gene_pos.keys(), key=str.lower)
            ambiguous_targets[target] = gene_pos
            ambiguous_genes |= set(genes)
    for target in ambiguous_targets:
        del all_targets[target]
    if ambiguous_targets:
        print("Will discard {} targets in amibiguous_targets.txt affecting {} not necessarily unique genes."
              .format(len(ambiguous_targets), len(ambiguous_genes)))
    else:
        print("There are no ambiguous targets, i.e., no target contains N.")
    print("Will output {} unique targets over {} not necessarily unique genes."
          .format(len(all_targets), len(all_files_lc)))
    l = max(len(v) for k,v in all_targets.items())
    for k, v in all_targets.items():
        if len(v) == l:
            break
    print("The most popular target {} is hit by {} genes.".format(k, l))
    output_targets(all_targets_index_path, all_targets)
    output_targets(ambiguous_targets_index_path, ambiguous_targets)
    return 0


def parse_args():
    parser = argparse.ArgumentParser()
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument("--targets-dir",
                             help="Directory containing input gene fastas.",
                             type=str)
    input_group.add_argument("--targets",
                             help="Fasta file containing target genes.",
                             type=argparse.FileType("r"),
                             metavar="file")
    parser.add_argument("--padding",
                        help="yaml or json file with padding info.",
                        type=argparse.FileType("r"),
                        metavar="file"
                        )
    return parser.parse_args()


def make_genes_and_identify_all_targets(files, padding=None):
    t = time.time()

    output_dir = build.genes_dir
    output_temp_dir = build.genes_temp_dir

    subprocess.check_call("rm -rf {}".format(output_temp_dir).split())
    subprocess.check_call("rm -rf {}".format(output_dir).split())

    os.makedirs(output_temp_dir)

    split_all(
        files,
        output_temp_dir,
        build.all_targets_path,
        build.ambiguous_targets_path,
        padding
    )

    print("Moving {} to {}.".format(output_temp_dir, output_dir))
    subprocess.check_call(["/bin/mv", output_temp_dir, output_dir])
    print("Completed make_genes_and_identify_all_targets in {:3.1f} seconds".format(time.time() - t))
    return 0


if __name__ == "__main__":
    args = parse_args()
    input_files= [args.targets.name]
    padding_file = set_to_none_if_padding_not_provided(args.padding)
    retcode = make_genes_and_identify_all_targets(padding=padding_file, files=input_files)
    sys.exit(retcode)
