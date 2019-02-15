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
from collections import defaultdict
import json
from Bio import SeqIO

import flash
import build
import aro
from common import MutationIndex
from padding import get_gene_to_padding


# It splits large fasta files into individual gene files, naming each output
# fasta file according to its gene's locus, resistance, and origin DB.
# Checks that genes do not accidentally clobber each other's files.
# Dedups genes that occur in multiple source files based on DNA sequences alone.
#
# TODO: There are a bunch of heuristics to come up with reasonable gene names
# for CARD.  These probably should be superceded by logic based on parsing
# the inputs/card/aro.obo file.  Parsing that file will be necessary anyway
# to correctly trace and generate antibiotic resistance tags.
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
# A gene is identified by its LOCUS, ORIGIN and ARO/resfinder filename.
# It is common for two gene keys to map to the same DNA, for example
#
#    vatA__L07778__ARO_3002840
#    vatA__L07778__macrolide
#


def long_word(s):
    return s.isalpha() and len(s) > 6


CARD_ORIGIN_ABBREVIATONS = {
    'nucleotide_fasta_protein_variant_model': 'protein_variant',
    'nucleotide_fasta_protein_homolog_model': 'protein_homolog',
    'nucleotide_fasta_protein_overexpression_model': 'protein_overexpression'
}


def extract_filename(path):
    return os.path.splitext(os.path.basename(path))[0]


def first_bacterial_gene_name(text):
    # Return the first token in text that looks like a bacterial gene name,
    # or at least like a protein, or anything biological...  Hacky but works
    # for 90% of the card protein_variant_model descriptions.
    for word in text.split():
        # See http://www.biosciencewriters.com/Guidelines-for-Formatting-Gene-and-Protein-Names.aspx
        if len(word) == 4 and word.isalnum() and (word[:3].lower() + word[-1:].upper()) == word:
            return word
    for word in text.split()[:5]:
        if len(word) == 3 and word.lower() == word and word not in ["and", "for"]:
            return word
    for word in text.split():
        if 2 < len(word) < 6 and (not word.isalpha() or word.lower() != word):
            return word
    return None


class GeneRecord(object):

    def __init__(self, str_seq, str_desc, input_path):
        self.name = None
        self.locus = None
        self.description = None
        self.aro = None
        self.inferred_aro = None
        self.resistances = None
        self.disambiguation = None
        self.key = None
        self.origin = None
        # str_seq and str_desc are from fasta
        self.str_seq = str_seq
        self.str_desc = str_desc
        if 'inputs/card' in input_path:
            self.parse_card(str_desc, input_path)
        elif 'inputs/resfinder' in input_path:
            self.parse_resfinder(str_desc, input_path)
        elif 'inputs/additional' in input_path:
            self.parse_additional(str_desc)
        else:
            self.parse_generic(str_desc)

    def parse_card(self, descr, input_path):
        self.origin = 'card'
        self.origin_short_filename = CARD_ORIGIN_ABBREVIATONS[extract_filename(input_path)]
        dl = descr.split("|")
        if (dl[0] != 'gb') or (len(dl) < 3):
            raise RuntimeError("Unexpected record description format.")
        self.locus = dl[1]
        if len(dl) > 4 and dl[4].startswith("ARO"):
            self.aro = dl[4][4:]
        self.name, self.description = dl[-1].split(' ', 1)
        if long_word(self.name):
            # heuristic fail, this "name" is not really a gene name
            # very typical for card protein_variant_models
            self.name = first_bacterial_gene_name(dl[-1])
            if self.name == None:
                # second heuristic fail
                # if you see this, please improve the heuristic
                self.name = "NoName" + self.aro
            self.description = dl[-1]
        self.resistances = []

    def parse_resfinder(self, descr, input_path):
        self.origin='resfinder'
        dl = descr.split("_", 2)
        self.name, self.locus = dl[0], dl[2]
        self.origin_short_filename = extract_filename(input_path).replace('-', '_')
        self.description = None
        # clean up so the resfinder genes appear next to the identical card genes
        # when sorted in lexicographic order
        if self.origin_short_filename == "beta_lactamase" and self.name.startswith("bla"):
            self.name = self.name[3:]
        # so amazingly clean in resfinder, though it doesn't address multiple resistances
        # and ignores the notes.txt file
        self.resistances = [self.origin_short_filename + "_resfinder"]
        self.disambiguation = dl[1]

    def parse_additional(self, descr):
        self.origin='additional'
        dl = descr.split("|")
        self.resistances = []
        raw_resistance = ""
        for part in dl:
            if part.startswith("flash_key:"):
                self.key = part.split("flash_key:", 1)[1]
            if part.startswith("flash_resistance:"):
                if raw_resistance:
                    raise RuntimeError("Multiple flash_resistance: fields per description.  Use comma separated values instead.")
                raw_resistance = part
                self.resistances.extend(part.split(":", 1)[1].split(","))
        if sanitize(self.key) != self.key:
            raise RuntimeError("flash_key '{}' contains disallowed characters".format(self.key))
        key_parts = self.key.split("__")
        if len(key_parts) != 3:
            raise RuntimeError("flash_key '{}' does not have 3 __-separated parts".format(self.key))
        if not key_parts[2].endswith("_additional") or key_parts[2].endswith("__additional"):
            raise RuntimeError("flash_key '{}' must end with <something>_additional".format(self.key))
        self.origin_short_filename = key_parts[2]
        self.name, self.locus = key_parts[0], key_parts[1]
        self.description = None
        self.resistances = sorted(r.strip() + "_additional" for r in set(self.resistances))
        # May relax this at some point
        if not any(self.origin_short_filename.startswith(r) for r in self.resistances):
            raise RuntimeError("Unconventional flash_resistance '{}' for flash_key '{}'".format(raw_resistance, self.key))

    def parse_generic(self, descr):
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


valid_chars = "-_.%s%s" % (string.ascii_letters, string.digits)

def sanitize(filename):
    f = filename.replace(' ', '_').replace("'", "p").replace("[", "__").replace("]", "__").replace("___", "__").replace(
        "___", "__").replace("___", "__").replace("/", "__").replace(":","__")
    if f.endswith("__"):
        f = f[:-2]
    return ''.join(c for c in f if c in valid_chars)

def shorten(filename):
    return "__".join(part for part in filename.split("__")[:4] if not part.startswith("from_card_")).replace("from_resfinder_", "")

def longer_filename(gene, use_disambiguation=False):
    disambiguated_locus = str(gene.locus)
    if gene.disambiguation != None and use_disambiguation:
        disambiguated_locus += "_{}".format(gene.disambiguation)
    result = "{gene_name}__{gene_locus}__from_{gene_origin}_{origin_short_filename}".format(
        gene_name=gene.name,
        gene_locus=disambiguated_locus,
        gene_origin=gene.origin,
        origin_short_filename=gene.origin_short_filename
    )
    if gene.aro != None:
        result += "__ARO_{}".format(gene.aro)
    if gene.description != None:
        result += "__{}".format(gene.description)
    return sanitize(result)

def filename(gene, use_disambiguation=False):
    if gene.key != None:
        return gene.key
    return shorten(longer_filename(gene, use_disambiguation))

def read_file(sequences, aro_genes, input):
    for record in SeqIO.parse(input, "fasta"):
        str_desc = str(record.description)
        str_seq = str(record.seq).upper()

        # Replace all ambiguous bases with N
        for c in ['U', 'Y', 'R', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V']:
            str_seq = str_seq.replace(c, 'N')

        # Remove '-', which indicates a deletion in an alignment
        str_seq = str_seq.replace('-', '')
        gr = GeneRecord(str_seq, str_desc, input)
        if gr.aro:
            aro_genes["ARO:" + gr.aro] = True
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


def output_unique_sequence(antibiotics_for_gene, genes_for_antibiotic, all_files_lc, all_targets, str_seq, gene_records, padding_seq, output_dir, mutation_index):
    canonical_gene = gene_records[0]
    canonical_key = filename(gene_records[0])
    inferred_aro = None
    inferred_resistances = set()
    aliases = []
    for gr in gene_records:
        gr.ofn = filename(gr)
        aliases.append(gr.ofn)
        if gr.resistances:
            inferred_resistances |= set(gr.resistances)
        if gr.aro:
            inferred_aro = gr.aro
            # We prefer the canonical key to include the ARO
            canonical_gene = gr
    inferred_resistances = sorted(inferred_resistances, key=str.lower)
    str_desc = "|".join(gr.str_desc.strip(">|")
        for gr in sorted(gene_records, key=lambda g: g.ofn.lower()))
    str_seq = canonical_gene.str_seq
    canonical_key = canonical_gene.ofn
    if canonical_key.lower() in all_files_lc:
        canonical_key = filename(gr, use_disambiguation=True)
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
        "flash_aro",
        "flash_aliases",
        "flash_mutation_ranges")
    if not inferred_resistances:
        str_desc += "|flash_resistance:flash_todo"
    else:
        str_desc += "|flash_resistance:{}".format(",".join(str(r) for r in inferred_resistances))
        antibiotics_for_gene[canonical_key] = inferred_resistances
        for ir in inferred_resistances:
            genes_for_antibiotic[ir].append(canonical_key)
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
    if inferred_aro != None:
        str_desc += "|flash_aro:{}".format(inferred_aro)
        multi_values["flash_mutation_ranges"].extend(mutation_index.mutations.get(str(inferred_aro), []))
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
    # no gene is this long
    infinity = 1000000000
    for i in flash.kmers_range(str_seq, 23):
        if 'G' == str_seq[i+21] == str_seq[i+22]:
            d = all_targets[str_seq[i:i+20]]
            d[canonical_key] = min(flash.cut_location((i, 'F')), d.get(canonical_key, infinity))
        if 'C' == str_seq[i] == str_seq[i+1]:
            d = all_targets[flash.reverse_complement(str_seq[i+3:i+23])]
            d[canonical_key] = min(flash.cut_location((i, 'R')), d.get(canonical_key, infinity))

def output_targets(output_path, targets_dict):
    with open(output_path, "w") as output_handle:
        sorted_targets = sorted(targets_dict.keys())
        for target in sorted_targets:
            output_handle.write(target + "\n    " +
                "\n    ".join((canonical_key + " " + str(targets_dict[target][canonical_key])) for canonical_key in sorted(targets_dict[target].keys(), key=str.lower)) + "\n\n")

def split_all(input_files, output_dir, all_targets_index_path, ambiguous_targets_index_path, padding_input_path, antibiotics_by_gene_path, genes_by_antibiotic_path, antibiotics_path):
    gene_to_padding = get_gene_to_padding(padding_input_path)
    mutation_index = MutationIndex()
    sequences = defaultdict(list)
    aro_genes = {}
    for input in input_files:
        read_file(sequences, aro_genes, input)
    all_terms = aro.parse_aro_obo("inputs/card/aro.obo")
    aro_antibiotics_for_gene = aro.infer_resistances(all_terms, aro_genes)
    for gene_records in sequences.values():
        for gr in gene_records:
            if gr.aro:
                gene_aro = "ARO:" + gr.aro
                assert gene_aro in aro_genes
                if gene_aro not in aro_antibiotics_for_gene:
                    print("WARNING:  Gene {} not found in ARO.OBO".format(gene_aro))
                else:
                    gr.resistances.extend(a[0] for a in aro_antibiotics_for_gene[gene_aro])
    all_targets = defaultdict(dict)
    all_files_lc = set()
    antibiotics_for_gene = {}
    genes_for_antibiotic = defaultdict(list)
    for str_seq, gene_records in sorted(sequences.items()):
        output_unique_sequence(antibiotics_for_gene, genes_for_antibiotic, all_files_lc, all_targets, str_seq, gene_records, gene_to_padding, output_dir, mutation_index)
    for ab in genes_for_antibiotic:
        genes_for_antibiotic[ab] = sorted(genes_for_antibiotic[ab], key=str.lower)
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
    output_genes_by_antibiotic(genes_by_antibiotic_path, genes_for_antibiotic)
    output_antibiotics(antibiotics_path, genes_for_antibiotic)
    output_antibiotics_by_gene(antibiotics_by_gene_path, antibiotics_for_gene)
    return 0

def output_antibiotics_by_gene(antibiotics_by_gene_path, antibiotics_for_gene):
    with open(antibiotics_by_gene_path, "w") as output_handle:
        for gene_key, inferred_resistances in sorted(antibiotics_for_gene.items(), key=lambda p: p[0].lower()):
            output_handle.write(
                gene_key +
                "\n    " +
                "\n    ".join(inferred_resistances) +
                "\n\n"
            )

def output_genes_by_antibiotic(genes_by_antibiotic_path, genes_for_antibiotic):
    with open(genes_by_antibiotic_path, "w") as output_handle:
        for antibiotic_fqn in sorted(genes_for_antibiotic.keys(), key=str.lower):
            output_handle.write(
                str(antibiotic_fqn).ljust(100) +
                str(len(genes_for_antibiotic[antibiotic_fqn])).rjust(20) +
                "\n    " +
                "\n    ".join(genes_for_antibiotic[antibiotic_fqn]) +
                "\n\n"
            )

def output_antibiotics(antibiotics_path, genes_for_antibiotic):
    with open(antibiotics_path, "w") as output_handle:
        for antibiotic_fqn in sorted(genes_for_antibiotic.keys(), key=str.lower):
            output_handle.write(
                str(antibiotic_fqn).ljust(100) +
                str(len(genes_for_antibiotic[antibiotic_fqn])).rjust(20) +
                "\n\n"
            )

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
                        metavar="file",
                        default=build.padding_input_path)
    parser.add_argument("--disable-git",
                        help="Do not add changed files to git.",
                        action='store_true')

    return parser.parse_args()

def get_files(args_obj):
    if args_obj.targets_dir:
        input_files = glob.glob(args.targets_dir + '/*.fasta')
    elif args_obj.targets:
        input_files = [args.targets.name]
    else:
        input_files = (
        glob.glob('inputs/card/*.fasta') +
        glob.glob('inputs/resfinder/*.fsa') +
        glob.glob('inputs/additional/*.fasta')
            )
    if not input_files:
        print("Could not find input files.")
        sys.exit(-2)
    return input_files

def make_genes_and_identify_all_targets(files,padding=None):
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
        padding,
        build.antibiotics_by_gene_path,
        build.genes_by_antibiotic_path,
        build.antibiotics_path
    )

    print("Moving {} to {}.".format(output_temp_dir, output_dir))
    subprocess.check_call(["/bin/mv", output_temp_dir, output_dir])
    print("Completed make_genes_and_identify_all_targets in {:3.1f} seconds".format(time.time() - t))
    return 0


if __name__ == "__main__":
    args = parse_args()
    input_files= get_files(args)

    retcode = make_genes_and_identify_all_targets(padding=args.padding.name,files=input_files)
    sys.exit(retcode)
