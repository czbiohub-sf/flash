import argparse
import sys
import json
from common import *


def parse_args():
    parser = argparse.ArgumentParser("Tool for visualizing CRISPR targets and"
    "cut sites in a gene. Sample command-line:\n\n"
    "python display_genes.py nalC__NC_002516__ARO_3000818 mecA__AB033763__beta_lactamase --library library.txt")
    parser.add_argument("genes",
                        help="flash_key of gene to be displayed.",
                        type=str,
                        nargs='+',)

    parser.add_argument("--library",
                        metavar="file",
                        type=argparse.FileType("r"),
                        required=False)

    return parser.parse_args()


def main():
    args = parse_args()

    library = None
    if args.library:
        def read_file(f): return [e.strip() for e in f.readlines()]
        library = read_file(args.library)

    for gene_key in args.genes:
        g = Gene(gene_key)

        if g.padding:
            with open('inputs/additional/padding.json', 'r') as fp:
                padding_seqs = json.load(fp)
            if g.name in padding_seqs:
                g.verify_padding(padding_seqs[g.name])

        g.load_targets("dna_good_5_9_18.txt")

        if library:
            g.cut_with_library(library)

        print("--------------------------------------")
        print("Displaying CRISPR targets for " + g.name + ":\n")

        print("Cut sites of potential guides are denoted by pipes (|). "
        "Mutation ranges are denoted by a yellow x. "
        "Color of | indicates state of the corresponding guide: pink if it "
        "would cover a mutation, red if it is in the library, and cyan otherwise.")
        g.display_gene_targets()
        print("")

        if library:
            print("Cut sites from the library are in red. "
            "Bases contained in expected paired-end 150 bp reads are green.")
            g.display_gene_cuts()
            print("")

if __name__ == "__main__":
    retcode = main()
    sys.exit(retcode)
