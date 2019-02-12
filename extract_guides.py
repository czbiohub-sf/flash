import argparse

import build
import common
from padding import get_gene_to_padding


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--library",
                        metavar="file",
                        type=argparse.FileType("r"),
                        required=True)
    parser.add_argument("--genes",
                        metavar="file",
                        type=argparse.FileType("r"),
                        required=True)
    parser.add_argument("--padding",
                        metavar="file",
                        type=argparse.FileType("r"),
                        default=build.padding_input_path)
    parser.add_argument("--max-cuts",
                        help="Max number of cuts to keep per gene.",
                        type=int)
    parser.add_argument("--output",
                        help="Path to the output file.",
                        type=argparse.FileType("w"),
                        required=True)

    return parser.parse_args()


def main():
    def read_file(f): return [e.strip() for e in f.readlines()]
    args = parse_args()
    gene_names = read_file(args.genes)
    library = read_file(args.library)
    gene_to_padding = get_gene_to_padding(args.padding.name)

    # Set up Gene objects.
    genes = [
        common.Gene(gname, gene_to_padding.get(gname)) for gname in gene_names]
    [g.load_targets("dna_good_5_9_18.txt") for g in genes]
    [g.cut_with_library(library) for g in genes]

    # select all guides that are part of the cuts
    guides = set()
    for gene in genes:
        if args.max_cuts:
            gene_guides = [
                e for g in genes for e in g.trim_library(args.max_cuts)]
        else:
            gene_guides = [t.guide for t in gene.cuts]

        for guide in gene_guides:
            if guide in library:
                guides.add(guide)

    for guide in sorted(guides):
        args.output.write("{}\n".format(guide))


if __name__ == "__main__":
    main()
