#!/usr/bin/env python3

import json
import sys

from common import *


#  Usage:
#
#      python3 display_genes.py < genes.txt
#      echo 'nalC__NC_002516__ARO_3000818' | python3 display_genes.py
#
#  Displays the sequence of a gene. Cyan bars represent cut sites from potential CRISPR targets.
#  Yellow x's represent SNPs.
#  Magenta bars represent CRISPR targets overlapping with SNPs.


def printf(*args, **kwargs):
    print(*args, flush=True, **kwargs)


def main(argv):
    retcode = -999

    printf("Reading gene names from stdin.\n")
    input_gene_names = set(s.strip() for s in sys.stdin.read().split()) - set([""])

    genes = [Gene(name) for name in input_gene_names]

    with open('inputs/additional/padding.json', 'r') as fp:
        padding_seqs = json.load(fp)

    for g in genes:
        g.load_targets("dna_good_5_9_18.txt")
        if g.name in padding_seqs:
            g.verify_padding(padding_seqs[g.name])
        else:
            assert g.padding == None

    for g in genes:
        printf(g.name)
        g.display_gene_targets()
        printf("")

        retcode = 1

    return retcode


if __name__ == "__main__":
    retcode = main(sys.argv)
    sys.exit(retcode)
