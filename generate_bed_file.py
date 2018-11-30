import argparse
import re
import sys

from Bio import SeqIO

sys.path.append('..')
import flash


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fasta", type=argparse.FileType("r"), required=True)
    parser.add_argument("--library", type=argparse.FileType("r"), required=True)
    parser.add_argument("--output", type=argparse.FileType("w"), required=True)
    return parser.parse_args()


def main(input_fasta, library, output):
    gene_names_to_seqs = {}
    for fasta in SeqIO.parse(input_fasta, format="fasta"):
        gene_names_to_seqs[fasta.id] = str(fasta.seq)
    library = {e.strip() for e in open(library).readlines()}

    with open(output, "w") as fh:    
        for gene_name, gene_seq in gene_names_to_seqs.items():
            for i in range(len(gene_seq) - 23):
                start = i
                end = i+23        
                kmer = gene_seq[start:end]
                if kmer.endswith('GG'):
                    guide = kmer[:20]
                    if guide in library:
                        fh.write("{}\t{}\t{}\t{}\n".format(gene_name,
                                                           start, end,
                                                           guide))   
                if kmer.startswith('CC'):
                    rc_guide = flash.reverse_complement(kmer)[:20]
                    if rc_guide in library:
                        fh.write("{}\t{}\t{}\t{}\n".format(gene_name,
                                                           start, end,
                                                           rc_guide))                   


if __name__ == '__main__':
    args = parse_args()
    main(
        args.input_fasta.name,
        args.library.name,
        args.output.name
    )
