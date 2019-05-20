import csv
import re
import os
import sys
from collections import namedtuple

from Bio import SeqIO
from colors import color

import build
import flash


READ_LENGTH = 150
IDEAL_CUTOFF = 200
OKAY_CUTOFF = 301
LONG_CUTOFF = 501


Target = namedtuple('Target', ['guide', 'cut'])


class Component(object):
    """A collection of genes for which a single library will be built.
    The library can be subsetted to the guides cutting a certain
    set of genes."""
    def __init__(self, genes):
        self.genes = genes
        self.library = None
        self.name = genes[0].name + "-comp-" + str(len(self.genes))


class FastaHeader(object):
    def __init__(self, header, gene_name):
        self.header = header
        self.gene_name = gene_name

        self.resistance = None
        self.mutation_ranges = []

        self.header_parts = str(self.header).split("|")
        self.parse_header()

    def parse_header(self):
        for part in self.header_parts:
            if part.startswith("flash_resistance:"):
                self.resistance = part.split(':')[1].split(',')
            if part.startswith("flash_mutation_ranges:"):
                self.mutation_ranges = []
                for rstr in part.split('flash_mutation_ranges:')[1].split(','):
                    rrng = MutationIndex.parse_mutation(rstr)
                    if type(rrng) == range:
                        self.mutation_ranges.append((rstr, rrng))
                    else:
                        print("ERROR: {}: Failed to parse mutation_range: {}". \
                                 format(self.gene_name, rstr))


class Gene(object):
    def __init__(self, name, padding=None, mutation_ranges=[]):
        self.name = name
        self.seq = None
        self.presence_absence = None
        self.mutation_ranges = mutation_ranges

        if padding:
            self.padding = (len(padding.prefix), len(padding.suffix))
        else:
            self.padding = None

        self.resistance = None

        self.load_fasta()

        self.targets = None # List of (guide, cut location)

        # Product of cutting with a given library

        self.cuts = None # Location of each cut in sequence

        self.fragments = None # List of (start, end) for each fragment

    def load_fasta(self):
        try:
            fasta_file = open(
                os.path.join(build.genes_dir,
                             "{}.fasta".format(self.name)),
            "r")
            record = list(SeqIO.parse(fasta_file, "fasta"))[0]

            # Typical record id from CARD:
            #
            #     >gb
            #     |KF730651
            #     |+
            #     |0-657
            #     |ARO:3002796
            #     |QnrS7 [Escherichia coli]
            #     |flash_resistance:fluoroquinolone_antibiotic
            #     |flash_key:QnrS7__KF730651__ARO_3002796
            #     |flash_padding:0_200
            #
            # (except all on one line)
            #
            fasta_header = FastaHeader(record.description, self.name)
            self.resistance = fasta_header.resistance
            self.mutation_ranges = fasta_header.mutation_ranges
            self.seq = record.seq
        except FileNotFoundError:
            print(self.name, " is missing a fasta file.")

    def load_targets(self, suffix="dna_good_5_9_18.txt"):
        "Typical suffix is dna_good_5_9_18.txt"
        try:
            f = open(
                os.path.join(build.gene_index_dir,
                             self.name,
                             suffix),
            'r')

            self.targets = []
            for line in f:
                (i, d) = line.strip().split()
                i = int(i)
                kmer = flash.forward_20mer_at(self.seq, i, d)
                self.targets.append(Target(kmer, flash.cut_location((i, d))))
            self.targets.sort(key=lambda item: item.cut)
        except FileNotFoundError:
            print(self.name, " is missing a target index file.")

    def verify_padding(self, padding_seq):
        assert padding_seq[0] in self.seq
        assert padding_seq[1] in self.seq
        assert self.padding == (len(padding_seq[0]), len(padding_seq[1]))

    def get_mutation_ranges(self):
        if self.padding is None:
            return self.mutation_ranges
        else:
            padded_mutation_ranges = []
            for m, ran in self.mutation_ranges:
                # If we have added sequence of length p before the official start of the
                # gene as padding, then the location of the mutation in self.seq is p
                # after where it would be in the unpadded gene.
                #
                # TODO:  Move all this processing to make_genes_and_identify_targets.py
                # and store the results in header metadata that can be loaded above.
                # That way more of the intermediate results can be reviewed and tracked.
                p = self.padding[0]
                padded_range = range(ran.start+p, ran.stop+p, ran.step)
                padded_mutation_ranges.append((m, padded_range))
            return padded_mutation_ranges

    def length(self):
        return len(self.seq)

    def has_snps(self):
        return len(self.mutation_ranges) > 0

    def generate_fragments_from_cuts(self):
        self.fragments = [
            (self.cuts[i].cut, self.cuts[i+1].cut)for i in range(len(self.cuts) - 1)]

    def target_overlaps_mutation(self, target):
        if target.guide == self.seq[target.cut - 17 :target.cut + 3]:
            return self.range_overlaps_mutation(range(target.cut - 17, target.cut + 6))
        else:
            return self.range_overlaps_mutation(range(target.cut - 6, target.cut + 17))

    def range_overlaps_mutation(self, ran):
        for mutation, snp_range in self.get_mutation_ranges():
            for snp_loc in snp_range:
                if snp_loc in ran:
                    return True
        return False

    def cut_with_library(self, library):
        #The convention is: the fragments are given by seq[cut_1, cut_2].
        self.cuts = []
        for i in range(len(self.seq) - 23):
            kmer = self.seq[i:i+23]
            if kmer.endswith('GG'):
                target = str(kmer[:20])
                if target in library:
                    cutloc = i + 17
                    self.cuts.append(Target(target, cutloc))
                    if self.range_overlaps_mutation(range(i, i+23)):
                        print('Warning: guide %s overlaps with a mutation.' % target)
            if kmer.startswith('CC'):
                target = str(flash.reverse_complement(kmer)[:20])
                if target in library:
                    cutloc = i + 6
                    self.cuts.append(Target(target, cutloc))
                    if self.range_overlaps_mutation(range(i, i+23)):
                        print('Warning: guide %s overlaps with a mutation.' % target)

        self.cuts = sorted(self.cuts, key=lambda x: x.cut)
        self.generate_fragments_from_cuts()

    def coverage(self, min_fragment_size, max_fragment_size):
        covered_bases = 0
        for fragment in self.fragments:
            length = fragment[1] - fragment[0]
            if length >= min_fragment_size and length <= max_fragment_size:
                covered_bases += min(300, length)
        return covered_bases

    def possible_fragments(self):
        return (self.targets[-1].cut - self.targets[0].cut)//IDEAL_CUTOFF

    def longest_possible_fragment(self):
        return self.targets[-1].cut - self.targets[0].cut

    def stats(self):
        short_fragments = 0
        ideal_fragments = 0
        okay_fragments = 0
        long_fragments = 0

        for fragment in self.fragments:
            length = fragment[1] - fragment[0]
            if length < IDEAL_CUTOFF:
                short_fragments += 1
            elif length < OKAY_CUTOFF:
                ideal_fragments += 1
            elif length < LONG_CUTOFF:
                okay_fragments += 1
            else:
                long_fragments += 1
        ideal_and_okay_fragments = ideal_fragments + okay_fragments

        return {
            "gene_length": self.length(),
            "cuts": len(self.cuts),
            "ideal_fragments": ideal_fragments,
            "okay_fragments": okay_fragments,
            "ideal_and_okay_fragments": ideal_and_okay_fragments,
            "coverage": self.coverage(IDEAL_CUTOFF, LONG_CUTOFF),
            "possible_fragments": self.possible_fragments()
        }

    def display_gene_targets(self):
        arr = []
        for i in range(len(self.seq)):
            arr.append(["white", self.seq[i]])

        for mutation, snp_range in self.get_mutation_ranges():
            for i in snp_range:
                arr[i][0] = 'yellow'
                arr[i][1] = 'x'

        # Reverse order so that the indices of already-inserted cuts don't push the later indices to
        # the wrong place.

        library_guides = {} if self.cuts is None else {cut.guide for cut in self.cuts}

        for target in self.targets[::-1]:
            if self.target_overlaps_mutation(target):
                arr.insert(target.cut, [170, '|'])
            elif target.guide in library_guides:
                arr.insert(target.cut, ['red', '|'])
            else:
                arr.insert(target.cut, ['cyan', '|'])

        for (i, (col, char)) in enumerate(arr):
            if i % 100 == 0:
                sys.stdout.write("\n")
            sys.stdout.write(color(char, bg = col))
        sys.stdout.flush()
        print()

    def display_gene_cuts(self):
        if self.cuts is None or self.fragments is None:
            print("Need to cut gene first.")
            return 0

        arr = []
        for i in range(len(self.seq)):
            arr.append(["white", self.seq[i]])

        for fragment in self.fragments:
            for i in range(fragment[0],
                           min(fragment[0]+READ_LENGTH, fragment[1])):
                arr[i][0] = 'green'
            for i in range(max(fragment[0], fragment[1]-READ_LENGTH),
                           fragment[1]):
                arr[i][0] = 'green'

        for _, snp_range in self.get_mutation_ranges():
            for i in snp_range:
                arr[i][0] = 'yellow'
                arr[i][1] = 'x'

        # Reverse order so that the indices of already-inserted cuts don't push the later indices
        # to the wrong place.
        for cut in self.cuts[::-1]:
            arr.insert(cut.cut, ['red', '|'])
            #arr[cut.cut][0] = 'red'
            #arr[cut.cut][1] = 'X'

        for (i, (col, char)) in enumerate(arr):
            if i % 100 == 0:
                sys.stdout.write("\n")
            sys.stdout.write(color(char, bg=col))
        sys.stdout.flush()
        print()

    def trim_library(self, n_cuts):
        # Trim the library to the extent possible,
        # keeping at least n_cuts segments per gene.
        # Returns the trimmed library.
        guides = set()

        if self.mutation_ranges:
            # we keep all cuts for genes with SNPs
            num_cuts = len(self.cuts)
        else:
            num_cuts = min(n_cuts, len(self.cuts))

        for cut in range(num_cuts):
            guides.add(self.cuts[cut].guide)
        return list(guides)


class MutationIndex(object):
    def __init__(self):
        self.mutations = {}

    @staticmethod
    def parse_mutation(m):
        s = re.search(r"^[A-Z-](\d+)[A-Z-]$", m)
        if s:
            # eg: E502Q
            b = int(s.group(1))
            return range(b*3-3, b*3)
        else:
            # eg: nt420+2:GG
            s1 = re.search(r"^nt(\d+)\+(\d+)", m)
            if s1:
                x = int(s1.group(1)) - 1
                y = int(s1.group(2))
                return range(x, x+y)
            else:
                # eg: Y99STOP
                s2 = re.search(r"^[A-Z-](\d+)(STOP|fs)", m)
                if s2:
                    b = int(s2.group(1))
                    return range(b*3-3, b*3)
                else:
                    # eg: +nt349:CACTG
                    s3 = re.search(r"^\+nt(\d+):([A-Z-]+)$", m)
                    if s3:
                        b = int(s3.group(1)) - 1
                        return range(b, b+2)
                    else:
                        print("Mutation not parsed: ", m)
                        return None
