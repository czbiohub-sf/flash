import itertools
import time
import sys
import requests
from collections import defaultdict


COMPLEMENTS = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', 'M': 'M', 'R': 'R', 'S': 'S', 'K': 'K', 'Y': 'Y'}


def kmers_range(seq, k):
    # last kmer is seq[len(seq) - k : len(seq)]
    # last value of i is len(seq) - k
    return range(0, len(seq) - k + 1)


def kmers(seq, k):
    seq = str(seq) # major speed up from this line when seq is a Bio.SeqIO object
    for i in kmers_range(seq, k):
        yield (i, seq[i:i+k])


def poor_structure(forward_kmer, need_explanation=False):
    assert len(forward_kmer) == 20
    kmer = forward_kmer
    freq = character_count(kmer)
    reasons = []
    gc_freq = freq['G'] + freq['C']
    if gc_freq < 5 or gc_freq > 15:
        if not need_explanation:
            return True
        reasons.append('gc_frequency')
    if longest_consecutive_run(kmer)>5:
        if not need_explanation:
            return True
        reasons.append('homopolymer>5')
    if longest_dinucleotide_run(kmer) > 3:
        if not need_explanation:
            return True
        reasons.append('dinucleotide_repeats>3')
    hairpin = find_hairpin(kmer)
    if hairpin:
        if not need_explanation:
            return True
        reasons.append('hairpin:' + hairpin)
    if not need_explanation:
        return False
    return reasons


def character_count(kmer):
    d = defaultdict(int)
    for c in kmer:
        d[c] += 1
    return d


def longest_consecutive_run(s):
    return longest_run(s)


def longest_dinucleotide_run(kmer):
    # This looks for patterns like GTGTGT that are called 'dinucleotide repeats'
    # https://en.wikipedia.org/wiki/Tandem_repeat
    # A dinucleotide repeat containing 2n or 2n+1 bases in the original string corresponds
    # to a "longest run of zeroes" consisting of 2n-2 or 2n-1 zeroes in the transformed
    # sequence, and yields result n in this function.  For example:
    #
    #           kmer => transformed sequence => longest run => result
    #           -------------------------------------------------------------------------
    #           xGzGTGTGTuTG => [2, 0, 3, 0, 0, 0, 0, 9,  0, 11] => [0, 0, 0, 0] => 3
    #           xyTGTGTGTuTG => [2, 3, 0, 0, 0, 0, 0, 9,  0, 11] => [0, 0, 0, 0, 0] => 3
    #
    if len(kmer) < 2:
        return 0
    lr = longest_run((0 if kmer[n] == kmer[n-2] else n) for n in range(2, len(kmer)))
    return 1 + (lr // 2)


def generate_5_partitions(k):
    # The standard way to partition an integer in combinatorics is by
    # reducing 5-partitions to 4-combinations.  The output of this
    # generator is the sequence of all 5-tuples with sum k.
    k4 = k + 4
    # Combinations are generated in lexicographic order, so that
    # 0 <= k0 < k1 < k2 < k3 < k4.
    for k0, k1, k2, k3 in itertools.combinations(range(k4), 4):
        p = (
            k0,
            k1 - k0 - 1,
            k2 - k1 - 1,
            k3 - k2 - 1,
            k4 - k3 - 1
        )
        # If you just add up the tuple elements defined above, the k0 ... k3
        # cancel out, and what remains is k4 - 4, which is equal to k.
        # That's why this way of partitioning k is so elegant.
        assert sum(p) == k
        yield p


def generate_hairpin_bounds(args):
    k, min_outer, min_inner = args
    for offset_left, outer_left, inner, outer_right, offset_right in generate_5_partitions(k):
        if outer_left == outer_right and outer_left >= min_outer and inner >= min_inner:
            yield (outer_left, inner, offset_left)


def find_hairpin(kmer, hairpin_bounds={}):
    k = 20
    min_outer = 5
    min_inner = 3
    #
    # One way to visualize what's happening here is to consider that we are
    # essentially generating all 20-character sequences that look like
    #
    #       ...oooooIIIIIooooo..
    #
    # with equal numbers of 'o' characters on each side (at least 5),
    # also with at least 3 'I' inner characters, and any number (including 0)
    # of '.' characters, so that the total string length is 20.
    #
    # The 'o' characters are called outer, the 'I' characters are called inner,
    # and the '.' characters are called offset.
    #
    # Combinatorially, these patterns are 5-partitions of 20 that satisfy the
    # additional equality and range constraints
    #
    #     outer_1 == outer_2
    #     5 <= outer_1
    #     3 <= inner
    #
    # The total number of l-partitions of k is Binomial(l + k - 1, l - 1).
    # In this case, Binomial(24, 4) = 10626, but only 70 satisfy the hairpin
    # constraints.  We cache those 70.
    #
    cache_key = (k, min_outer, min_inner)
    if cache_key not in hairpin_bounds:
        hairpin_bounds[cache_key] = list(generate_hairpin_bounds(cache_key))
    for outer, inner, offset in hairpin_bounds[cache_key]:
        k0 = offset
        k1 = offset + outer
        k2 = offset + outer + inner
        k3 = offset + outer + inner + outer
        left = kmer[k0:k1]
        right = kmer[k2:k3]
        if complementary_pattern(left, right, max(outer-1, min_outer)):
            return ("-" * offset) + left + ("-" * inner) + right + ((20 - k3) * "-")
    return ""


def complementary_pattern(left, right, min_matches):
    if len(left) != len(right):
        return False
    rc = reverse_complement(right)
    diff_count = 0
    for i in range(len(left)):
        if left[i] != rc[i]:
            diff_count += 1
    return diff_count <= len(left) - min_matches


def reverse_complement(kmer):
    return ''.join(COMPLEMENTS[x] for x in reversed(kmer))


def longest_run(l):
    last = None
    c = 0
    m = 0
    for chr in l:
        if chr != last:
            last = chr
            m = max(c, m)
            c = 1
        else:
            c +=1
    m = max(c, m)
    return m


def cut_location(guide):
    (i, d) = guide
    if d == 'F':
        return i + 17
    elif d == 'R':
        return i + 6


def kmer_at(seq, i, k):
    return str(seq[i:i+k])


def forward_20mer_at(seq, i, d):
    kmer = kmer_at(seq, i, 23)
    if d == 'R':
        kmer = reverse_complement(kmer)
    return kmer[:20]


def test_longest_dinucleotide_run(challenge, desired_response):
    response = longest_dinucleotide_run(challenge)
    print(challenge, response)
    assert response == desired_response


if __name__ == '__main__':
     print("Testing hairpin logic.")
     x = list(generate_hairpin_bounds((20, 5, 3)))
     for outer, inner, offset_left in x:
        offset_right = 20 - offset_left - outer - inner - outer
        pattern = '.' * offset_left + 'o' * outer + 'I' * inner + 'o' * outer + '.' * offset_right
        print("{:2d} {:2d} {:2d} {:2d} {:2d}    {}".format(
            offset_left, outer, inner, outer, offset_right, pattern))
     assert len(x) == 70
     print("Hairpin logic test passed.")
     print("Testing longest dinucleotide repeats.")
     test_longest_dinucleotide_run("x", 0)
     test_longest_dinucleotide_run("xy", 1)
     test_longest_dinucleotide_run("xyxzzz", 1)
     test_longest_dinucleotide_run("xzxzzzzzx", 2)
     test_longest_dinucleotide_run("xyzGTGTGTuvw", 3)
     test_longest_dinucleotide_run("xGzGTGTGTuTG", 3)
     test_longest_dinucleotide_run("xyTGTGTGTuTG", 3)
     print("Longest dinucleotide repeats test passed.")
