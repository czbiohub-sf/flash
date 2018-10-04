import sys
import pytest

sys.path.append('..')
from common import Gene


def test_targets_are_in_cut_order():
    gene = Gene("test_gene_with_targets")
    gene.load_targets("dna_good_5_9_18.txt")
    expected_target_cuts = [
        21,
        39,
        40,
        74,
        170,
        201,
        208,
        208,
        209,
        210,
        211
    ]

    assert [t.cut for t in gene.targets] == expected_target_cuts
