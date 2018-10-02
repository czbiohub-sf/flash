import sys
import pytest

sys.path.append('..')
from common import Gene


def test_targets_are_in_cut_order():
    gene = Gene("test_gene_with_targets")
    gene.load_targets("dna_good_5_9_18.txt")
    
    assert [t.cut for t in gene.targets] == [25, 35]
