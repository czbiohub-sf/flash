import sys
import pytest

sys.path.append('..')
from common import Gene


def test_padding_is_inserted():
    # TODO: Implement.
    pass


def test_padding_is_parsed_corrently():
    gene = Gene("test_gene_with_padding")

    assert gene.padding == (5, 3)
    gene.verify_padding(("AAAAA", "TTT"))
