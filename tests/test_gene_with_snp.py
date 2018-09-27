import sys
import pytest

sys.path.append('..')
from common import Gene


def test_nt_snps_are_in_correct_location():
    gene = Gene("test_gene_with_snp")
    mutation_ranges = gene.mutation_ranges
    expected_mutation = [('Y151V', range(450, 453))]

    assert mutation_ranges == expected_mutation
