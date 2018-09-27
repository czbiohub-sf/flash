import sys
import pytest

sys.path.append('..')
from common import Gene


def test_aa_snps_are_in_correct_location():
    gene = Gene("test_gene_with_aa_snp")
    mutation_ranges = gene.mutation_ranges
    expected_mutation = [('Y151V', range(450, 453))]

    assert mutation_ranges == expected_mutation


def test_nt_snps_are_in_correct_location():
    gene = Gene("test_gene_with_nt_snp")
    mutation_ranges = gene.mutation_ranges
    expected_mutation = [('nt1+2:AT', range(0, 2))]

    assert mutation_ranges == expected_mutation


def test_nt_insert_are_in_correct_location():
    gene = Gene("test_gene_with_insert_nt_snp")
    mutation_ranges = gene.mutation_ranges
    expected_mutation = [('+nt1:A', range(0, 2))]

    assert mutation_ranges == expected_mutation    

def test_stop_codon_in_correct_location():
    gene = Gene("test_gene_with_stop_codon")
    mutation_ranges = gene.mutation_ranges
    expected_mutation = [('Y151STOP', range(450, 453))]

    assert mutation_ranges == expected_mutation    
    
