# FLASHit Design

## Main Function

flashit.py

Make a FLASH library.

Arguments:
    genes (fasta): file containing all target genes
    capture_regions (bed): file indicating nucleotides in target to capture
    variable_regions (bed): file indicating additional nucleotides to avoid w/ CRISPR
    offtarget (fasta): fasta file(s) containing sequence not to cut
    use-cache (bool): if true, use same offtargets as last time (if already computed)

    config (yaml): describes parameters like read length, max fragment size,
    GC cutoffs for structural filter, etc.

Output:
    library (txt): list of generated guides
    cuts (bed): locations of cuts, regions targeted by guides

## Structure

All intermediates will go in `tmp` directory unless otherwise specified.
This will allow for caching, but not require checking the files in/out of github.
