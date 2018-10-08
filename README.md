_The code version associated with the [bioArxiv preprint](https://www.biorxiv.org/content/early/2018/09/27/426338)
is [tag 1.0](https://github.com/czbiohub/flash/releases/tag/v1.0)_

# FLASHit

FLASH (Finding Low Abundance Sequences by Hybridization) is a method
using the targeted nuclease Cas9 to enrich for an arbitrary set of sequences
in a complex DNA library. FLASH-NGS achieves up to 5 orders of magnitude of
enrichment, sub-attomolar detection of genes of interest,
and minimal background from non-FLASH-derived reads.

FLASHit is open-source software for designing FLASH experiments.
It allows users to design guide RNA libraries for any set of target genes.
The guides are guaranteed to have good structure, and to avoid cutting a set of off-target genes,
such as human DNA in the case of a library designed to target resistant genes in pathogenic bacteria.

The problem of constructing a library efficiently, using the fewest guides
necessary to get the desired coverage of all target genes, is solved
as a mixed integer program optimization. For details, see the [formal description](docs/Guide_design.pdf).

The inputs consist of:

* fasta files of genes to be targeted
* fasta files of genes to be avoided

The output is:

* a list of RNA guides

## Prerequisites

1) Install the [Anaconda](https://www.anaconda.com/) package manager for python.

2) Clone the FLASH repository from GitHub

	`git clone https://github.com/czbiohub/flash.git`

	`cd flash`

3) Install the dependencies with conda.

	`conda env create -f environment.yml`

	`source activate flash`

3) License the Gurobi optimizer.

Get a license from [Gurobi](https://user.gurobi.com/download/licenses/free-academic). Academic licenses are free.

	grbgetkey KEYHERE

4) Install [GO](https://golang.org/doc/install).

## Workflows

This section will cover the most common use cases of creating your own library and creating the AMR library from the paper.

### Creating your own library

To generate a FLASH library targeting genes from a single `fasta` file and avoiding
the human genome and a reference _E. Coli_ strain, run

`make library TARGETS=[input.fasta] OUTPUT=library.txt`.

For example,

`make library TARGETS=tests/inputs/colistin.fasta OUTPUT=colistin_library.txt`

### Creating the AMR library

To reproduce the library construction from the paper, after installing the
prerequisites below, run

`make amr_library`.

This will produce 2 files in `generated_files/untracked`:

* `library.txt` (all optimized guides for the full AMR gene set)
* `amr_library.txt` (optimized guides restricted to the 127 genes used in the paper)

## Formatting Genes

Input genes to FLASHit should be placed in one or more fasta files. If multiple
fasta files are used, they should be in a single subdirectory of `inputs`.

The header of each gene can contain arbitrary key:value pairs separated by pipes (|).

Some keys have special meaning in FLASHit.

`flash_key`: A unique ID for the gene. If unspecified, it will be automatically
generated from the fasta header before the first pipe.
`flash_mutation_ranges`: Locations of known mutations in each gene. These sites
will be avoided as CRISPR targets, and the optimizer will attempt to generate
fragments that contain each of those loci.

Example:

```
>NCTC_8325|ARO:3003312|Staphyloccocus_aureus_NCTC_8325_parC|flash_key:parC__NCTC_8325__fluoroquinolone_additional|flash_resistance:fluoroquinolone|flash_mutation_ranges:nt240+12
GTGAGTGAAATAATTCAAGATTTATCACTTGAAGATGTTTTAGGTGATCGCTTTGGAAGATATAGTAAATATATTATTCAAGAGCGTGCATTGCCAGATGTTCGTGATGGTTTAAAACCAGTACAACGTCGTATTTTATATGCAATGTATTCAAGTGGTAATACACACGATAAAAATTTCCGTAAAAGTGCGAAAACAGTCGGTGATGTTATTGGTCAATATCATCCACATGGAGAC...
```

### Mutations

The guides may be designed to capture certain mutations present in
each gene. Guides that would hybridize to those variable regions are forbidden,
and we ensure that each variable region is contained in a fragment. (The default
	assumption is that reads are paired-end 150. We ensure that
	the mutation is contained in the first or last 150 bp of the fragment.) This will allow
you to identify which allele of a gene present in each sample.

The mutations are encoded in the fasta headers of each gene as follows.

* Amino Acid Change

`A50Q` indicates that an alanine at position 50 in the
amino acid translation of the sequence may be mutated to a glutamine. Only the
position (50 here) is used by the optimizer. If you do not know the result of
the mutation, you may use an arbitrary letter.

* Nucleic Acid Substitution/Deletion

`nt420+4:GGAT` indicates that the 4 bases beginning at position 420 (420-423),
 may be mutated to `GGAT`. The target is optional (`nt420+4` would also work).

* Nucleic Acid Insertion

`+nt349:CACTG` indicates that `CACTG` may be inserted in between bases 349 and 350.
The target is optional (`+nt349` would also work).

These are included in a comma-separated field. For example,

`flash_mutation_ranges: A50Q,Y400L,+nt349:CACTG,+nt100,+nt120:CG,nt420+4`

The optimizer will return an error if there is a gene for which one of the
mutations cannot be captured in a fragment. This happens most often for SNPs
near the end of genes, when there is no legal target between the SNP and
the end of the gene. If the neighboring genomic context is known,
it can be added by hand by pre/post-pending sequence to the gene (and adjusting
the locations for the SNPs accordingly). We organize this through
the use of padding in the case of the AMR guide set:
the file `padding.json` contains the additional
sequence for those genes requiring it.

## Visualization

To inspect a given gene--its list of potential
CRISPR cut sites, its SNPs, and the cuts made by a given library--run the
`display_genes.py` script:

`python display_genes.py nalC__NC_002516__ARO_3000818 mecA__AB033763__beta_lactamase [--library library.txt]`

## Additional details about the build stages

These steps are handled by the workflows mentioned above. This section gives additional details about the individual steps.

### Build gene files

This step creates a standardized fasta file for each gene, with metadata
tags embedded in the header. These files are deposited in `generated_files/genes`.
This is automatic for genes coming from CARD and Resfinder;
for information about how to format your own
genes appropriately consult the Gene Formatting section below.

```
make build_gene_files
```

### Compile offtargets

This step extracts and collates all CRISPR targets from the offtarget files.
For the human genome, this takes about 2 minutes and produces a ~4 GB txt file.

```
make generate_offtargets
```

### Build gene indices

This step finds all the CRISPR targets in the input genes, filtered for structure and offtargets.

```
make build_indices
```

### Run optimizer

This step finds the optimal library of guides for targeting the input genes.

`make optimizer'

To search for guides extending a given library, include that library as an additional input.
To ensure that certain guides are not included in your library, include a list of those guides as an additional input.

`python optimizer.py --extend [library.txt] --exclude [badguides.txt] --output [extended_library.txt]`

### Extract guides for a certain set of genes

`python extract_guides.py [library.txt] [gene_list.txt]`

## Git

For the paper, we keep many of the intermediate files in github for easy versioning.
All files in 'generated_files/under_version_control' are automatically committed
when generated. To review what has changed, run `git diff HEAD`.

If you would like to work with a different collection of genes,
SNPs, or offtargets, we recommend working in a new git branch. (To create a new branch,
run `git checkout -b BRANCHNAME`).

## AMR Use Case

FLASHit has been used to design guides to target antimicrobial resistance (AMR)
genes. Those genes were gathered from:

1) The Comprehensive Antibiotic Resistance Database ([CARD](https://card.mcmaster.ca/)) from McMaster University.
Those genes are in `inputs/card`.
2) [Resfinder](https://cge.cbs.dtu.dk/services/ResFinder/) from the Center for Genomic Epidemiology.
Those genes are in `inputs/resfinder`.
3) Additional resistance genes gathered from the literature.
Those genes are in `inputs/additional`.

AMR genes in the FLASH pipeline are canonically named with unique keys
that look like this

	GES-11__FJ854362__ARO_3002340
	ermD__M29832__macrolide
	fdg2__NC_000962.3__delamanid_additional

Each such key consists of 3 fields separated by double underscore.

The first field is a gene name, in the above example GES-11 or ermD or fdg2.

The second field is a LOCUS or accession number of some sort.

The third field is either

   -- an ARO number for CARD genes, or

   -- a Resfinder file name (usually an antibiotic name) for Resfinder genes, or

   -- a tag usually consisting of an antibiotic name and the word "additional",
	  for genes coming from the `inputs/additional/*.fasta` files

 These keys are inferred for genes from CARD and Resfinder, or user-supplied
 for genes from inputs/additional.   The rules are strictly enforced, so if
 you are adding manually a gene under inputs/additional, and you neglect to
 specify a valid unique key or resistance tag, your gene will be rejected.

 Here is an example.  Suppose you want to add the gene

 >NC_000962.3:490783-491793_fdg1 Mycobacterium tuberculosis H37Rv, complete genome
GTGGCTGAACTGAAGCTAGGTTACAAAGCATCGGCCGAACAATTCGCACCGCGCGAGCTCGTCGAACTAG
CCGTCGCCGCCGAAGCCCACGGCATGGACAGCGCGACCGTCAGCGACCATTTTCAGCCTTGGCGCCACCA
...

and you have the above in file inputs/additional/myantibiotic.fasta.  It will
be rejected unless you extend the header with two more tags, like so

>NC_000962.3:490783-491793_fdg1 Mycobacterium tuberculosis H37Rv, complete genome|flash_key:fdg1__NC_000962.3__myantibiotic_additional|flash_resistance:myantibiotic
GTGGCTGAACTGAAGCTAGGTTACAAAGCATCGGCCGAACAATTCGCACCGCGCGAGCTCGTCGAACTAG
CCGTCGCCGCCGAAGCCCACGGCATGGACAGCGCGACCGTCAGCGACCATTTTCAGCCTTGGCGCCACCA
...

"myantibiotic" should be replaced with the specific antibiotic that
this gene confers resistance to, and should match the file name as
well as the flash_resistance key added to the header above.

Note the addition of the flash_key with value

	fdg1__NC_000962.3__myantibiotic_additional

where again myantibiotic should be replaced with the name of the specific
antibiotic for this gene.
