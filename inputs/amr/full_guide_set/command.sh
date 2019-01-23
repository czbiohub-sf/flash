#!/bin/bash

set -e

time make library_including TARGETS=/home/neevor/genes.fasta \
			    OUTPUT=library_from_all_fasta.txt \
			    INCLUDE=/scratch/all_ordered_guides.txt
