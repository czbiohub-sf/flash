#! /bin/bash

RNA_FASTA=$1
sed '/^[^>]/ y/uU/tT/' $RNA_FASTA
