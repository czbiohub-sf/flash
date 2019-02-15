#!/usr/bin/env bash

set -e

TARGETS_FASTA=$1
EXPECTED_GUIDES=$2

pushd ../
make library TARGETS=tests/$TARGETS_FASTA OUTPUT=tests/library.txt 1>/dev/null
popd

diff library.txt $EXPECTED_GUIDES || (echo "Test failed, guides don't match!" && exit 1)
echo "Tests passed!"
rm library.txt
