#!/bin/bash

echo "Running tigmint with linked reads..."

../../bin/tigmint-make tigmint draft=test_contig reads=test_linkedreads

echo "Running tigmint with long reads..."

../../bin/tigmint-make tigmint-long draft=test_contig_long reads=test_longreads span=auto G=1000000 dist=-1

echo "Done! Compare your output files to those in the 'expected_outputs' directory to ensure that the run was successful."
