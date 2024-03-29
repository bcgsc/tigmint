#!/bin/bash

echo "Please ensure that tigmint-make is in your PATH!"

echo "Running tigmint with linked reads..."

tigmint-make -B tigmint draft=test_contig reads=test_linkedreads

echo "Running tigmint with long reads..."

tigmint-make -B tigmint-long draft=test_contig_long reads=test_longreads span=10 G=7000 dist=auto

echo "Done! Compare your output files to those in the 'expected_outputs' directory to ensure that the run was successful."
