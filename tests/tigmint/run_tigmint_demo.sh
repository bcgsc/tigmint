#!/bin/bash

echo "Running tigmint-make tigmint..."

../../bin/tigmint-make tigmint draft=test_contig reads=test_linkedreads

echo "Done! Compare your output files to those in the 'expected_outputs' directory to ensure that the run was successful."
