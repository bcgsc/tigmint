#!/bin/bash

echo "Running tigmint-make tigmint-long-cut..."

../../bin/tigmint-make tigmint-long draft=test_contig reads=test_longreads span=2

echo "Done! Compare your output files to those in the 'expected_outputs' directory to ensure that the run was successful."
