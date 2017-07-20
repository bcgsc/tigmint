#!/usr/bin/env Rscript

library(dplyr)
library(readr)

args = commandArgs(trailingOnly=TRUE)
reads_filename <- args[1]
molecules_filename <- args[2]

# Read data
reads_orig <- read_tsv(reads_filename) %>% filter(!is.na(BX))

# Count number of reads per barcode
reads_per_barcode <- reads_orig %>% count(BX) %>% rename(Reads = n)

# Select barcodes with 4 or more reads
reads_per_barcode_threshold <- 4

reads <- reads_orig %>%
	semi_join(filter(reads_per_barcode, Reads >= reads_per_barcode_threshold), by = "BX")

# Group reads into molecules
distance_threshold <- 50000

molecules_orig <- reads %>%
	arrange(Rname, BX, Pos) %>%
	mutate(MI = cumsum(BX != lag(BX, default = "NA") | Rname != lag(Rname) | Pos - lag(Pos) >= distance_threshold))

# Count number of reads per molecule
reads_per_molecule <- molecules_orig %>% count(MI) %>% rename(Reads = n)

# Select molecules with 4 or more reads
reads_per_molecule_threshold <- 4

molecules <- molecules_orig %>%
	semi_join(filter(reads_per_molecule, Reads >= reads_per_molecule_threshold), by = "MI") %>%
	mutate(MI = as.integer(factor(MI)))

# Write molecule identifiers to a TSV file
molecules %>% write_tsv(molecules_filename)

# Number of barcodes and molecules
cat(
	"Barcodes", sum(reads_per_barcode$Reads >= reads_per_barcode_threshold),
	"Molecules", max(molecules$MI))
