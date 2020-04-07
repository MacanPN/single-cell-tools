#!/usr/bin/Rscript

library('tidyverse')
library('fs')
library('rprojroot')
library('glue')
library('janitor')
library('zeallot')
library('seuratTools')

expression_table <- read_csv("resources/example_input_files/transcripts.tpm_census_matrix-comma-delimited.csv