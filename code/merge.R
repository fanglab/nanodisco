#!/usr/bin/env Rscript

# Example
# nanodisco merge -d /home/nanodisco/analysis/difference_subset -o /home/nanodisco/analysis -b EC_subset

suppressMessages(library(optparse))

# Parsing arguments
option_list <- list(
	make_option(c("-d", "--path_diff_data"), type="character", default=NULL, help="Path to current differences directory (*.rds produced from nanodisco difference)", metavar="<path>"),
	make_option(c("-o", "--path_output"), type="character", default="./", help="Path to output directory (default is ./)", metavar="<path>"),
	make_option(c("-b", "--base_name"), type="character", default="results", help="Base name for outputing results (e.g. Ecoli_K12; default is 'results')", metavar="<name>")
)

default_usage <- c("nanodisco merge -d <path_difference> -o <path_output> -b <analysis_name>")
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# Load parameters and functions
path_script <- "/home/nanodisco/code/"
source(paste0(path_script,"analysis_functions.R"))
suppressMessages(load.libraries.merge())

## Check input parameters
opt <- check.input.merge(opt)

path_diff_data <- opt$path_diff_data
path_output <- paste0(gsub("/$","",opt$path_output),"/") # Make sure it's ending with a slash
base_name <- opt$base_name

# Create output file full path
path_file <- paste0(path_output, base_name, "_difference.RDS")

final_dataset <- save.merged.stat.chunks(path_diff_data, path_file) # final_dataset is saved internally
