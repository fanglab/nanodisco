#!/usr/bin/env Rscript

# Example
# nanodisco refine -p 4 -b test_EC -m GATC,CCWGG,AACNNNNNNGTGC,GCACNNNNNNGTT -M all -d /home/nanodisco/dataset/EC_difference.RDS -o /home/nanodisco/analysis -r /home/nanodisco/reference/Ecoli_K12_MG1655_ATCC47076.fasta

suppressMessages(library(optparse))

# Parsing arguments
option_list <- list(
	make_option(c("-p", "--nb_threads"), type="integer", default=1, help="Number of threads to use", metavar="<integer>"),
	make_option(c("-b", "--base_name"), type="character", default="results", help="Base name for outputing results (e.g. Ecoli_K12; default is 'results')", metavar="<name>"),
	make_option(c("-d", "--path_diff_data"), type="character", default=NULL, help="Path to current differences file (*.RDS produced from nanodisco difference)", metavar="<path>"),
	make_option(c("-o", "--path_output"), type="character", default="./", help="Path to output directory (default is ./)", metavar="<path>"),
	make_option(c("-m", "--list_motif"), type="character", default=NULL, help="Comma separated list of motifs following IUPAC nucleotide code (e.g. GATC,CCWGG)", metavar="<motif1,motif2,...>"),
	make_option(c("-M", "--candidate_motif"), type="character", default=NULL, help="Comma separated list of motifs following IUPAC nucleotide code or 'all' (e.g. GATC,CCWGG)", metavar="<motif1,motif2,...>|all"),
	make_option(c("-r", "--genome"), type="character", default=NULL, help="Path to reference genome (.fasta)", metavar="<path>")
)

default_usage <- c("nanodisco motif -p <nb_threads> -b <analysis_name> -m <motif1,motif2,...> -M <motif1,motif2,...|all> -d <path_difference> -o <path_output> -r <path_fasta>")
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# Load parameters and functions
path_script <- "/home/nanodisco/code/"
source(paste0(path_script,"analysis_functions.R"))
suppressMessages(load.libraries.refine())

## Check input parameters
opt <- check.input.refine(opt)

nb_threads <- opt$nb_threads
base_name <- opt$base_name
path_diff_data <- opt$path_diff_data
genome <- normalizePath(opt$genome) # Need full path after setwd()
path_output <- paste0(gsub("/$","",opt$path_output),"/") # Make sure it's ending with a slash
list_motif <- opt$list_motif # Potential conversion from comma separated list to vector made in check.input.motif
candidate_motif <- opt$candidate_motif # Potential conversion from comma separated list to vector made in check.input.motif

print_message("Load supplied current differences")
difference_data <- readRDS(path_diff_data)

print_message("Generate refine plots")
stifle <- foreach(motif_of_interest=candidate_motif) %do% {
	no_return <- refine.motif(motif_of_interest, list_motif[! list_motif %in% motif_of_interest], difference_data, genome, nb_threads, seq_params, iupac_nc, path_output, base_name, FALSE)

	return(NA)
}

print_message("Done")
