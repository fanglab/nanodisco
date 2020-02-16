#!/usr/bin/env Rscript

# Example
# nanodisco filter_profile -p 40 -r /home/nanodisco/reference/metagenome.fasta -d /home/nanodisco/dataset/metagenome_difference.RDS -f /home/nanodisco/dataset/selected_features_MGM1.RDS -b MGM1 -o ./

suppressMessages(library(optparse))

# Parsing arguments
option_list <- list(
	make_option(c("-p", "--nb_threads"), type="integer", default=1, help="Number of threads to use", metavar="<integer>"),
	make_option(c("-r", "--metagenome"), type="character", default=NULL, help="Path to reference metagenome (.fasta)", metavar="<path>"),
	make_option(c("-d", "--path_difference"), type="character", default=NULL, help="Path to current differences file (*.RDS produced from nanodisco difference)", metavar="<path>"),
	make_option(c("-f", "--path_feature"), type="character", default=NULL, help="Path to selected features file (*.RDS produced from nanodisco selected_feature)", metavar="<path>"),
	make_option(c("-b", "--base_name"), type="character", default="results", help="Base name for outputing results (e.g. Ecoli_K12; default is 'results')", metavar="<character>"),
	make_option(c("-o", "--path_output"), type="character", default="./", help="Path to output directory (default is ./)", metavar="<path>"),
	make_option(c("-c", "--min_pos_coverage"), type="integer", default=10, help="Minimum coverage/number of current values needed at given position for methylation feature computation (default is 10)", metavar="<integer>")
)

default_usage <- c("nanodisco filter_profile -p <nb_threads> -r <path_fasta> -d <path_difference> -f <path_feature> -b <analysis_name> -o <path_output>")
opt_parser <- OptionParser(usage=default_usage, option_list=option_list)
opt <- parse_args(opt_parser)

# Load parameters and functions
source("/home/nanodisco/code/metagenome_analysis_functions.R")
suppressMessages(load.libraries.metagenome())

## Check input parameters
opt <- check.input.filter(opt)

nb_threads <- opt$nb_threads
metagenome <- opt$metagenome
path_difference <- opt$path_difference
path_feature <- opt$path_feature
base_name <- opt$base_name
path_output <- paste0(gsub("/$","",opt$path_output),"/") # Make sure it's ending with a slash
min_pos_coverage <- opt$min_pos_coverage

# Option for development only
genomic_range_expected_signal <- NA
target_length <- NA

print_message("Prepare default metagenome annotation")

# Prepare default metagenome annotation
metagenome_annotation <- readDNAStringSet(metagenome)
metagenome_annotation <- data.frame(contig=names(metagenome_annotation), length=width(metagenome_annotation), id=NA) # Simple metagenome contig annotation

print_message("Load supplied current differences")

# Load supplied current differences
metagenome_diff_data <- readRDS(path_difference) # Read saved full current differences dataset

print_message("Load supplied selected features")

# Load supplied current differences
selected_features <- readRDS(path_feature)

print_message("Retrieve contigs coverage information")

# Retrieve contigs coverage information
contig_coverage <- attr(selected_features, "contig_coverage")

print_message("Compute informative methylation features on all contigs")

# Define subset of motif to process
list_motifs_filtered <- unique(as.character(selected_features$motif))

# Compute features all contigs; ~2 hour on 2 threads for the example (high memory usage)
methylation_profile <- score.metagenome.motifs(metagenome_annotation, metagenome, "real", list_motifs_filtered, "dist_real", metagenome_diff_data, genomic_range_expected_signal, min_pos_coverage, "any", target_length, iupac_nc, nb_threads) # Can be avoided, see below

print_message("Select significant features")

# Select significant features only
methylation_profile_filtered <- filtering.sigOnly.features(methylation_profile, selected_features, "all") # Not necessary when using list of de novo discovered motifs.	
attr(methylation_profile_filtered, "contig_coverage") <- contig_coverage # Add coverage information as attribute. Needed for binning and avoid later confusion.

print_message("Save methylation profile matrix")

# Save methylation profile matrix
saveRDS(methylation_profile_filtered, file=paste0(path_output,"methylation_profile_",base_name,".RDS"))

print_message("Done")




