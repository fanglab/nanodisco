#!/usr/bin/env Rscript

# Example
# nanodisco binning -r /home/nanodisco/reference/metagenome.fasta -s /home/nanodisco/analysis/methylation_profile_MGM1.RDS -b MGM1 -o ./

suppressMessages(library(optparse))

# Parsing arguments
option_list <- list(
	make_option(c("-r", "--metagenome"), type="character", default=NULL, help="Path to reference metagenome (.fasta)", metavar="<path>"),
	make_option(c("-s", "--path_profile"), type="character", default=NULL, help="Path to methylation profile file (*.RDS produced from nanodisco profile)", metavar="<path>"),
	make_option(c("-b", "--base_name"), type="character", default="results", help="Base name for outputing results (e.g. Ecoli_K12; default is 'results')", metavar="<character>"),
	make_option(c("-o", "--path_output"), type="character", default="./", help="Path to output directory (default is ./)", metavar="<path>"),
	make_option(c("--min_motif_occ"), type="integer", default=5, help="Minimum number of motif occurrence to conserve entry in the methylation profile matrix (default is 5)", metavar="<integer>"),
	make_option(c("--min_contig_len"), type="integer", default=25000, help="Minimum contig length to conserve entry in the methylation profile matrix (default is 25000)", metavar="<integer>"),
	make_option(c("--contig_weight_unit"), type="integer", default=50000, help="Weight unit (bp) used for additional exageration in binning (default is 50000)", metavar="<integer>"),
	make_option(c("--max_relative_weight"), type="double", default=0.05, help="Maximum relative weight a contig can have, weighting ceiling (default is 0.05)", metavar="<double>"),
	make_option(c("--tsne_perplexity"), type="integer", default=30, help="t-SNE perplexity parameter (default is 30)", metavar="<integer>"),
	make_option(c("--tsne_max_iter"), type="integer", default=2500, help="t-SNE maximum iteration parameter (default is 2500)", metavar="<integer>"),
	make_option(c("--tsne_seed"), type="integer", default=101, help="Seed set before t-SNE processing (default is 101, set.seed function)", metavar="<integer>"),
	make_option(c("--rdm_seed"), type="integer", default=2, help="Seed used for random number generation in missing value filling (default is 2, set.seed function)", metavar="<integer>")
)

default_usage <- c("nanodisco binning -r <path_fasta> -s <path_profile> -b <analysis_name> -o <path_output> [+ advanced parameters]")
opt_parser <- OptionParser(usage=default_usage, option_list=option_list)
opt <- parse_args(opt_parser)

# Load parameters and functions
source("/home/nanodisco/code/metagenome_analysis_functions.R")
suppressMessages(load.libraries.metagenome())

## Check input parameters
opt <- check.input.binning(opt)

metagenome <- opt$metagenome
path_profile <- opt$path_profile
base_name <- opt$base_name
path_output <- paste0(gsub("/$","",opt$path_output),"/") # Make sure it's ending with a slash
min_motif_occ <- opt$min_motif_occ
min_contig_len <- opt$min_contig_len
contig_weight_unit <- opt$contig_weight_unit
max_relative_weight <- opt$max_relative_weight
tsne_perplexity <- opt$tsne_perplexity
tsne_max_iter <- opt$tsne_max_iter
tsne_seed <- opt$tsne_seed
rdm_seed <- opt$rdm_seed

print_message("Prepare default metagenome annotation")

# Prepare default metagenome annotation
metagenome_annotation <- readDNAStringSet(metagenome)
metagenome_annotation <- data.frame(contig=names(metagenome_annotation), length=width(metagenome_annotation), id=NA, stringsAsFactors=TRUE) # Simple metagenome contig annotation

print_message("Load supplied methylation profile matrix")

# Load supplied methylation profile matrix
methylation_profile <- readRDS(path_profile)

print_message("Retrieve contigs coverage information")

# Retrieve contigs coverage information
contig_coverage <- attr(methylation_profile, "contig_coverage")

print_message("Perform binning with dimentionality reduction")

# Perform binning with dimentionality reduction
methylation_binning <- tsne.motifs.score(metagenome, methylation_profile, "real", min_motif_occ, min_contig_len, "pseudoZero", metagenome_annotation, "id", contig_weight_unit, contig_coverage, max_relative_weight, tsne_perplexity, tsne_max_iter, rdm_seed, tsne_seed)

print_message("Save binning results")

# Save binning results
saveRDS(methylation_binning, file=paste0(path_output,"methylation_binning_",base_name,".RDS"))

print_message("Done")

