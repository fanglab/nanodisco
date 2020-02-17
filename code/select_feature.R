#!/usr/bin/env Rscript

# Example
# nanodisco select_feature -p 40 -r /home/nanodisco/reference/metagenome.fasta -s /home/nanodisco/dataset/methylation_profile_MGM1_auto.RDS -b MGM1 -o ./

suppressMessages(library(optparse))

# Parsing arguments
option_list <- list(
	make_option(c("-p", "--nb_threads"), type="integer", default=1, help="Number of threads to use", metavar="<integer>"),
	make_option(c("-r", "--metagenome"), type="character", default=NULL, help="Path to reference metagenome (.fasta)", metavar="<path>"),
	make_option(c("-s", "--path_profile"), type="character", default=NULL, help="Path to methylation profile file (*.RDS produced from nanodisco profile)", metavar="<path>"),
	make_option(c("-b", "--base_name"), type="character", default="results", help="Base name for outputing results (e.g. Ecoli_K12; default is 'results')", metavar="<character>"),
	make_option(c("-o", "--path_output"), type="character", default="./", help="Path to output directory (default is ./)", metavar="<path>"),
	make_option(c("--fsel_min_cov"), type="integer", default=10, help="Minimum average coverage to consider a contig for feature selection (default is 10)", metavar="<integer>"),
	make_option(c("--fsel_min_motif_occ"), type="integer", default=20, help="Minimum number of motif occurrences in a contig for feature selection (default is 20)", metavar="<integer>"),
	make_option(c("--fsel_min_signal"), type="double", default=1.5, help="Absolute threshold for considering a feature informative (default is 1.5)", metavar="<float>"),
	make_option(c("--fsel_min_contig_len"), type="integer", default=NULL, help="Minimum length to consider a contig for feature selection (default use value from methylation profile)", metavar="<integer>")
)

default_usage <- c("nanodisco select_feature -p <nb_threads> -r <path_fasta> -s <path_profile> -b <analysis_name> -o <path_output>")
opt_parser <- OptionParser(usage=default_usage, option_list=option_list)
opt <- parse_args(opt_parser)

# Load parameters and functions
source("/home/nanodisco/code/metagenome_analysis_functions.R")
suppressMessages(load.libraries.metagenome())

print_message("Check input files and parameters")

## Check input parameters
opt <- check.input.feature(opt)

nb_threads <- opt$nb_threads
metagenome <- opt$metagenome
path_profile <- opt$path_profile
base_name <- opt$base_name
path_output <- paste0(gsub("/$","",opt$path_output),"/") # Make sure it's ending with a slash
fsel_min_cov <- opt$fsel_min_cov
fsel_min_motif_occ <- opt$fsel_min_motif_occ
fsel_min_signal <- opt$fsel_min_signal
fsel_min_contig_len <- opt$fsel_min_contig_len # or attr(methylation_profile, "min_contig_len"); Handled in check.input.feature()

# Option for development only
genomic_range_expected_signal <- NA
target_length <- NA

print_message("Prepare default metagenome annotation")

# Prepare default metagenome annotation
metagenome_annotation <- readDNAStringSet(metagenome)
metagenome_annotation <- data.frame(contig=names(metagenome_annotation), length=width(metagenome_annotation), id=NA) # Simple metagenome contig annotation

print_message("Load methylation profile")

# Load supplied methylation profile matrix
methylation_profile <- readRDS(path_profile)

print_message("Retrieve contigs coverage information")

# Retrieve contigs coverage information and minimum contig length
contig_coverage <- attr(methylation_profile, "contig_coverage")

# Select longer contigs (>= --fsel_min_contig_len)
long_metagenome <- gsub(".fasta$",paste0("_",fsel_min_contig_len,"bp.fasta"),metagenome)
if(!file.exists(long_metagenome)){select.contigs.motifs.filtering(metagenome, fsel_min_contig_len, long_metagenome)} # Create subset of metagenome
sequence_long_metagenome <- readDNAStringSet(long_metagenome)
metagenome_long_annotation <- subset(metagenome_annotation, contig %in% names(sequence_long_metagenome))
metagenome_long_annotation$length <- width(sequence_long_metagenome)

# Identify contigs with insufficient coverage (< fsel_min_cov)
contigs_low_cov <- as.character(subset(contig_coverage, contig_length>=fsel_min_contig_len & (avg_cov.dataset_A<fsel_min_cov | avg_cov.dataset_B<fsel_min_cov))$chr)

print_message("Select informative features")

# Select informative features; ~30 min for the example
features <- select.features(methylation_profile, metagenome_long_annotation, fsel_min_contig_len, fsel_min_motif_occ, fsel_min_signal)
# Remove features from contigs with low coverage
features_noLowCov <- ignore.features.from.bin(features, contigs_low_cov, metagenome_annotation)
# Filter features from bipartite motifs with fewer than 2 significants features
features_noLowCov_2posBipartite <- filtering.sig.bipartite.motifs(features_noLowCov, "at_least", 2)
# Filter features from bipartite motifs overlapping non-bipartite motifs
features_noLowCov_2posBipartite_noOverlapBi <- filtering.overlapping.motifs(features_noLowCov_2posBipartite, c(4,5,6), "bipartite", nb_threads)

selected_features <- features_noLowCov_2posBipartite_noOverlapBi
attr(selected_features, "contig_coverage") <- contig_coverage

print_message("Save selected methylation features")

# Save selected features
saveRDS(selected_features, file=paste0(path_output,"selected_features_",base_name,".RDS"))

print_message("Done")
