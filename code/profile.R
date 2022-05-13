#!/usr/bin/env Rscript

# Example
# nanodisco profile -p 40 -r /home/nanodisco/reference/metagenome.fasta -d /home/nanodisco/dataset/metagenome_difference.RDS -w /home/nanodisco/dataset/metagenome_WGA.cov -n /home/nanodisco/dataset/metagenome_NAT.cov --motifs_file /home/nanodisco/dataset/list_de_novo_discovered_motifs.txt -b MGM1 -o ./

suppressMessages(library(optparse))

# Parsing arguments
option_list <- list(
	make_option(c("-p", "--nb_threads"), type="integer", default=1, help="Number of threads to use", metavar="<integer>"),
	make_option(c("-r", "--metagenome"), type="character", default=NULL, help="Path to reference metagenome (.fasta)", metavar="<path>"),
	make_option(c("-d", "--path_difference"), type="character", default=NULL, help="Path to current differences file (*.RDS produced from nanodisco difference)", metavar="<path>"),
	make_option(c("-w", "--path_cov_wga"), type="character", default=NULL, help="Path to WGA sample coverage", metavar="<path>"),
	make_option(c("-n", "--path_cov_nat"), type="character", default=NULL, help="Path to native sample coverage", metavar="<path>"),
	make_option(c("-a","--auto"), type="character", default=NULL, help="Compute methylation profile from predefined common motifs followed by filtering (automated binning; all|4mer|5mer|6mer|noBi)", metavar="<all|4mer|5mer|6mer|noBi>"),
	make_option(c("-m", "--list_motif"), type="character", default=NULL, help="Comma separated list of motifs following IUPAC nucleotide code (e.g. GATC,CCWGG)", metavar="<motif1,motif2,...>"),
	make_option(c("--motifs_file"), type="character", default=NULL, help="Path to file with list of motifs (one per line) following IUPAC nucleotide code", metavar="<path>"),
	make_option(c("-b", "--base_name"), type="character", default="results", help="Base name for outputing results (e.g. Ecoli_K12; default is 'results')", metavar="<character>"),
	make_option(c("-o", "--path_output"), type="character", default="./", help="Path to output directory (default is ./)", metavar="<path>"),
	make_option(c("-c", "--min_pos_coverage"), type="integer", default=10, help="Minimum coverage/number of current values needed at given position for methylation feature computation (default is 10)", metavar="<integer>"),
	make_option(c("--min_contig_len"), type="integer", default=100000, help="Minimum length to consider a contig for feature selection (default is 100000)", metavar="<integer>")
)

default_usage <- c("nanodisco profile -p <nb_threads> -r <path_fasta> -d <path_difference> -w <path_WGA_cov> -n <path_NAT_cov> -b <analysis_name> -o <path_output> (-a <all|4mer|5mer|6mer|noBi> || -m <motif1,motif2,...> || --motifs_file <path_motif>)")
opt_parser <- OptionParser(usage=default_usage, option_list=option_list)
opt <- parse_args(opt_parser)

# Load parameters and functions
source("/home/nanodisco/code/metagenome_analysis_functions.R")
suppressMessages(load.libraries.metagenome())

## Check input parameters
opt <- check.input.profile(opt)

nb_threads <- opt$nb_threads
metagenome <- opt$metagenome
path_difference <- opt$path_difference
path_cov_wga <- opt$path_cov_wga
path_cov_nat <- opt$path_cov_nat
list_motif <- str_split(opt$list_motif,",",simplify=TRUE)[1,] # Convert comma separated list of motifs into a vector
base_name <- opt$base_name
path_output <- paste0(gsub("/$","",opt$path_output),"/") # Make sure it's ending with a slash
min_pos_coverage <- opt$min_pos_coverage
min_contig_len <- opt$min_contig_len

# Option for development only
genomic_range_expected_signal <- NA
target_length <- NA

print_message("Prepare default metagenome annotation")

# Prepare default metagenome annotation
metagenome_annotation <- readDNAStringSet(metagenome)
metagenome_annotation <- data.frame(contig=names(metagenome_annotation), length=width(metagenome_annotation), id=NA, stringsAsFactors=TRUE) # Simple metagenome contig annotation

print_message("Load supplied current differences")

# Load supplied current differences
metagenome_diff_data <- readRDS(path_difference) # Read saved full current differences dataset

print_message("Load contigs coverage information")

# Load contigs coverage information
contig_coverage <- prepare.contig.coverage(path_cov_wga, path_cov_nat)

## Compute methylation profile matrix
# automated (-a/--auto): compute feature from predefined common motifs on subset of contigs
# motif driven (-m/--list_motif/--motifs_file): only compute feature from requested motifs
if(!is.null(opt$auto)){ # if -a/--auto is defined (then -m/--list_motif/--motifs_file is not)
	print_message(paste0("Prepare subset of contigs (>=",min_contig_len," bp)"))

	# Select longer contigs (>= --min_contig_len)
	long_metagenome <- gsub(".fasta$",paste0("_",min_contig_len,"bp.fasta"),metagenome)
	select.contigs.motifs.filtering(metagenome, min_contig_len, long_metagenome) # Create subset of metagenome
	sequence_long_metagenome <- readDNAStringSet(long_metagenome)
	metagenome_long_annotation <- subset(metagenome_annotation, contig %in% names(sequence_long_metagenome))
	metagenome_long_annotation$length <- width(sequence_long_metagenome)

	print_message("Compute methylation features on subset of contigs")

	# Score motifs in long contigs only
	methylation_profile <- score.metagenome.motifs(metagenome_long_annotation, long_metagenome, "real", list_motif, "dist_real", metagenome_diff_data, genomic_range_expected_signal, min_pos_coverage, "any", target_length, iupac_nc, nb_threads)
	attr(methylation_profile, "min_contig_len") <- min_contig_len # Add contig length threshold used.
}else{
	print_message("Compute features from supplied methylation motifs on all contigs")

	# Compute features all contigs
	methylation_profile <- score.metagenome.motifs(metagenome_annotation, metagenome, "real", list_motif, "dist_real", metagenome_diff_data, genomic_range_expected_signal, min_pos_coverage, "any", target_length, iupac_nc, nb_threads)
}
attr(methylation_profile, "contig_coverage") <- contig_coverage # Add coverage information as attribute. Needed for binning and avoid later confusion.

print_message("Save methylation profile matrix")

# Save methylation profile matrix
saveRDS(methylation_profile, file=paste0(path_output,"methylation_profile_",base_name,".RDS"))

print_message("Done")



