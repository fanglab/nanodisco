#!/usr/bin/env Rscript

# Example
# nanodisco motif -p 40 -b test_EC -d /home/nanodisco/dataset/EC_difference.RDS -o /home/nanodisco/analysis -r /home/nanodisco/reference/Ecoli_K12_MG1655_ATCC47076.fasta

suppressMessages(library(optparse))

# Parsing arguments
option_list <- list(
	make_option(c("-p", "--nb_threads"), type="integer", default=1, help="Number of threads to use", metavar="<integer>"),
	make_option(c("-b", "--base_name"), type="character", default="results", help="Base name for outputing results (e.g. Ecoli_K12; default is 'results')", metavar="<name>"),
	make_option(c("-d", "--path_diff_data"), type="character", default=NULL, help="Path to current differences file (*.RDS produced from nanodisco difference)", metavar="<path>"),
	make_option(c("-o", "--path_output"), type="character", default="./", help="Path to output directory (default is ./)", metavar="<path>"),
	make_option(c("-m", "--list_motif"), type="character", default=NULL, help="Comma separated list of motifs following IUPAC nucleotide code (e.g. GATC,CCWGG)", metavar="<motif1,motif2,...>"),
	make_option(c("-r", "--genome"), type="character", default=NULL, help="Path to reference genome (.fasta)", metavar="<path>"),
	make_option(c("-t", "--threshold"), type="double", default=NA, help="Smoothed peaks p-values threshold for sequence selection (if double: peaks > <threshold> or if NA: top <nb_peaks> only, default is NA)", metavar="<double>"),
	make_option(c("-c", "--list_contig"), type="character", default=NULL, help="Comma separated list of contigs (e.g. contig_1,contig_3)", metavar="<contig_1,contig_3,...>"),
	make_option(c("--contigs_file"), type="character", default=NULL, help="Path to file with list of contigs (one per line)", metavar="<path>"),
	make_option(c("-a","--auto"), type="character", default=FALSE, action="store_true", help="Disable manual motif discovery procedure (not recommended; default is FALSE)"),
	make_option(c("--score_threshold"), type="double", default=2, help="Threshold used in motif refinement (default is 2)", metavar="<double>"),
	make_option(c("--nb_peaks"), type="integer", default=2000, help="Number of sequence with p-value peaks to keep for each round (default is 2000)", metavar="<integer>"),
	make_option(c("--stat_type"), type="character", default="u_test_pval", help="Select which type of p-value sources used (default is u_test_pval)", metavar="<character>"),
	make_option(c("--smooth_func"), type="character", default="mySumlog", help="Function to use for p-values smoothing (default is mySumlog)", metavar="<character>"),
	make_option(c("--smooth_win_size"), type="integer", default=5, help="Window size used for smoothing p-values (default is 5)", metavar="<integer>"),
	make_option(c("--peak_win_size"), type="integer", default=2, help="Window size used for p-values peaks detection (default is 2)", metavar="<integer>")
)

default_usage <- c("nanodisco motif -p <nb_threads> -b <analysis_name> -d <path_difference> -o <path_output> -r <path_fasta>")
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# Load parameters and functions
path_script <- "/home/nanodisco/code/"
source(paste0(path_script,"analysis_functions.R"))
suppressMessages(load.libraries.motif())

## Check input parameters
opt <- check.input.motif(opt)

nb_threads <- opt$nb_threads
base_name <- opt$base_name
path_diff_data <- opt$path_diff_data
genome <- opt$genome
path_output <- paste0(gsub("/$","",opt$path_output),"/") # Make sure it's ending with a slash
list_motif <- opt$list_motif # Potential conversion from comma separated list to vector made in check.input.motif
smooth_win_size <- opt$smooth_win_size
peak_win_size <- opt$peak_win_size
stat_type <- opt$stat_type
nb_peaks <- opt$nb_peaks
smooth_func <- opt$smooth_func
threshold <- opt$threshold
automated <- opt$auto # Default is FALSE
score_threshold <- opt$score_threshold # Default is FALSE
list_contig <- opt$list_contig

print_message("Prepare output folder")
if(is.na(threshold)){
	output_fasta_name=paste0("motif_detection/motifs_22bp_top_",nb_peaks,"_peaks.fasta")
}else{
	output_fasta_name=paste0("motif_detection/motifs_22bp_",paste0("th",threshold),"_",nb_peaks,"_peaks.fasta")
}

# Aggregate motif detection parameters
detection_params <- list(
	smooth_win_size=smooth_win_size, peak_win_size=peak_win_size, stat_val=stat_type, nb_peaks=nb_peaks, smooth_func=smooth_func,
	threshold=threshold, seq_params=seq_params, automated=automated, score_threshold=score_threshold, output_fasta_name=output_fasta_name,
	meme_main_dir="motif_detection/", meme_analysis_dir=paste0("meme_",base_name,"/"), nbCPU=nb_threads
)

print_message("Load supplied current differences")

difference_data <- readRDS(path_diff_data)

if(!is.na(list_contig)){
	print_message("Select subset of contigs")

	difference_data <- subset(difference_data, contig %in% list_contig)
}

print_message("Detect motifs")

discovered_motifs <- wrapper.motif.detection(difference_data, genome, path_output, detection_params, iupac_nc, list_motif) # list_motif can be ignore if NULL
cat(paste0(c("Final results:\n",paste0(discovered_motifs, collapse=","),"\n"), collapse=""))

print_message("Done")
