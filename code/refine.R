#!/usr/bin/env Rscript

# Example
# nanodisco refine -p 4 -b test_EC -m GATC,CCWGG,AACNNNNNNGTGC,GCACNNNNNNGTT -M all -d /home/nanodisco/dataset/EC_difference.RDS -o /home/nanodisco/analysis -r /home/nanodisco/reference/Ecoli_K12_MG1655_ATCC47076.fasta

load.libraries.refine <- function(){
	library(stringr)
	library(doMC)
	library(Biostrings)
	library(GenomicRanges)
	library(plyr)
	library(dplyr)
	library(ggplot2)
	library(RColorBrewer)
}

check.input.refine <- function(opt){
	# Check if current differences file exist.
	if(is.null(opt$path_diff_data)){
		cat("Parameter -d/--path_diff_data is missing. Please provide the path to a current differences file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$path_diff_data)){
			cat(paste0("Current differences file doesn't exist (",opt$path_diff_data,").\n"))
			cat("Please check -d/--path_diff_data parameter and/or run nanodisco difference.\n")
			quit(save="no", status=4)
		}
	}
	# Check if output directory exist and create if not.
	if(!dir.exists(opt$path_output)){
		dir.create(opt$path_output, recursive=TRUE)
	}
	# Check list of motif
	if(!is.null(opt$list_motif)){
		# Check list of motif formating
		if(!grepl(pattern="^[ACGTRYSWKMBDHVN,]+$", x=opt$list_motif)){
			cat(paste0("Unknown character found in comma separated list of motifs (",opt$list_motif,").\n"))
			cat("Please check -m/--list_motif parameter. Only the following characters are recognized: 'ACGTRYSWKMBDHVN,' (e.g. GATC,CCWGG).\n")
			quit(save="no", status=4)
		}else{
			opt$list_motif <- str_split(opt$list_motif, ",", simplify=TRUE)[1,]
		}
	}else{
		cat("Parameter -m/--list_motif is missing. Please provide a comma separated list of motifs following IUPAC nucleotide code (e.g. GATC,CCWGG).\n")
		quit(save="no", status=3)
	}
	# Check list of motif
	if(!is.null(opt$candidate_motif)){
		# Check list of motif formating
		if(opt$candidate_motif=="all"){
			opt$candidate_motif <- opt$list_motif
		}else if(!grepl(pattern="^[ACGTRYSWKMBDHVN,]+$", x=opt$candidate_motif)){
			cat(paste0("Unknown character found in comma separated list of motifs (",opt$candidate_motif,").\n"))
			cat("Please check -M/--candidate_motif parameter. Only the following characters are recognized: 'ACGTRYSWKMBDHVN,' (e.g. GATC,CCWGG).\n")
			quit(save="no", status=4)
		}else{
			opt$candidate_motif <- str_split(opt$candidate_motif, ",", simplify=TRUE)[1,]
		}
	}else{
		cat("Parameter -M/--candidate_motif is missing. Please provide a comma separated list of motifs following IUPAC nucleotide code or 'all' (e.g. GATC,CCWGG).\n")
		quit(save="no", status=3)
	}
	# Check if reference genome file exist.
	if(is.null(opt$genome)){
		cat("Parameter -r/--genome is missing. Please provide the path to a reference genome file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$genome)){
			cat(paste0("Reference genome file doesn't exist (",opt$genome,").\n"))
			cat("Please check -r/--genome parameter. Path to reference genome (.fasta or .fa).\n")
			quit(save="no", status=4)
		}
	}

	return(opt)
}

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
