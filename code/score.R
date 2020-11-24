#!/usr/bin/env Rscript

# Example
# nanodisco score -r /home/nanodisco/reference/Ecoli_K12_MG1655_ATCC47076.fasta -d /home/nanodisco/dataset/EC_difference.RDS -b Ecoli -o analysis/Ecoli_motifs -m GATC,CCWGG,GCACNNNNNNGTT,AACNNNNNNGTGC

suppressMessages(library(optparse))

# Parsing arguments
option_list <- list(
	make_option(c("-p", "--nb_threads"), type="integer", default=1, help="Number of threads to use (default is 1)", metavar="<integer>"),
	make_option(c("-d", "--path_diff_data"), type="character", default=NULL, help="Path to current differences file (*.RDS produced from nanodisco difference)", metavar="<path>"),
	make_option(c("-r", "--path_reference"), type="character", default=NULL, help="Path to reference metagenome (.fasta)", metavar="<path>"),
	make_option(c("-b", "--base_name"), type="character", default="results", help="Base name for outputing results (e.g. Ecoli_K12; default is 'results')", metavar="<character>"),
	make_option(c("-o", "--path_output"), type="character", default="./", help="Path to output directory (default is ./)", metavar="<path>"),
	make_option(c("-m", "--list_motif"), type="character", default=NULL, help="Comma separated list of motifs following IUPAC nucleotide code (e.g. GATC,CCWGG)", metavar="<motif1,motif2,...>"),
	make_option(c("-c", "--list_contig"), type="character", default=NULL, help="Comma separated list of contigs (e.g. contig_1,contig_3)", metavar="<contig_1,contig_3,...>"),
	make_option(c("--contigs_file"), type="character", default=NULL, help="Path to file with list of contigs (one per line)", metavar="<path>")
)

default_usage <- c("nanodisco score -p <nb_threads> -b <analysis_name> -d <path_difference> -o <path_output> -m <motif1,motif2,...> -r <path_fasta>")
opt_parser <- OptionParser(usage=default_usage, option_list=option_list)
opt <- parse_args(opt_parser);

# Load parameters and functions
source("/home/nanodisco/code/metagenome_analysis_functions.R")
suppressMessages(load.libraries.score())

## Check input parameters
check.input.score <- function(opt){
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
	if(is.null(opt$list_motif)){
		cat("Parameter -m/--list_motif is missing. Please provide a comma separated list of motifs (e.g. GATC,CCWGG).\n")
		quit(save="no", status=3)
	}else{
		# Check list of motif formating
		if(!grepl(pattern="^[ACGTRYSWKMBDHVN,]+$", x=opt$list_motif)){
			cat(paste0("Unknown character found in comma separated list of motifs (",opt$list_motif,").\n"))
			cat("Please check -m/--list_motif parameter. Only the following characters are recognized: 'ACGTRYSWKMBDHVN,' (e.g. GATC,CCWGG).\n")
			quit(save="no", status=4)
		}
	}
	# Check if reference genome file exist.
	if(is.null(opt$path_reference)){
		cat("Parameter -r/--path_reference is missing. Please provide the path to a reference genome file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$path_reference)){
			cat(paste0("Reference genome file doesn't exist (",opt$path_reference,").\n"))
			cat("Please check -r/--path_reference parameter. Path to reference genome (.fasta or .fa).\n")
			quit(save="no", status=4)
		}
	}
	# Check if a subset of contigs should be analyzed.
	if(is.null(opt$list_contig) && is.null(opt$contigs_file)){
		tmp_reference <- readDNAStringSet(opt$path_reference)
		if(length(tmp_reference) > 1){
			cat("Parameters -c/--list_contig or --contigs_file were not supplied. All contigs will be processed.\n")
		}
		opt$list_contig <- str_split(names(tmp_reference), " ", simplify=TRUE)[,1]
	}else{
		if(!is.null(opt$list_contig) && !is.null(opt$contigs_file)){
			cat("Parameters -c/--list_contig and --contigs_file cannot be both supplied. Please select one of the two options.\n")
			quit(save="no", status=4)
		}else{
			if(!is.null(opt$list_contig)){
				# TODO could check if contig exist in metagenome.fasta and/or in methylation profile file
			}else if(!is.null(opt$contigs_file)){
				if(!file.exists(opt$contigs_file)){
					cat(paste0("List of contigs file doesn't exist (",opt$contigs_file,").\n"))
					cat("Please check --contigs_file parameter. Path to file with of contigs (one per line).\n")
					quit(save="no", status=4)
				}else{
					contigs_from_file <- read.table(opt$contigs_file, stringsAsFactors=FALSE)
					# TODO could check if contig exist in metagenome.fasta and/or in methylation profile file
					opt$list_contig <- paste0(contigs_from_file$V1, collapse=",")
				}
			}
		}
	}

	return(opt)
}
opt <- check.input.score(opt)

path_diff_data <- opt$path_diff_data
path_reference <- opt$path_reference
base_name <- opt$base_name
path_output <- paste0(gsub("/$","",opt$path_output),"/") # Make sure it's ending with a slash
list_motif <- str_split(opt$list_motif, ",", simplify=TRUE)[1,]
list_contig <- opt$list_contig
nb_threads <- opt$nb_threads

print_message("Load supplied current differences")

# Load supplied current differences
difference_data <- readRDS(path_diff_data) # Read saved full dataset

list_motif <- unique(list_motif) # Remove duplicate if used multiple lists

print_message("Compute scores for each motif occurrence in each requested contig")

registerDoMC(nb_threads)
scored_motif <- foreach(contig_name=list_contig, .combine=rbind) %do% {
	contig_scored_motif <- foreach(motif_to_score=list_motif, .combine=rbind) %dopar% {
		background_motifs <- list_motif[!grepl(motif_to_score, list_motif)]
		subset_scored_motif <- score.contig.motif(motif_to_score, background_motifs, contig_name, difference_data, path_reference, NA, "top_win", FALSE)
		subset_scored_motif$motif <- as.character(subset_scored_motif$motif)

		return(subset_scored_motif)
	}

	return(contig_scored_motif)
}
registerDoSEQ()
scored_motif$motif <- as.factor(scored_motif$motif)
scored_motif$mean_diff <- NULL
scored_motif$distance <- NULL
scored_motif$N_wga <- NULL
scored_motif$N_nat <- NULL
scored_motif$cov_wga <- round(scored_motif$cov_wga, 2)
scored_motif$cov_nat <- round(scored_motif$cov_nat, 2)
scored_motif$score <- round(scored_motif$score, 4)

print_message("Writing results")

file_name <- paste0("Motifs_scores_",base_name,".tsv")
output_file_name <- paste0(path_output,file_name)

write.table(scored_motif, output_file_name, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

print_message("Done")



