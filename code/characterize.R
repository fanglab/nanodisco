#!/usr/bin/env Rscript

# Example
# nanodisco characterize -p 20 -b results -d /home/nanodisco/dataset/EC_difference.RDS -o /home/nanodisco/analysis/class -m GATC,CCWGG,GCACNNNNNNGTT,AACNNNNNNGTGC -t nn,rf,knn -r /home/nanodisco/reference/Ecoli_K12_MG1655_ATCC47076.fasta

suppressMessages(library(optparse))

# Parsing arguments
option_list <- list(
	make_option(c("-p", "--nb_threads"), type="integer", default=1, help="Number of threads to use", metavar="<integer>"),
	make_option(c("-b", "--base_name"), type="character", default="results", help="Base name for outputing results (e.g. Ecoli_K12; default is 'results')", metavar="<name>"),
	make_option(c("-d", "--path_diff_data"), type="character", default=NULL, help="Path to current differences file (*.RDS produced from nanodisco difference)", metavar="<path>"),
	make_option(c("-o", "--path_output"), type="character", default="./", help="Path to output directory (default is ./)", metavar="<path>"),
	make_option(c("-m", "--list_motif"), type="character", default=NULL, help="Comma separated list of motifs following IUPAC nucleotide code (e.g. GATC,CCWGG)", metavar="<motif1,motif2,...>"),
	make_option(c("-t", "--type_model"), type="character", default=NULL, help="Comma separated list of model to apply (nn: neural network, rf: random forest, or knn: k-nearest neighbor; e.g. nn,rf)", metavar="<models>"),
	make_option(c("-r", "--genome"), type="character", default=NULL, help="Path to reference genome (.fasta)", metavar="<path>"),
	make_option(c("-c", "--list_contig"), type="character", default=NULL, help="Comma separated list of contigs (e.g. contig_1,contig_3)", metavar="<contig_1,contig_3,...>"),
	make_option(c("--contigs_file"), type="character", default=NULL, help="Path to file with list of contigs (one per line)", metavar="<path>")
)

default_usage <- c("nanodisco characterize -p <nb_threads> -b <analysis_name> -d <path_difference> -o <path_output> -m <motif1,motif2,...> -t <type_model> -r <path_fasta>")
opt_parser <- OptionParser(usage=default_usage, option_list=option_list)
opt <- parse_args(opt_parser)

# Load parameters and functions
path_models <- "/home/nanodisco/models/"
path_script <- "/home/nanodisco/code/"
source(paste0(path_script,"analysis_functions.R"))
suppressMessages(load.libraries.characterize())

## Check input parameters
opt <- check.input.characterize(opt)

nb_threads <- opt$nb_threads
base_name <- opt$base_name
path_diff_data <- opt$path_diff_data
path_output <- paste0(gsub("/$","",opt$path_output),"/") # Make sure it's ending with a slash
list_motif <- str_split(opt$list_motif, ",", simplify=TRUE)[1,]
list_model <- str_split(opt$type_model, ",", simplify=TRUE)[1,]
genome <- opt$genome
list_contig <- opt$list_contig

print_message("Load supplied current differences")

# Load supplied current differences
difference_data <- readRDS(path_diff_data) # Read saved full dataset
motif_summary <- data.frame(motif=list_motif, mod_type=NA, mod_pos=NA, stringsAsFactors=FALSE)

print_message("Check current differences file version")

# Check current differences version
model_version <- check.model.version(difference_data)

if(!is.na(list_contig)){
	print_message("Select subset of contigs")

	difference_data <- subset(difference_data, contig %in% list_contig)
}

print_message("Determine motif signature center")

# Determine motif signature center; approximation of modification position
motif_center_summary <- find.signature.center(difference_data, motif_summary, genome, nb_threads, seq_params, iupac_nc, 6, 6, paste0(base_name,"_example"), path_output)
motif_center_summary$col_motif <- as.vector(iwanthue(nrow(motif_center_summary))) # Add color to motifs

# Classify motif with requested model(s)
stifle <- foreach(model_name=list_model) %do% {
	print_message(paste0("Classify motif(s) with ",model_name," model"))
	
	# Load classifier models
	model <- readRDS(file=paste0(path_models,"final_model_",model_version$basecaller,"_",model_version$version,"_",model_name,".RDS"))

	# Predict modification type and position
	classification_results <- classify.detected.motifs(difference_data, base_name, motif_center_summary, model, genome, 0, TRUE, iupac_nc, nb_threads)

	# Plot predictions as heatmap
	draw.classification.results(classification_results, paste0(base_name,"_",model_name,"_model"), path_output)

	# Write predictions as tsv
	write.best.classification.results(classification_results, paste0(base_name,"_",model_name,"_model"), path_output)

	return(NA)
}

print_message("Done")
