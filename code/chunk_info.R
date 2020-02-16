#!/usr/bin/env Rscript

# Example
# nanodisco chunk_info -r /home/nanodisco/reference/Ecoli_K12_MG1655_ATCC47076.fasta
# nanodisco chunk_info -r /home/nanodisco/reference/Ecoli_K12_MG1655_ATCC47076.fasta -t CP014225.1:25050-84950

suppressMessages(library(optparse))

option_list <- list(
	make_option(c("-r", "--genome"), type="character", default=NULL, help="Path to reference genome (.fasta)", metavar="<path>"),
	make_option(c("-s", "--chunk_size"), type="integer", default=5000, help="Size of chunk to use (default is 5000)", metavar="<integer>"),
	make_option(c("-t", "--target_region"), type="character", default=NULL, help="Specify a genomic region whose chunks need to be processed (e.g. chr1:2500-85000)", metavar="<contig:start-end>")
)

default_usage <- c("nanodisco chunk_info -r <path_fasta> [-t <contig:start-end> -s <integer>]")
opt_parser <- OptionParser(usage=default_usage, option_list=option_list)
opt <- parse_args(opt_parser)

# Load parameters and functions
path_script <- "/home/nanodisco/code/"
source(paste0(path_script,"difference_functions.R"))
suppressMessages(load.libraries.chunk_info())

## Check input parameters
opt <- check.input.chunk_info(opt)

genome <- opt$genome
chunk_size <- opt$chunk_size
target_region <- opt$target_region

# Define chunks coordinate
chunk_information <- generate.genome.chunks.information(genome, chunk_size)

if(is.null(target_region)){
	response <- paste0("Number of chunks: ",nrow(chunk_information),"\n")
}else{
	coordinate <- as.data.frame(matrix(str_split(gsub("^(.*):([0-9]+)-([0-9]+)$","\\1 \\2 \\3",target_region), " ", simplify=TRUE)[1,], ncol=3), stringsAsFactors=FALSE) # Dup
	colnames(coordinate) <- c("contig","start","end") # Dup
	coordinate$start <- as.numeric(coordinate$start) # Dup
	coordinate$end <- as.numeric(coordinate$end) # Dup

	first_chunk <- subset(chunk_information, contig_name==coordinate$contig & start <= coordinate$start & start > coordinate$start - chunk_size)$chunk_id
	last_chunk <- subset(chunk_information, contig_name==coordinate$contig & end >= coordinate$end & end < coordinate$end + chunk_size)$chunk_id

	if(length(first_chunk)==0 | length(last_chunk)==0){
		cat(paste0("None of the chunks match provided genomic region (",target_region,").\n"))
		cat("Please check -t/--target_region parameter. Specify a genomic region as <contig:start-end> whose chunk(s) need to be processed (e.g. chr1:2500-85000).\n")
		quit(save="no", status=4)
	}

	response <- paste0("First chunk: ",first_chunk,"\nLast chunk: ",last_chunk,"\n")
}

# Print response
cat(response)
