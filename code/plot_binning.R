#!/usr/bin/env Rscript

# Example
# nanodisco plot_binning -r /home/nanodisco/reference/metagenome.fasta -u /home/nanodisco/analysis/methylation_binning_MGM1.RDS -b MGM1 -o ./ [-a /home/nanodisco/reference/motif_binning_annotation.RDS --MGEs_file /home/nanodisco/dataset/list_MGE_contigs.txt]

suppressMessages(library(optparse))

# Parsing arguments
option_list <- list(
	make_option(c("-r", "--metagenome"), type="character", default=NULL, help="Path to reference metagenome (.fasta)", metavar="<path>"),
	make_option(c("-u", "--path_methylation_binning"), type="character", default=NULL, help="Path to methylation binning file (*.RDS produced from nanodisco binning)", metavar="<path>"),
	make_option(c("-b", "--base_name"), type="character", default="results", help="Base name for outputing results (e.g. Ecoli_K12; default is 'results')", metavar="<character>"),
	make_option(c("-o", "--path_output"), type="character", default="./", help="Path to output directory (default is ./)", metavar="<path>"),
	make_option(c("-a", "--path_annotation"), type="character", default=NULL, help="Path to contig annotation", metavar="<path>"),
	make_option(c("-c", "--list_MGE_contig"), type="character", default=NULL, help="Comma separated list of MGE contigs (e.g. contig_1,contig_3)", metavar="<contig_1,contig_3,...>"),
	make_option(c("--MGEs_file"), type="character", default=NULL, help="Path to file with list of MGE contigs (one per line)", metavar="<path>"),
	make_option(c("--xlim"), type="character", default=NULL, help="Optional x-axis zooming (e.g. -5:10)", metavar="<x1:x2>"),
	make_option(c("--ylim"), type="character", default=NULL, help="Optional y-axis zooming (e.g. -10:9)", metavar="<y1:y2>"),
	make_option(c("--min_contig_len"), type="integer", default=25000, help="Minimum length for plotting contigs (default is 25000 bp)", metavar="<integer>"),
	make_option(c("--split_fasta"), type="character", default="no", help="Split reference metagenome into binned fasta ('yes' split from annotation, 'default'|'<integer,integer>' split from dbscan cluster analysis", metavar="<no|yes|default|integer,integer>")
)

default_usage <- c("nanodisco plot_binning -r <path_fasta> --mb <path_methylation_binning> -b <analysis_name> -o <path_output> [+ advanced parameters]")
opt_parser <- OptionParser(usage=default_usage, option_list=option_list)
opt <- parse_args(opt_parser);

# Load parameters and functions
source("/home/nanodisco/code/metagenome_analysis_functions.R")
suppressMessages(load.libraries.metagenome())

## Check input parameters
opt <- check.input.plot.binning(opt)

metagenome <- opt$metagenome
path_methylation_binning <- opt$path_methylation_binning
base_name <- opt$base_name
path_output <- paste0(gsub("/$","",opt$path_output),"/") # Make sure it's ending with a slash
path_annotation <- opt$path_annotation
type_annotation <- opt$type_annotation
list_MGE_contig <- str_split(opt$list_MGE_contig, ",", simplify=TRUE)[1,]
new_xlim <- str_split(opt$xlim, ":", simplify=TRUE)[1,] # NA by default
new_ylim <- str_split(opt$ylim, ":", simplify=TRUE)[1,] # NA by default
min_contig_len <- opt$min_contig_len
split_fasta <- opt$split_fasta
param_dbscan <- opt$param_dbscan

print_message("Prepare default metagenome annotation")

# Load metagenome information
metagenome_annotation <- readDNAStringSet(metagenome)
metagenome_annotation <- data.frame(contig=names(metagenome_annotation), length=width(metagenome_annotation), id=NA, stringsAsFactors=TRUE) # Simple metagenome contig annotation

# Load methylation binning results
methylation_binning <- readRDS(path_methylation_binning)

# Load additional annotation (optional)
if(is.null(path_annotation)){
	if(split_fasta){
		print_message("Detection of potential bins")

		methylation_binning_annotated <- find.tsne.clusters(methylation_binning, param_dbscan[["set_eps"]], param_dbscan[["set_minPts"]])
		binning_annotation <- data.frame(contig=methylation_binning_annotated$contig, id=methylation_binning_annotated$id, stringsAsFactors=TRUE)
	}else{
		# No annotation provided
		binning_annotation <- metagenome_annotation
		binning_annotation$id <- as.factor("No annotation")
	}
}else{
	print_message("Load additional annotation")

	if(type_annotation=="is_rds"){
		binning_annotation <- readRDS(path_annotation)
		colnames(binning_annotation) <- c("contig","id")
	}else if(type_annotation=="is_txt"){
		binning_annotation <- read.table(path_annotation, stringsAsFactors=FALSE)
		colnames(binning_annotation) <- c("contig","id")
	}
}
motif_binning_legend <- generate.color.palette(binning_annotation)

print_message("Plot binning")

# Plot binning
gp_motif_binning <- plot.tsne.motifs.score(methylation_binning, binning_annotation, base_name, motif_binning_legend, list_MGE_contig, new_xlim, new_ylim, min_contig_len, path_output)

if(split_fasta){
	print_message("Generate binned fasta files")

	write.binned.fasta(metagenome, metagenome_annotation, binning_annotation, base_name, path_output)
}

print_message("Done")

