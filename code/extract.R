#!/usr/bin/env Rscript

# /home/nanodisco/code/extract.R -i dataset/EC_WGA -o analysis/preprocessed_subset -b EC_WGA -p 40 -c 5000 -s fa
# /home/nanodisco/code/extract.R -i dataset/EC_NAT -o analysis/preprocessed_subset -b EC_NAT -p 40 -c 5000 -s fa

load.libraries.extract <- function(){
	library("rhdf5")
	library("foreach")
	library("doMC")
	library("optparse")
	library("progress")
	library("Biostrings")
	library("stringr")
}

extract.basecall.version <- function(path_first_fast5, path_basecalling){
	f5_data <- h5readAttributes(path_first_fast5, path_basecalling) # Extract read data (fastq, move, and trace)
	
	if(length(f5_data)==0){
		# This is likely a live base called read with MinKNOW
		f5_data <- h5readAttributes(path_first_fast5, "/UniqueGlobalKey/tracking_id") # Extract read data (fastq, move, and trace)

		basecaller <- "Guppy_live"
		bc_version <- f5_data$guppy_version
	}else{
		# This is likely an "offline" base calling
		if(grepl("Albacore",f5_data$name)){
			basecaller <- "Albacore"
			bc_version <- f5_data$version
		}else if(grepl("Guppy",f5_data$name)){
			basecaller <- "Guppy"
			bc_version <- f5_data$version
		}else{
			basecaller <- "Unknown"
			bc_version <- f5_data$version
		}
	}

	return(data.frame(basecaller=basecaller, version=bc_version))
}

find.basecall.versions <- function(path_fast5_dir){
	# List fast5 files
	list_fast5_files <- list.files(path_fast5_dir, pattern="*.fast5", recursive=TRUE)

	# Extract only 1st path
	path_first_fast5 <- paste0(gsub("/$","",normalizePath(path_fast5_dir)),"/",head(list_fast5_files,1))
	f5_content <- h5ls(path_first_fast5)

	if(any(grepl("^/read_.*",f5_content$group, perl=TRUE)==TRUE)){
		# This is a multi-read fast5 file
		path_first_read <- f5_content$group[grepl("^/read_.*",f5_content$group, perl=TRUE)][1]
		
		list_basecalling <- gsub("^/read_.*/Analyses/(.*)/BaseCalled_template","\\1",subset(f5_content, grepl(path_first_read, group) & name=="Fastq")$group)
		available_versions <- foreach(basecall=list_basecalling, .combine=rbind) %do% {
			path_basecalling <- paste0(path_first_read,"/Analyses/",basecall)

			available_version <- extract.basecall.version(path_first_fast5, path_basecalling)
			available_version$basecall_group <- basecall

			return(available_version)
		}
	}else if(any(grepl("^/Analyses",f5_content$group, perl=TRUE)==TRUE)){
		# This is NOT a multi-read fast5 file
		list_basecalling <- gsub("^/Analyses/(.*)/BaseCalled_template","\\1",subset(f5_content, name=="Fastq")$group)
		available_versions <- foreach(basecall=list_basecalling, .combine=rbind) %do% {
			path_basecalling <- paste0("/Analyses/",basecall)

			available_version <- extract.basecall.version(path_first_fast5, path_basecalling)
			available_version$basecall_group <- basecall

			return(available_version)
		}
	}

	return(available_versions)
}

check.input.extract <- function(opt){
	# Check input directory
	if(is.null(opt$path_input)){
		cat("Parameter -i/--path_input is missing. Please provide the path to an input directory containing fast5.\n")
		quit(save="no", status=3)
	}else{
		# Check if input directory exist and contain fast5s
		if(dir.exists(opt$path_input)){
			nb_fast5 <- length(list.files(paste0(opt$path_input), pattern="*.fast5", recursive=TRUE))
			if(nb_fast5==0){
				cat(paste0("No fast5 files were found in input (",normalizePath(opt$path_input),").\n"))
				cat("Please check -i/--path_input.\n")
				quit(save="no", status=4)
			}
		}else{
			cat("-i/--path_input directory doesn't exist.\n")
			quit(save="no", status=4)
		}
	}
	# Check output directory
	if(is.null(opt$path_output)){
		cat("Parameter -o/--path_output is missing. Please provide the path to an output directory.\n")
		quit(save="no", status=3)
	}else{
		# Create output directory if necessary
		if(!dir.exists(opt$path_output)){
			dir.create(opt$path_output, recursive=TRUE)
		}
	}
	# Check base name for sequence output
	if(is.null(opt$base_name)){
		cat("Parameter -b/--base_name is missing. Please provide the sequence file output name (<base_name>.<seq_type>).\n")
		quit(save="no", status=3)
	}
	# Control requested output file type (fasta or fastq)
	if(opt$seq_type %in% c("fq","fastq")){
		opt$seq_type <- "fastq"
	}else if(opt$seq_type %in% c("fa","fasta")){
		opt$seq_type <- "fasta"
	}else{
		cat("Parameter -s/--seq_type is not recognized. Please select between fastq or fasta (short option fq or fa).\n")
		quit(save="no", status=4)
	}
	# Control requested basecalling version available
	if(opt$basecall_version %in% c("default")){
		opt$basecall_group <- "Basecall_1D_000" # Use the default location
	}else{
		if(!grepl(":",opt$basecall_version)){
			cat("Parameter --basecall_version is not recognized. Please provide the basecaller name and version (e.g. Guppy:3.2.4).\n")
			quit(save="no", status=4)
		}
		basecall_version <- strsplit(opt$basecall_version,":")[[1]]
		available_versions <- find.basecall.versions(opt$path_input)
		matching_version <- subset(available_versions, basecaller==basecall_version[1] & grepl(paste0("^",basecall_version[2]), version))
		
		if(nrow(matching_version)==1){
			opt$basecall_group <- matching_version$basecall_group
		}else{
			cat("Parameter --basecall_version don't match available basecalling. Please provide the basecaller name and version (e.g. Guppy:3.2.4).\n")
			cat("Available basecalling are displayed below.\n")
			print(available_versions)
			quit(save="no", status=4)
		}
	}
	
	return(opt)
}

print_message <- function(message){
	# Print message to terminal
	cat(paste0("[",Sys.time(),"] ",message,".\n"))
}

make.progress.bar <- function(chunk_list_fast5_files, nbCPU){
	pb <- progress_bar$new(format=" Processed fast5 [:bar] :percent eta: :eta (elapsed: :elapsedfull)", total=length(chunk_list_fast5_files)/nbCPU, show_after=0)

	return(pb)
}

initialize.progress.bar <- function(pb){
	# Forces progress bar to be displayed
	pb$tick(0)
}

progress.tracker <- function(pb, chunk_idx, nbCPU){
	if(chunk_idx %% nbCPU == 0){
		pb$tick()
	}
}

terminate.progress.bar <- function(pb, chunk_list_fast5_files, nbCPU){
	# Forces progress bar to complete before removed, not critical
	pb$tick(length(chunk_list_fast5_files)/nbCPU)
}

extract.sequence <- function(path_fast5, sample_name, path_output, nbCPU, chunk_size, seq_type, basecall_group){
	list_fast5_files <- list.files(paste0(path_fast5), pattern="*.fast5", recursive=TRUE)

	chunk_list_fast5_files <- split(list_fast5_files, ceiling(seq_along(list_fast5_files)/chunk_size))

	# Define which row to keep from the default fastq information present in fast5
	if(seq_type=="fastq"){
		row_to_keep <- "2:4"
	}else if(seq_type=="fasta"){
		row_to_keep <- "2"
	}

	pb <- make.progress.bar(chunk_list_fast5_files, nbCPU)
	print_message("Extract sequences from fast5")
	initialize.progress.bar(pb)
	registerDoMC(cores=nbCPU)
	sequences <- foreach(chunk_idx=seq(1,length(chunk_list_fast5_files))) %dopar% {
		fast5_files <- chunk_list_fast5_files[[chunk_idx]]

		subset_sequence <- foreach(f5_file=fast5_files) %do% {
			f5_file_path <- normalizePath(paste0(path_fast5, f5_file))
			f5_content <- h5ls(f5_file_path)
			if(any(grepl(paste0("^/read_.*/Analyses/",basecall_group,"/BaseCalled_template"),f5_content$group, perl=TRUE)==TRUE)){
				# This is a multi-read fast5 file
				list_fast5_reads <- unique(f5_content$group[grepl(paste0("^/read_.*/Analyses/",basecall_group,"/BaseCalled_template"),f5_content$group, perl=TRUE)])
				list_linear_sequence <- foreach(f5_read_data=list_fast5_reads) %do% { # , .combine=c
					f5_data <- h5read(f5_file_path, f5_read_data) # Extract read data (fastq, move, and trace)
					detail_fastq <- strsplit(f5_data$Fastq,"\n")[[1]] # Isolate fastq
					read_name <- gsub("^@",">",strsplit(detail_fastq[1]," ")[[1]][1]) # Edit read name to match indexed formating
					read_name_indexed <- paste0(read_name," ",gsub(".fast5","",basename(f5_file_path))," ",f5_file_path) # Assemble read name and index
					sequence_noname <- paste0(detail_fastq[eval(parse(text=row_to_keep))], collapse="\n") # Aggrerate sequence data without name (fasta or fastq)

					subset_linear_sequence <- paste0(read_name_indexed,"\n",sequence_noname,"\n") # Assemble read information as fasta or fastq

					return(subset_linear_sequence)
				}
				linear_sequence <- do.call(paste0,list_linear_sequence) # Put all fasta/fastq entries in same string

			}else if(any(grepl(paste0("^/Analyses/",basecall_group,"/BaseCalled_template"),f5_content$group, perl=TRUE)==TRUE)){
				# This is NOT a multi-read fast5 file
				f5_data <- h5read(f5_file_path,paste0("/Analyses/",basecall_group,"/BaseCalled_template")) # Extract read data (fastq, move, and trace)
				detail_fastq <- strsplit(f5_data$Fastq,"\n")[[1]] # Isolate fastq
				read_name <- gsub("^@",">",strsplit(detail_fastq[1]," ")[[1]][1]) # Edit read name to match indexed formating
				read_name_indexed <- paste0(read_name," ",gsub(".fast5","",basename(f5_file_path))," ",f5_file_path) # Assemble read name and index
				sequence_noname <- paste0(detail_fastq[eval(parse(text=row_to_keep))], collapse="\n") # Aggrerate sequence data without name (fasta or fastq)

				linear_sequence <- paste0(read_name_indexed,"\n",sequence_noname) # Assemble read information as fasta or fastq
			}else{
				linear_sequence <- "uncalled"
			}

			return(linear_sequence)
		}

		progress.tracker(pb, chunk_idx, nbCPU)

		return(do.call(rbind,subset_sequence))
	}
	registerDoSEQ()
	terminate.progress.bar(pb, chunk_list_fast5_files, nbCPU)
	sequences <- do.call(rbind,sequences)

	# Remove uncalled reads and print warning
	isnot_basecalled <- grepl("uncalled",sequences,fixed=T)
	if(any(isnot_basecalled)){
		warning(paste0(sum(isnot_basecalled==TRUE)," reads weren't basecalled."))
		sequences <- sequences[!isnot_basecalled]
	}

	write.table(sequences, file=paste0(path_output,sample_name,".",seq_type), quote=F, row.names=F, col.names=F)
}

suppressMessages(library(optparse))

option_list <- list(
	make_option(c("-i", "--path_input"), type="character", default=NULL, help="Path to input directory (recursively look for *.fast5 files)", metavar="<path>"),
	make_option(c("-o", "--path_output"), type="character", default=NULL, help="Path to output directory (keep trailing slash)", metavar="<path>"),
	make_option(c("-b", "--base_name"), type="character", default=NULL, help="Sequence file output name (<base_name>.<seq_type>)", metavar="<name>"),
	make_option(c("-p", "--nb_threads"), type="integer", default=1, help="Number of threads (default is 1)", metavar="<integer>"),
	make_option(c("-s", "--seq_type"), type="character", default=NULL, help="Type of sequence to extract, fastq/fa or fastq/fq", metavar="<fq/fa/fastq/fasta>"),
	make_option(c("-c", "--nb_chunks"), type="integer", default=40, help="Number of reads per chunks (default is 40; if single-fast5s then best nb_chunks >= nb_threads, if multi-fast5 then best nb_chunks == nb_threads)", metavar="<integer>"),
	make_option(c("--basecall_version"), type="character", default="default", help="Basecalling version (when multiple ones available)", metavar="<basecaller:version>")
)

default_usage <- c("%prog -i <fast5_directory> -o <path_output> -b <sample_name> -p <nb_threads> -c <nb_reads_per_chunks> -s <fa|fq>")
opt_parser <- OptionParser(usage=default_usage, option_list=option_list)
opt <- parse_args(opt_parser)

suppressMessages(load.libraries.extract())

opt <- check.input.extract(opt)

path_input <- paste0(gsub("/$","",normalizePath(opt$path_input)),"/") # Deal with relative path, stay usable afetr setwd()
base_name <- opt$base_name
path_output <- paste0(gsub("/$","",normalizePath(opt$path_output)),"/") # Deal with relative path, stay usable afetr setwd()
nb_threads <- opt$nb_threads
nb_chunks <- opt$nb_chunks
seq_type <- opt$seq_type
basecall_group <- opt$basecall_group

extract.sequence(path_input, base_name, path_output, nb_threads, nb_chunks, seq_type, basecall_group)
