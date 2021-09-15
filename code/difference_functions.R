load.libraries <- function(){
	library(optparse)
	library(data.table)
	library(doMC)
	library(foreach)
	library(rhdf5)
	library(GenomicRanges)
	library(GenomicAlignments) # Load GenomicRanges Biostrings Rsamtools
	library(Rsamtools)
	library(Biostrings) # In Rsamtools
	library(plyr)
	library(dplyr)
	library(tidyr)
	library(Matrix)
	library(MASS) # For rlm; Masked select
	library(stringr)
}

load.libraries.chunk_info <- function(){
	library(optparse)
	library(foreach)
	library(stringr)
}

print_message <- function(message){
        # Print message to terminal
        cat(paste0("[",Sys.time(),"] ",message,".\n"))
}

check.input.chunk_info <- function(opt){
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
	# Check if genomic region is not recognized
	if(!is.null(opt$target_region)){
		if(!grepl("^(.*):([0-9]+)-([0-9]+)$", opt$target_region)){
			cat(paste0("Genomic region is not recognized (",opt$target_region,").\n"))
			cat("Please check -t/--target_region parameter. Specify a genomic region as <contig:start-end> whose chunks need to be processed (e.g. chr1:2500-85000).\n")
			quit(save="no", status=4)
		}else{
			coordinate <- as.data.frame(matrix(str_split(gsub("^(.*):([0-9]+)-([0-9]+)$","\\1 \\2 \\3",opt$target_region), " ", simplify=TRUE)[1,], ncol=3), stringsAsFactors=FALSE) # Dup
			colnames(coordinate) <- c("contig","start","end") # Dup
			coordinate$start <- as.numeric(coordinate$start) # Dup
			coordinate$end <- as.numeric(coordinate$end) # Dup

			if(coordinate$start > coordinate$end){
				cat(paste0("Genomic region is not recognized, start and end seems inverted (",opt$target_region,").\n"))
				cat("Please check -t/--target_region parameter. Specify a genomic region as <contig:start-end> whose chunks need to be processed (e.g. chr1:2500-85000).\n")
				quit(save="no", status=4)

			}
		}
	}

	return(opt)
}

generate.reverse.complement.genome <- function(genome, reverse_complement_genome){
	rev_comp_genome_sequences <- reverseComplement(readDNAStringSet(genome))
	writeXStringSet(rev_comp_genome_sequences, file=reverse_complement_genome, width=80)
}

check.reference <- function(genome){
	if(!file.exists(genome)){
		stop(paste0("Reference genome (",genome,") is missing."))
	}else{
		reverse_complement_genome <- gsub(".fasta",".rev_comp.fasta",genome)
		if(!file.exists(reverse_complement_genome)){
			print_message(paste0("Generating reverse complement for ",genome))
			generate.reverse.complement.genome(genome, reverse_complement_genome)
		}
		# TODO add check for special characters in references names
	}
}

extract.basecall.version <- function(path_first_fast5, path_basecalling, path_first_read=NA){
	f5_data <- h5readAttributes(path_first_fast5, path_basecalling) # Extract read data (fastq, move, and trace)
	
	if(length(f5_data)==0){
		# This is likely a live base called read with MinKNOW
		f5_data <- h5readAttributes(path_first_fast5, "/UniqueGlobalKey/tracking_id") # Extract read data (fastq, move, and trace)
		if(is.null(f5_data$guppy_version)){
			basecaller <- "Unknown"
			bc_version <- "0.0.0"
		}else{
			basecaller <- "Guppy_live"
			bc_version <- f5_data$guppy_version
		}
	}else{
		# This is likely an "offline" base calling
		if(grepl("Albacore",f5_data$name)){
			basecaller <- "Albacore"
			bc_version <- f5_data$version
		}else if(grepl("Guppy",f5_data$name)){
			basecaller <- "Guppy"
			bc_version <- f5_data$version
		}else if(grepl("MinKNOW-Live-Basecalling",f5_data$name)){
			if(is.na(path_first_read)){ # If not multi-read fast5
				f5_data <- h5readAttributes(path_first_fast5, "/UniqueGlobalKey/tracking_id") # Extract read data (fastq, move, and trace)
			}else{
				f5_data <- h5readAttributes(path_first_fast5, paste0(path_first_read, "/tracking_id")) # Extract read data (fastq, move, and trace)
			}

			basecaller <- "Guppy_live"
			bc_version <- f5_data$guppy_version
		}else{
			basecaller <- "Unknown"
			bc_version <- f5_data$version
		}
	}

	return(data.frame(basecaller=basecaller, version=bc_version, stringsAsFactors=FALSE))
}

find.basecall.versions <- function(path_indexed_fasta){
	# Read fasta file containing path to fast5 in header
	read_name_nat <- names(readDNAStringSet(path_indexed_fasta))

	# Extract only 1st path
	path_first_fast5 <- head(str_split(read_name_nat, " ", simplify=TRUE)[,3], 1)
	f5_content <- h5ls(path_first_fast5)

	if(any(grepl("^/read_.*",f5_content$group, perl=TRUE)==TRUE)){
		# This is a multi-read fast5 file
		path_first_read <- f5_content$group[grepl("^/read_.*",f5_content$group, perl=TRUE)][1]
		
		list_basecalling <- gsub("^/read_.*/Analyses/(.*)/BaseCalled_template","\\1",subset(f5_content, grepl(path_first_read, group) & name=="Fastq")$group)
		available_versions <- foreach(basecall=list_basecalling, .combine=rbind) %do% {
			path_basecalling <- paste0(path_first_read,"/Analyses/",basecall)

			available_version <- extract.basecall.version(path_first_fast5, path_basecalling, path_first_read)
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
	available_versions$simplified <- paste0(available_versions$basecaller, "_", available_versions$version)

	return(available_versions)
}

check.basecall.version <- function(path_input, sample_name_nat, sample_name_wga, basecall_version){
	# Get basecaller version from fast5 files
	available_versions_wga <- find.basecall.versions(paste0(path_input,sample_name_wga,".fasta"))
	available_versions_nat <- find.basecall.versions(paste0(path_input,sample_name_nat,".fasta"))

	# Compare basecalling version
	common_basecalling <- intersect(available_versions_wga$simplified, available_versions_nat$simplified)
	if(length(common_basecalling)==0){
		cat("Base calling was performed with different versions.\n")
		cat(paste0("Native reads:\n"))
		print(available_versions_nat)
		cat(paste0("WGA reads:\n"))
		print(available_versions_wga)
		cat("We highly recommend to use the same base caller version as signal properties can differ.\n")
	}else if(length(common_basecalling)>1){
		available_versions_wga <- subset(available_versions_wga, simplified %in% common_basecalling)
		available_versions_nat <- subset(available_versions_nat, simplified %in% common_basecalling)
		if(basecall_version=="default"){
			default_basecalling_index <- match(common_basecalling[1],available_versions_nat$simplified)
			default_basecalling_version <- c(available_versions_nat$basecaller[default_basecalling_index], available_versions_nat$version[default_basecalling_index])
			cat("Multiple base calling versions are available.\n")
			cat(paste0("The first matching version will be considered for tracking purpose (i.e. ",paste0(default_basecalling_version, collapse=":"),").\n"))
			cat("The desired base calling version can be set with --basecall_version in \"nanodisco preprocess\" and \"nanodisco difference\".\n")

			return(c(available_versions_nat$basecaller[default_basecalling_index], available_versions_nat$version[default_basecalling_index]))
		}else{
			if(!grepl(":",basecall_version)){
				cat("Parameter --basecall_version is not recognized. Please provide the basecaller name and version (e.g. Guppy:3.2.4).\n")
				quit(save="no", status=4)
			}

			basecall_version <- strsplit(basecall_version,":")[[1]]
			matching_version <- subset(available_versions_nat, basecaller==basecall_version[1] & grepl(paste0("^",basecall_version[2]), version))
			
			if(nrow(matching_version)==1){
				return(c(matching_version$basecaller, matching_version$version))
			}else{
				cat("Parameter --basecall_version don't match available basecalling. Please provide the basecaller name and version (e.g. Guppy:3.2.4).\n")
				cat("Available basecalling are displayed below.\n")
				print(subset(available_versions_nat, select=c("basecaller","version")))
				quit(save="no", status=4)
			}
		}
	}else if(length(common_basecalling)==1){
		default_basecalling_index <- match(common_basecalling[1],available_versions_nat$simplified)

		return(c(available_versions_nat$basecaller[default_basecalling_index], available_versions_nat$version[default_basecalling_index]))
	}
}

# tmp
error.handling <- function(error_id){
	if(error_id==1){
		# TODO Handle by R
	}else if(error_id==2){
		print_message("     A reference fasta file is missing") # TODO handle by R chartr("ATGC","TACG",s)
	}
}

# Processing
clean.temporary.files <- function(path_output, idx_chunk){
	print_message("  Remove temporary files")
	# Remove temp dirs
	# list_directories <- list.dirs()
	# unlink(list_directories[grep(paste0("tmp.",idx_chunk), list_directories)], recursive=TRUE)
	
	# Remove temp files
	list_tmp_files <- list.files(path=path_output, pattern=paste0("tmp[.]",idx_chunk,"[.]"))
	if(length(list_tmp_files)>0){ # Attempt to remove path_output if no matches
		stifle <- file.remove(paste0(path_output,list_tmp_files))
	}
}

relative.path <- function(x, y){
	if(x==y){
		# Nothing to change
		if(substr(x, 1, 1)=="/"){
			return("")
		}else{
			return("./")
		}
	}

	split_x <- strsplit(x,"/")[[1]]
	split_y <- strsplit(y,"/")[[1]]

	len_x <- length(split_x)
	len_y <- length(split_y)
	max_len <- min(len_x, len_y)
	path_diff <- split_x[seq(1,max_len)]==split_y[seq(1,max_len)] # Compare path with same level to find where it's differ
	path_diff <- c(path_diff, rep(FALSE,max(len_x - max_len, len_y - max_len))) # Add differences for uncomparable dir
	step_back <- rep("../",table(path_diff)[["FALSE"]])
	if(is.na(split_x[!path_diff])){ # x shorter than y
		path_x_relative_to_y <- step_back
	}else{
		path_branch <- paste(split_x[!path_diff], collapse="/")
		path_x_relative_to_y <- paste0(paste0(step_back, collapse=""), path_branch, "/", collapse="")
	}

	return(path_x_relative_to_y)
}

generate.genome.chunks.information <- function(genome, chunk_size){
	genome_index <- read.table(gsub(".fasta",".fasta.fai",genome),stringsAsFactors=FALSE)
	list_contig_names <- genome_index[,1]

	# Define chunks id across contigs
	starting_chunk <- 1
	nb_prev_chunks <- 0
	contigs_chunks <- foreach(contig=list_contig_names, .combine=rbind) %do% {
		contig_size <- genome_index[which(genome_index[,1]==contig),2]
		nb_chunks_contig <- ceiling(contig_size/chunk_size)

		starting_chunk <- starting_chunk + nb_prev_chunks

		list_contig_chunks_id <- seq(starting_chunk, starting_chunk + nb_chunks_contig - 1)
		list_contig_chunks_starts <- seq(1, contig_size, chunk_size)
		if(contig_size>chunk_size){
			if(contig_size %% chunk_size==0){ # Only full chunks
				list_contig_chunks_ends <- seq(chunk_size, contig_size, chunk_size)
			}else{
				list_contig_chunks_ends <- c(seq(chunk_size, contig_size, chunk_size), contig_size)
			}
		}else{
			list_contig_chunks_ends <- contig_size
		}
		nb_prev_chunks <- nb_chunks_contig

		contig_chunks <- data.frame(contig_name=as.factor(contig), start=list_contig_chunks_starts, end=list_contig_chunks_ends, chunk_id=list_contig_chunks_id)

		return(contig_chunks)
	}

	return(contigs_chunks)
}

generate.detail.chunk.mapping <- function(chunk_fast5_info){
	detail_chunk_mapping <- chunk_fast5_info %>%
		dplyr::rename(start=pos) %>%
		mutate(end=start + cigarWidthAlongReferenceSpace(cigar) - 1) %>%
		dplyr::select(-cigar)
	detail_chunk_mapping <- droplevels(detail_chunk_mapping)

	return(detail_chunk_mapping)
}

compute.chunk.coverage <- function(detail_chunk_mapping){
	gr_mapping <- GRanges(
		seqnames=detail_chunk_mapping$rname,
		ranges=IRanges(detail_chunk_mapping$start,detail_chunk_mapping$end),
		strand=detail_chunk_mapping$strand
	)
	mapping_coverage <- foreach(curr_strand=as.character(unique(detail_chunk_mapping$strand)), .combine=rbind) %do% {
		strand_mapping_range <- subset(detail_chunk_mapping, strand==curr_strand)
		strand_mapping_range <- range(strand_mapping_range$start,strand_mapping_range$end)
		strand_mapping_coverage <- coverage(gr_mapping[strand(gr_mapping)==curr_strand])
		strand_mapping_coverage <- data.frame(
			rname=as.character(unique(detail_chunk_mapping$rname)),
			pos=seq(1,strand_mapping_range[2]),
			strand=curr_strand,
			coverage=sapply(strand_mapping_coverage, as.numeric)[,1]
		)

		return(strand_mapping_coverage)
	}

	return(mapping_coverage)
}

# For debug only
plot.chunk.coverage <- function(detail_chunk_mapping, cov_th, chunk_start=NA, chunk_end=NA){
	mapping_coverage <- compute.chunk.coverage(detail_chunk_mapping)

	gp <- ggplot(mapping_coverage) +
		geom_point(aes(x=pos, y=coverage, col=strand)) +
		geom_hline(yintercept=cov_th)

	if(!is.na(chunk_start) & !is.na(chunk_end)){
		gp <- gp +
			coord_cartesian(xlim=c(chunk_start, chunk_end))
	}else{
		subset_mapping_coverage <- subset(mapping_coverage, coverage > 0)
		range_mapping_coverage <- range(subset_mapping_coverage$pos)
		gp <- gp +
			coord_cartesian(xlim=c(range_mapping_coverage[1], range_mapping_coverage[2]))

	}
	print(gp)
}

prepare.index <- function(sample_name, path_input, path_output, genome, chunk_size, list_chunks){
	print_message(paste0("  Prepare index for ",sample_name,""))

	#Define subsets of read by genomic mapping position
	print_message(paste0("  Extract read mapped on chunks for ",sample_name,""))
	bamFile <- paste0(path_input,sample_name,".sorted.bam")

	contigs_chunks <- generate.genome.chunks.information(genome, chunk_size) # Can be long for metagenome; TODO move outside to do only once

	# Remove out of bound chunk indexes
	list_chunks <- list_chunks[list_chunks %in% seq(1,max(contigs_chunks$chunk_id))]

	chunks_fast5_info <- foreach(chunk=list_chunks, .combine=rbind) %do% {
		chunk_contig_name <- as.character(subset(contigs_chunks, chunk_id==chunk)$contig_name) # as.character not necessary
		chunk_start <- subset(contigs_chunks, chunk_id==chunk)$start
		chunk_end <- subset(contigs_chunks, chunk_id==chunk)$end

		# Retrieve mapping information on chunk
		range <- GRanges(seqnames=chunk_contig_name, ranges=IRanges(chunk_start,chunk_end))
		param <- ScanBamParam(which=range, what=c("qname","rname","pos","strand","cigar","qwidth"))
		bam <- scanBam(bamFile, param=param) # Can fail if not correctly indexed
		chunk_fast5_info <- as.data.frame(bam[[1]]) # Fail if no reads mapped to chunk

		if(nrow(chunk_fast5_info)>0){ # If mapped reads
			cov_th <- 1000
			detail_chunk_mapping <- generate.detail.chunk.mapping(chunk_fast5_info)
			extreme_coverage_region <- compute.chunk.coverage(detail_chunk_mapping) %>%
				filter(coverage>cov_th)
			if(nrow(extreme_coverage_region)>0){
				extreme_coverage_region <- extreme_coverage_region %>%
					group_by(rname, strand) %>% # TODO warning when empty
					arrange(rname, strand, pos) %>%
					mutate(step_prev=pos - lag(pos), step_next=lead(pos) - pos) %>%
					mutate(step_prev=ifelse(is.na(step_prev),2,step_prev)) %>% # Replace NA to filter
					mutate(step_next=ifelse(is.na(step_next),2,step_next)) %>% # Replace NA to filter
					filter(xor(step_prev!=1, step_next!=1)) # | would keep single pos; xor because single pos high coverage is rare, likely only when avg. cov at 1000 (and not a problem).
			}
			nb_region <- nrow(extreme_coverage_region)/2

			if(nb_region>0){
				# Downsample reads from repeat-like mapping (regional extreme coverage, ~ duplicate reads)
				full_region_reads <- foreach(region=seq(1,nb_region), .combine=c) %do% {
					region_rname <- as.character(extreme_coverage_region$rname[region*2])
					region_start <- extreme_coverage_region$pos[region*2-1]
					region_end <- extreme_coverage_region$pos[region*2]
					region_strand <- as.character(extreme_coverage_region$strand[region*2])
					region_length <- region_end - region_start - 1

					# Contig information for region
					subset_contigs_chunks <- subset(contigs_chunks, contig_name==region_rname)
					contig_start <- min(subset_contigs_chunks$start)
					contig_end <- max(subset_contigs_chunks$end)

					if(region_length<200){ # Avoid spurious regions
						write(paste0("	No regional downsampling for ",sample_name," chunk #",chunk,": region too short (",region_rname,":",region_start,"-",region_end,",",region_strand,"; ",region_length," bp)."), stderr()) # Write message in stderr
						
						region_reads_qname <- ""
					}else if(region_length>10000){ # Avoid downsampling when region longer then 10kb; however small plasmid could be missed
						write(paste0("	No regional downsampling for ",sample_name," chunk #",chunk,": region too long (",region_rname,":",region_start,"-",region_end,",",region_strand,"; ",region_length," bp)."), stderr()) # Write message in stderr

						region_reads_qname <- ""
					}else if(region_start==chunk_start & region_end==chunk_end){ # If downsampling cover the whole chunk
						strand_chunk_fast5_info <- subset(chunk_fast5_info, strand==region_strand)
						not_strand_chunk_fast5_info <- subset(chunk_fast5_info, strand!=region_strand)
						nb_reads_downsample <- nrow(strand_chunk_fast5_info) - region_length
						set.seed(101)
						strand_chunk_fast5_info <- strand_chunk_fast5_info %>%
							sample_n(nb_reads_downsample, replace=FALSE)

						region_reads_qname <- as.character(unique(strand_chunk_fast5_info$qname))
						write(paste0("	Random downsampling: ",length(region_reads_qname)," reads from ",sample_name," chunk #",chunk," (",region_rname,":",region_start,"-",region_end,",",region_strand,"; ",region_length," bp)."), stderr()) # Write message in stderr
					}else if(region_start==contig_start & region_end==contig_end){ # If downsampling cover the whole chromosome
						strand_chunk_fast5_info <- subset(chunk_fast5_info, strand==region_strand)
						not_strand_chunk_fast5_info <- subset(chunk_fast5_info, strand!=region_strand)
						nb_reads_downsample <- nrow(strand_chunk_fast5_info) - region_length
						set.seed(101)
						strand_chunk_fast5_info <- strand_chunk_fast5_info %>%
							sample_n(nb_reads_downsample, replace=FALSE)

						region_reads_qname <- as.character(unique(strand_chunk_fast5_info$qname))
						write(paste0("	Random contig downsampling: ",length(region_reads_qname)," reads from ",sample_name," chunk #",chunk," (",region_rname,":",region_start,"-",region_end,",",region_strand,"; ",region_length," bp)."), stderr()) # Write message in stderr
					}else{
						detail_chunk_mapping <- generate.detail.chunk.mapping(chunk_fast5_info)

						region_reads <- subset(detail_chunk_mapping, start >= region_start & end <= region_end & strand==region_strand)
						region_reads_qname <- as.character(unique(region_reads$qname))

						write(paste0("	Regional downsampling: ",length(region_reads_qname)," reads from ",sample_name," chunk #",chunk," (",region_rname,":",region_start,"-",region_end,",",region_strand,"; ",region_length," bp)."), stderr()) # Write message in stderr
					}
					
					return(region_reads_qname)
				}
				# Downsample region by removing reads
				chunk_fast5_info <- subset(chunk_fast5_info, !qname %in% unique(full_region_reads))

				detail_chunk_mapping <- generate.detail.chunk.mapping(chunk_fast5_info)
				extreme_coverage_region <- compute.chunk.coverage(detail_chunk_mapping) %>%
					filter(coverage>cov_th)
				if(nrow(extreme_coverage_region)>0){
					extreme_coverage_region <- extreme_coverage_region %>%
						group_by(rname, strand) %>% # TODO warning when empty
						arrange(rname, strand, pos) %>%
						mutate(step_prev=pos - lag(pos), step_next=lead(pos) - pos) %>%
						mutate(step_prev=ifelse(is.na(step_prev),2,step_prev)) %>% # Replace NA to filter
						mutate(step_next=ifelse(is.na(step_next),2,step_next)) %>% # Replace NA to filter
						filter(xor(step_prev!=1, step_next!=1)) # | would keep single pos; xor because single pos high coverage is rare, likely only when avg. cov at 1000 (and not a problem).
				}
				nb_region <- nrow(extreme_coverage_region)/2

				# Downsample reads (localized extreme coverage)
				if(nb_region>0){
					complete_list_reads_qname <- ""
					while(nb_region>0){
						region <- 1 # Always process 1st region until sorted out
						region_rname <- as.character(extreme_coverage_region$rname[region*2])
						region_start <- extreme_coverage_region$pos[region*2-1]
						region_end <- extreme_coverage_region$pos[region*2]
						region_strand <- as.character(extreme_coverage_region$strand[region*2])
						region_length <- region_end - region_start - 1

						# Compute mapping start position to identify region with extreme coverage
						summary_chunk_fast5_info <- chunk_fast5_info %>%
							group_by(rname, pos, strand) %>%
							summarize(nb_start=n(), .groups="drop_last") %>%
							filter(rname==region_rname & strand==region_strand & pos >= region_start & pos <= region_end) %>%
							arrange(desc(nb_start), pos)

						if(max(summary_chunk_fast5_info$nb_start)>5){
							nb_start_th <- 5
						}else{
							nb_start_th <- max(summary_chunk_fast5_info$nb_start) - 1
						}

						# Select position with more read start to downsample
						if(nb_start_th==0){ # Remove progressively, maximum 5 by 5.
							pos_to_downsample <- summary_chunk_fast5_info$pos[seq(1,min(10,length(summary_chunk_fast5_info$pos)))]
							region_reads <- chunk_fast5_info %>%
								group_by(rname, pos, strand) %>%
								filter(rname == summary_chunk_fast5_info$rname[1] & strand == summary_chunk_fast5_info$strand[1]) %>%
								filter(pos %in% pos_to_downsample)
						}else{
							region_reads <- chunk_fast5_info %>%
								group_by(rname, pos, strand) %>%
								mutate(nb_start=n()) %>%
								filter(rname == region_rname & strand == region_strand & pos >= region_start & pos <= region_end) %>%
								filter(nb_start > nb_start_th)
						}
						region_reads_qname <- as.character(unique(region_reads$qname))
						complete_list_reads_qname <- c(complete_list_reads_qname, region_reads_qname) # Store full list for statistics

						# Downsample region by removing reads
						chunk_fast5_info <- subset(chunk_fast5_info, !qname %in% unique(region_reads_qname))

						# Update coverage
						detail_chunk_mapping <- generate.detail.chunk.mapping(chunk_fast5_info)
						extreme_coverage_region <- compute.chunk.coverage(detail_chunk_mapping) %>%
							filter(coverage>cov_th)
						if(nrow(extreme_coverage_region)>0){
							extreme_coverage_region <- extreme_coverage_region %>%
								group_by(rname, strand) %>% # TODO warning when empty
								arrange(rname, strand, pos) %>%
								mutate(step_prev=pos - lag(pos), step_next=lead(pos) - pos) %>%
								mutate(step_prev=ifelse(is.na(step_prev),2,step_prev)) %>% # Replace NA to filter
								mutate(step_next=ifelse(is.na(step_next),2,step_next)) %>% # Replace NA to filter
								filter(xor(step_prev!=1, step_next!=1)) # | would keep single pos; xor because single pos high coverage is rare, likely only when avg. cov at 1000 (and not a problem).
						}
						nb_region <- nrow(extreme_coverage_region)/2
					}
					
					write(paste0("	Localized downsampling: ",length(complete_list_reads_qname)," reads from ",sample_name," chunk #",chunk,"."), stderr()) # Write message in stderr
				}
			}else{
				# write(paste0("	No downsampling for ",sample_name," chunk #",chunk,"."), stderr()) # Write message in stderr
			}
		}
		chunk_fast5_info$rname <- NULL
		chunk_fast5_info$pos <- NULL
		chunk_fast5_info$qwidth <- NULL
		chunk_fast5_info$cigar <- NULL
		chunk_fast5_info <- droplevels(chunk_fast5_info)

		if(nrow(chunk_fast5_info)>0){ # If remaining mapped reads
			chunk_fast5_info$chunk_id <- chunk
			chunk_fast5_info$contig_name <- as.factor(chunk_contig_name)
		}else{
			chunk_fast5_info <- data.frame(qname=c(NA,NA), strand=c("+","-"), chunk_id=chunk, contig_name=chunk_contig_name) # Dummy return; Keep strand for mapvalues 
		}

		return(chunk_fast5_info)
	}

	# Remove unmapped read, not necessary?
	chunks_fast5_info <- subset(chunks_fast5_info, strand %in% c("+","-"))
	chunks_fast5_info <- droplevels(chunks_fast5_info)
	idx_strand_avail <- match(levels(chunks_fast5_info$strand), c("+","-"))
	chunks_fast5_info$strand <- mapvalues(chunks_fast5_info$strand, from=c("+","-")[idx_strand_avail], to=c("fwd","rev")[idx_strand_avail])

	# Create index from read name to fast5 file
	print_message(paste0("  Link fast5 to fasta for ",sample_name,""))
	fasta <- readDNAStringSet(paste0(path_input,sample_name,".fasta"))
	index_fast5ToFasta <- do.call(rbind,strsplit(names(fasta)," "))[,c(1,2,3)]
	colnames(index_fast5ToFasta) <- c("read_name","read_id","fast5_path")
	index_fast5ToFasta <- as.data.frame(index_fast5ToFasta, stringsAsFactors=FALSE)

	# Correct path fast5 files in fasta with relative path between data and path_output if not full path
	if(!grepl("^/",index_fast5ToFasta$fast5_path[1])){
		path_correction <- relative.path(path_input, path_output)
		index_fast5ToFasta$fast5_path <- paste0(path_correction, index_fast5ToFasta$fast5_path)
	}
	names(fasta) <- paste(index_fast5ToFasta$read_name, index_fast5ToFasta$read_id, index_fast5ToFasta$fast5_path)

	# Filter out reads outside list_chunks
	list_read_name <- unique(chunks_fast5_info$qname)
	list_fast5_files <- subset(index_fast5ToFasta, select=-c(read_id), read_name %in% list_read_name) # Drop read_id
	fasta_names <- paste(list_fast5_files$read_name, gsub(".fast5","",basename(as.character(list_fast5_files$fast5_path))), list_fast5_files$fast5_path)
	fasta <- fasta[names(fasta) %in% fasta_names]

	# Create DNAStringSet subset; avoid near-complete fasta cached
	tmp_fasta_filename <- paste0(path_output,"tmp.",paste(range(list_chunks), collapse="_"),".",sample_name,".fasta")
	writeXStringSet(fasta, file=tmp_fasta_filename, width=80)
	fasta <- readDNAStringSet(tmp_fasta_filename)

	return(list(chunk=chunks_fast5_info, index=list_fast5_files, fasta=fasta))
}

handle.empty.chunk <- function(idx_chunk, index_wga, index_nat, sample_name_wga, sample_name_nat, genome, chunk_size, path_output, inSilico, bc_version=NA){
	wga_chunk_with_data <- unique(index_wga$chunk$chunk_id)
	nat_chunk_with_data <- unique(index_nat$chunk$chunk_id)

	if(!idx_chunk %in% wga_chunk_with_data | !idx_chunk %in% nat_chunk_with_data){
		if(!idx_chunk %in% wga_chunk_with_data & idx_chunk %in% nat_chunk_with_data){
			print_message(paste0("No mapped reads for ",sample_name_wga," chunk #",idx_chunk))
		}else if(idx_chunk %in% wga_chunk_with_data & !idx_chunk %in% nat_chunk_with_data){
			print_message(paste0("No mapped reads for ",sample_name_nat," chunk #",idx_chunk))
		}else{
			print_message(paste0("No mapped reads for ",sample_name_wga," & ",sample_name_nat," chunk #",idx_chunk))
		}

		# Create empty stat file to keep track of processed chunks
		trash <- create.empty.stat_data.file(genome, chunk_size, path_output, idx_chunk, inSilico, bc_version)

		return(FALSE)
	}else{
		return(TRUE)
	}
}

prepare.input.data <- function(index_sample, idx_chunk, path_output, sample_name, processed_reads, corr_type, nb_threads){
	print_message(paste0("  Preparing ",sample_name," input data for chunk #",idx_chunk))
	reads_chunk <- subset(index_sample[["chunk"]], chunk_id==idx_chunk)
	read_name_toProcess <- as.character(subset(reads_chunk, ! qname %in% processed_reads)$qname)
	read_name_toKeep <- as.character(subset(reads_chunk, qname %in% processed_reads)$qname)

	stifle <- foreach(dir=levels(reads_chunk$strand)) %do% {
		list_fast5_files <- subset(index_sample[["index"]], read_name %in% subset(reads_chunk, ! qname %in% processed_reads & strand==dir)$qname)

		if(corr_type=="nanopolish"){
			fasta_names <- paste(list_fast5_files$read_name, gsub(".fast5","",basename(as.character(list_fast5_files$fast5_path))), list_fast5_files$fast5_path)

			fasta_subset <- index_sample[["fasta"]][names(index_sample[["fasta"]]) %in% fasta_names]
			path_fasta <- paste0(path_output,"tmp.",idx_chunk,".",sample_name,".",dir,".fasta")
			writeXStringSet(fasta_subset, file=path_fasta, width=80)
		}else if(corr_type=="nanoraw"){
			path_f5_dir <- paste0(path_output,"tmp.",idx_chunk,".",sample_name,".",dir,".fast5") # TODO maybe add fwd/rev
			system(paste("mkdir",path_f5_dir))
			registerDoMC(nb_threads)
			stifle <- foreach(f5_file=list_fast5_files$fast5_path) %dopar% { # NEW
				f5_file_path <- normalizePath(f5_file)

				system(paste0("ln -s ",f5_file_path," ",path_f5_dir,"/"))

				return(NA)
			}
			registerDoSEQ()
		}

		return(NA)
	}

	return(list(toKeep=read_name_toKeep, toProcess=read_name_toProcess))
}

correct.data <- function(corr_type, path_output, idx_chunk, sample_name_wga, sample_name_nat, genome, chunk_size, path_script, nb_threads, sig_norm, map_type){
	print_message("  Correcting mapping")

	if(corr_type=="nanopolish"){
		corrected_wga <- correct.event.data(path_output, idx_chunk, sample_name_wga, genome, chunk_size, path_script, nb_threads, sig_norm, map_type) #np normalized
		corrected_nat <- correct.event.data(path_output, idx_chunk, sample_name_nat, genome, chunk_size, path_script, nb_threads, sig_norm, map_type) #np normalized
	}else if(corr_type=="nanoraw"){ # Not used
		# corrected_wga <- correct.raw.data(path_output, idx_chunk, sample_name_wga, genome, path_script, nb_threads) #median normalized
		# corrected_nat <- correct.raw.data(path_output, idx_chunk, sample_name_nat, genome, path_script, nb_threads) #median normalized
	}
	corrected_data <- rbind(corrected_wga,corrected_nat)

	if(length(unique(corrected_data$data_type))==2){ # If data from both type; TODO modify for in silico
		return(corrected_data)
	}else{
		return(NULL)
	}
}

correct.event.data <- function(path_output, idx_chunk, sample_name, genome, chunk_size, path_script, nb_threads, sig_norm, map_type){
	contigs_chunks <- generate.genome.chunks.information(genome, chunk_size)
	chunk_contig_name <- as.character(subset(contigs_chunks, chunk_id==idx_chunk)$contig_name)

	path_fasta <- paste0(path_output,"tmp.",idx_chunk,".",sample_name,".fasta")

	# Realign events with error "handleling"
	res <- NULL
	attempt <- 0
	while(is.null(res)){
		attempt <- attempt + 1
		if(attempt > 2){
			print_message(paste0("  Chunk #",idx_chunk," not processed"))
			clean.temporary.files(path_output, idx_chunk)

			stop("Correction failed.")
		}else{
			gc()
			res <- tryCatch({
				system(paste0(path_script,"realign_events.sh ",genome," ",path_fasta," ",nb_threads," ",idx_chunk," ",sig_norm," ",map_type), intern=TRUE, ignore.stderr=TRUE)
			}, warning = function(w) {
				print(w)
				print_message(paste0("  Failed correction: ",idx_chunk))
				error.handling(as.integer(gsub("status ","",str_match(w$message, "status [0-9]+")[1,1])))

				return(NULL)
			}, error = function(e) {
				print(e)
				print_message(paste0("  Failed correction: ",idx_chunk))

				return(NULL)
			})
		}
	}

	# Load corrected events
	path_fwd <- gsub(".fasta",".fwd.eventalign",path_fasta)
	realign_fwd <- fread(path_fwd, header=T, sep="\t", stringsAsFactors=TRUE, showProgress=FALSE) #showProgress=FALSE
	realign_fwd$contig <- as.factor(realign_fwd$contig)

	path_rev <- gsub(".fasta",".rev.eventalign",path_fasta)
	realign_rev <- fread(path_rev, header=T, sep="\t", stringsAsFactors=TRUE, showProgress=FALSE) #showProgress=FALSE
	realign_rev$contig <- as.factor(realign_rev$contig)
	
	# Discard unconcordant mapped contigs; Can happenned with multi-mapped reads
	realign_fwd <- subset(realign_fwd, contig==chunk_contig_name)
	realign_rev <- subset(realign_rev, contig==chunk_contig_name)

	# TODO older version of Nanopolish can report same event multiple time. Considere filtering for uniqueness.

	if(nrow(realign_fwd)>0){ # Do not proceed if no data
		realign_fwd$dir <- as.factor(rep("fwd",nrow(realign_fwd)))
	}
	if(nrow(realign_rev)>0){ # Do not proceed if no data
		if(sig_norm=="revc"){
			# Convert position to follow original order
			list_contig_shortnames <- sapply(strsplit(names(readDNAStringSet(genome))," "), `[`, 1) # Cannot match names end because only first word is ID
			current_contig_shortnames <- gsub("\\|","\\\\|",paste0("^",as.character(unique(realign_rev$contig)),"$")) # Escape | TODO others still problematic
			length_current_contig <- width(readDNAStringSet(genome))[grepl(current_contig_shortnames,list_contig_shortnames)]

			realign_rev$position <- ((length_current_contig - realign_rev$position) + 1) - 7
		}else if(sig_norm=="ori"){
			# Convert reference_kmer to complement because np only based on forward sequence
			realign_rev$reference_kmer <- as.character(reverseComplement(DNAStringSet(realign_rev$reference_kmer)))
			realign_rev$reference_kmer <- as.factor(realign_rev$reference_kmer)
		}
		realign_rev$dir <- as.factor(rep("rev",nrow(realign_rev)))
	}

	# Merge directionnal dataset and handle missing data 
	if(nrow(realign_fwd)>0 & nrow(realign_rev)>0){
		realign <- rbind(realign_fwd, realign_rev)
		rm(realign_fwd, realign_rev)
		gc()
	}else if(nrow(realign_fwd)>0 & nrow(realign_rev)==0){
		realign <- realign_fwd
		rm(realign_fwd, realign_rev)
		gc()
	}else if(nrow(realign_fwd)==0 & nrow(realign_rev)>0){
		realign <- realign_rev
		rm(realign_fwd, realign_rev)
		gc()
	}else if(nrow(realign_fwd)==0 & nrow(realign_rev)==0){

		return(NULL) # Handled in main loop
	}

	realign$data_type <- as.factor(rep(sample_name,nrow(realign)))
	realign$contig <- droplevels(realign$contig)

	return(realign)
}

extract.corrected.events <- function(f5_file, corrgroup_name, strand_type){
	tmp_corrected_events <- h5read(f5_file,paste0("/Analyses/",corrgroup_name,"/BaseCalled_",strand_type))$Events
	nb_events <- nrow(tmp_corrected_events)

	f5_structure <- h5ls(f5_file)
	path_uuid <- f5_structure$group[grep("/Raw/Reads/Read_", f5_structure$group)]
	tmp_uuid <- h5readAttributes(f5_file,path_uuid)
	tmp_corrected_events$read_name <- rep(tmp_uuid$read_id,nb_events)
	tmp_corrected_events$strand <- substr(strand_type,1,1)
	# TODO remove base?
	colnames(tmp_corrected_events)[1] <- "mean" #Not robust
	colnames(tmp_corrected_events)[2] <- "stdv" #Not robust

	tmp_attr <- h5readAttributes(f5_file,paste0("/Analyses/",corrgroup_name,"/BaseCalled_",strand_type,"/Alignment"))
	H5close() # TODO change place?
	tmp_corrected_events$dir <- ifelse(tmp_attr$mapped_strand=="+",ifelse(strand_type=="template","fwd","rev"),ifelse(strand_type=="template","rev","fwd"))
	offset <- ifelse(tmp_attr$mapped_strand=="+",-2,-3) # -4 -3
	event_index <- seq(0,nb_events-1)
	if(tmp_attr$mapped_strand=="-"){
		event_index <- rev(event_index)
	}
	tmp_corrected_events$corr_event_index <- event_index
	tmp_corrected_events$position <- tmp_attr$mapped_start + tmp_corrected_events$corr_event_index + offset # offset to have the same alignment than nanopolish

	return(tmp_corrected_events)
}

filter.mapped.reads <- function(corrected_data, min_read_length){
	corrected_data <- corrected_data %>%
		group_by(read_name, data_type, dir) %>%
		filter(length(unique(position))>=min_read_length)

	return(corrected_data)
}

ont.normalize <- function(corrected_data, name_mean_col, t_model, normalized){
	corrected_data <- merge(corrected_data,subset(t_model,select=c("kmer","level_mean")),by.x=c("reference_kmer"),by.y=c("kmer"))

	names(corrected_data)[names(corrected_data)==name_mean_col] <- "mean" # Set back?

	if(normalized==1){
		set.seed(101) # Maybe not necessary
		fit <- lm(level_mean ~ mean, data=corrected_data)
	}else if(normalized==2){
		set.seed(101) # Maybe not necessary
		fit <- rlm(level_mean ~ mean, data=corrected_data, maxit=20) # TODO maxit > 20, but should be enough in most cases
		# Can fail to converge but partial fitting results are still returned. Bad event level will be handled in OR step.
		# TODO catch warning and wrote summary (e.g. nb reads where normalization failed to converge) instead.
		# TODO could also remove those reads from downstream analysis.
	}
	ev_shift <- fit$coefficients[[1]]
	ev_scale <- fit$coefficients[[2]]
	corrected_data$norm_mean <- corrected_data$mean*ev_scale + ev_shift 

	return(subset(corrected_data, select=-c(level_mean)))
}

normalize.data.serial <- function(corrected_data, name_mean_col, idx_chunk, t_model, normalized){
	corrected_data <- tryCatch({
			ddply(corrected_data, .(read_name), ont.normalize, name_mean_col=name_mean_col, t_model=t_model, normalized=normalized, .parallel=FALSE)
		}, warning=function(w) {
			print(w) # TODO remove from final release.
			
			if(grepl("'rlm' failed to converge",w$message)){
				# Rare warning but should not impact results too much. Likely small or mismapped reads or too many difference in reference. 
				# TODO Could increase maxit=20 but likely not fixing issue
				# Not completely fixed with minReadLength
				
				cpt_warning <- 0 
				withCallingHandlers({  # Will silenced 'rlm' warning
						corrected_data <- ddply(corrected_data, .(read_name), ont.normalize, name_mean_col=name_mean_col, t_model=t_model, normalized=normalized, .parallel=FALSE)
					}, warning=function(w) if(grepl("'rlm' failed to converge",w$message)){
						cpt_warning <<- cpt_warning + 1
						invokeRestart("muffleWarning")
					}
				)
				write(paste0("	Normalization did not reach convergence for ",cpt_warning," read(s) on chunk #",idx_chunk,"."), stderr()) # Write message in stderr

				return(corrected_data)
			}else{
				print_message("Unexpected warning")

				return(corrected_data) # Will break next step
			}
		}, error=function(e) {
			print(e)
			print_message("Unexpected error")

			return(corrected_data) # Will break next step
		}
	)

	return(corrected_data)
}

normalize.data.parallel <- function(corrected_data, name_mean_col, idx_chunk, t_model, normalized, nb_threads){
	corrected_data <- tryCatch({
		registerDoMC(cores=nb_threads)
		ddply(corrected_data, .(read_name), ont.normalize, name_mean_col=name_mean_col, t_model=t_model, normalized=normalized, .parallel=TRUE)
	}, warning=function(w) {
		print(w)
		if(grepl("scheduled core[s]? .*encountered errors in user code",w$message)){
			registerDoSEQ()
			corrected_data <- normalize.data.serial(corrected_data, name_mean_col, idx_chunk, t_model, normalized)

			return(corrected_data)
		}else if(grepl("'rlm' failed to converge",w$message)){
			# Rare warning but should not impact results too much. Likely small or mismapped reads or too many difference in reference. 
			# TODO Could increase maxit=20 but likely not fixing issue
			# Not completely fixed with minReadLength
			
			# TODO can I avoid recomputing here?
			corrected_data <- tryCatch({ # Will silenced 'rlm' warning
				cpt_warning <- 0 
				withCallingHandlers(
						corrected_data <- ddply(corrected_data, .(read_name), ont.normalize, name_mean_col=name_mean_col, t_model=t_model, normalized=normalized, .parallel=TRUE)
					, warning=function(w) if(grepl("'rlm' failed to converge",w$message)){
						cpt_warning <<- cpt_warning + 1
						invokeRestart("muffleWarning")
					}
				)
				write(paste0("	Normalization did not reach convergence for ",cpt_warning," read(s) on chunk #",idx_chunk,"."), stderr()) # Write message in stderr
				
				return(corrected_data)
			}, warning=function(w) {
				print(w)
				if(grepl("scheduled core[s]? .*encountered errors in user code",w$message)){
					registerDoSEQ()
					corrected_data <- normalize.data.serial(corrected_data, name_mean_col, idx_chunk, t_model, normalized)

					return(corrected_data)
				}else{
					print_message("Unexpected warning")

					return(corrected_data) # Will break next step
				}
			}, error=function(e) {
				print_message("Unexpected error")

				return(corrected_data) # Will break next step
			})

			return(corrected_data)
		}else{
			# print_message("Too large for preschedule")
			opts <- list(.options.multicore=list(preschedule=FALSE))
			corrected_data <- ddply(corrected_data, .(read_name), ont.normalize, name_mean_col=name_mean_col, t_model=t_model, normalized=normalized, .parallel=TRUE, .paropts=opts)
		
			return(corrected_data)
		}
	}, error=function(e) {
		registerDoSEQ()
		print(e)
		if(grepl("Cannot allocate memory",e$message)){
			corrected_data <- normalize.data.serial(corrected_data, name_mean_col, idx_chunk, t_model, normalized)

			return(corrected_data)
		}else{
			print_message("Unexpected error")

			return(corrected_data) # Will break next step
		}
	})
	registerDoSEQ()

	return(corrected_data)
}

normalize.data <- function(corrected_data, name_mean_col, idx_chunk, nb_threads, t_model, normalized){
	print_message("  Normalization")

	gc()
	if(nb_threads>1){
		corrected_data <- normalize.data.parallel(corrected_data, name_mean_col, idx_chunk, t_model, normalized, nb_threads)
	}else{
		corrected_data <- normalize.data.serial(corrected_data, name_mean_col, idx_chunk, t_model, normalized)
	}

	return(corrected_data)
}

remove.outlier <- function(corrected_data, idx_chunk, chunk_size, genome, IQR_factor){
	print_message("  Removing potential outliers")
	contigs_chunks <- generate.genome.chunks.information(genome, chunk_size)
	chunk_contig_name <- as.character(subset(contigs_chunks, chunk_id==idx_chunk)$contig_name)
	chunk_start <- subset(contigs_chunks, chunk_id==idx_chunk)$start
	chunk_end <- subset(contigs_chunks, chunk_id==idx_chunk)$end

	# TODO capping? stat test are robust to them?
	trimmed_data <- subset(corrected_data, contig==chunk_contig_name & position>chunk_start & position<=chunk_end) %>%
		group_by(position, dir, data_type) %>%
		mutate(iqr=IQR(norm_mean), q25=quantile(norm_mean)[[2]], q75=quantile(norm_mean)[[4]]) %>%
		filter(norm_mean >= q25 - (IQR_factor * iqr) & norm_mean <= q75 + (IQR_factor * iqr)) %>%
		dplyr::select(-c(iqr, q25, q75))
	trimmed_data <- as.data.frame(trimmed_data)

	corrected_data <- rbind(subset(corrected_data, contig==chunk_contig_name & (position<=chunk_start | position>chunk_end)), trimmed_data)

	return(corrected_data)
}

create.empty.stat_data.file <- function(genome, chunk_size, path_output, idx_chunk, inSilico, bc_version=NA){
	contigs_chunks <- generate.genome.chunks.information(genome, chunk_size)
	chunk_contig_name <- as.character(subset(contigs_chunks, chunk_id==idx_chunk)$contig_name)
	chunk_start <- subset(contigs_chunks, chunk_id==idx_chunk)$start

	if(inSilico=="no"){
		empty_stat_data <- data.frame(
			contig=chunk_contig_name, position=chunk_start, dir="rev", strand="t", # Keep dummy value to have consitant factor later
			N_wga=0, N_nat=0, mean_diff=NA, t_test_pval=NA, u_test_pval=NA
		)
	}else{
		empty_stat_data <- data.frame(
			contig=chunk_contig_name, position=chunk_start, data_type="No_data", dir="rev", strand="t", # Keep dummy value to have consitant factor later
			N_val=0, mean_diff=NA, t_test_pval=NA, u_test_pval=NA
		)
	}

	# Add information about base caller version
	if(!any(is.na(bc_version))){
		attr(empty_stat_data, "basecaller") <- bc_version[1]
		attr(empty_stat_data, "version") <- bc_version[2]
	}

	saveRDS(empty_stat_data,file=paste0(path_output,"chunk.",idx_chunk,".difference.rds"))

	return(empty_stat_data)
}

scoring.position <- function(sub_corrected_data, inSilico, inSilico_model, min_coverage){
	# Record min output
	contig_tested <- unique(sub_corrected_data$contig)
	pos_tested <- unique(sub_corrected_data$position)
	dir_tested <- unique(sub_corrected_data$dir)

	sub_res <- sub_corrected_data %>%
		group_by(contig, position, data_type, dir, strand, reference_kmer) %>%
		summarize(means=list(norm_mean), .groups="drop_last") %>%
		group_by(contig, position, data_type, dir, strand) %>%
		filter(length(unlist(means))>=min_coverage) %>% # default min 5 value by groups
		group_by(contig, position, dir, strand)

	if((inSilico=="no" & nrow(sub_res)<2) | (inSilico!="no" & nrow(sub_res)==0)){ # When there is not enough values per dataset
		if(inSilico=="no"){
			res <- data.frame(
				contig=contig_tested, position=pos_tested, 
				dir=dir_tested, strand="t", N_wga=0, N_nat=0,
				mean_diff=NA, t_test_pval=NA, u_test_pval=NA
			)
		}else{
			res <- data.frame(
				contig=contig_tested, position=pos_tested, data_type="Low_data",
				dir=dir_tested, strand="t", N_val=0, # random dir
				mean_diff=NA, t_test_pval=NA, u_test_pval=NA
			)
		}
	}else if(inSilico=="no"){ # TODO
		res <- sub_res %>%
			group_by(contig, position, dir, strand) %>%
			filter(length(data_type)>1) %>%
			spread(data_type, means) %>%
			mutate(N_wga=length(unlist(wga)), N_nat=length(unlist(nat))) %>%
			mutate(mean_diff=mean(unlist(nat))-mean(unlist(wga))) %>%
			mutate(t_test_pval=ifelse(N_wga>=5 & N_nat>=5, tryCatch(t.test(unlist(wga), unlist(nat), var.equal=TRUE)$p.value, error=function(e) NA), NA)) %>% # try needed because ifelse evaluate both results
			mutate(u_test_pval=ifelse(N_wga>=5 & N_nat>=5, try(wilcox.test(unlist(wga), unlist(nat), correct=FALSE)$p.value, silent=TRUE), NA)) %>%				# Old Nanopolish can report same event multiple time (bug) creating error with t.test
			# mutate(LLR=compute.LLR(wga,nat)) %>% # Error "L-BFGS-B needs finite values of 'fn'"
			dplyr::select(-c(wga, nat, reference_kmer))
	}else if(inSilico=="ONT"){
		sub_res <- merge(sub_res, inSilico_model, by.x="reference_kmer", by.y="kmer") # TODO include before but modify group_by?

		res <- sub_res %>%
			group_by(contig, position, data_type, dir, strand, level_mean) %>%
			mutate(N_val=length(unlist(means))) %>%
			mutate(mean_diff=mean(unlist(means))-level_mean) %>%
			mutate(t_test_pval=ifelse(N_val>=5, tryCatch(t.test(unlist(means), mu=level_mean, var.equal=TRUE)$p.value, error=function(e) NA), NA)) %>%
			mutate(u_test_pval=ifelse(N_val>=5, try(wilcox.test(unlist(means), mu=level_mean, correct=FALSE)$p.value, silent=TRUE), NA)) %>%
			ungroup() %>%
			dplyr::select(c(contig, position, data_type, dir, strand, N_val, mean_diff, t_test_pval, u_test_pval))
	}

	return(res)
}

generate.kmer <- function(s_seq, kmer_len){
	seq_length <- nchar(s_seq)

	return(substring(s_seq, 1:(seq_length-kmer_len+1), kmer_len:seq_length))
}

compute.statistic <- function(corrected_data, idx_chunk, chunk_size, sample_name_wga, sample_name_nat, path_output, reads_toSave, nb_threads, genome, inSilico, inSilico_model, min_coverage, bc_version=NA){
	print_message("  Computing stats by genomic position")
	contigs_chunks <- generate.genome.chunks.information(genome, chunk_size)
	chunk_contig_name <- as.character(subset(contigs_chunks, chunk_id==idx_chunk)$contig_name)
	chunk_start <- subset(contigs_chunks, chunk_id==idx_chunk)$start
	chunk_end <- subset(contigs_chunks, chunk_id==idx_chunk)$end

	if(inSilico!="no"){
		if(inSilico=="ONT"){
			# handled in scoring.position
			inSilico_model <- inSilico_model # it's t_model
		}
	}else{
		inSilico_model <- NA
	}

	chunk_corrected_data <- subset(corrected_data, contig==chunk_contig_name & position>=chunk_start & position<=chunk_end)
	if(nrow(chunk_corrected_data)>0){ # If reads mapped on chunk; multimapped reads might be lost after remapping
		chunk_corrected_data$data_type <- mapvalues(chunk_corrected_data$data_type, c(sample_name_wga,sample_name_nat), c("wga","nat")) # Reverse not needed because not returned
		chunk_stat_data <- NULL

		tryCatch({
				gc()
				registerDoMC(nb_threads) # Warning if nb_threads==1 but ok
				chunk_stat_data <- ddply(chunk_corrected_data, .(position, dir), scoring.position, inSilico=inSilico, inSilico_model=inSilico_model, min_coverage=min_coverage, .parallel=TRUE)
			}
		)
		registerDoSEQ()

		if(is.null(chunk_stat_data)){ # When parallel processing failed
			tryCatch({ # Can still exceed memory
					gc()
					chunk_stat_data <- ddply(chunk_corrected_data, .(position, dir), scoring.position, inSilico=inSilico, inSilico_model=inSilico_model, min_coverage=min_coverage, .parallel=FALSE)
				}, error = function(e) {
					print_message(paste0("  Failed scoring chunk #",idx_chunk," (e)"))
					print(e)
					clean.temporary.files(path_output, idx_chunk)
					stop(e)
				}
			)
		}
	}else{
		chunk_stat_data <- create.empty.stat_data.file(genome, chunk_size, path_output, idx_chunk, inSilico, bc_version)
	}
	
	# Add information about base caller version
	if(!any(is.na(bc_version))){
		attr(chunk_stat_data, "basecaller") <- bc_version[1]
		attr(chunk_stat_data, "version") <- bc_version[2]
	}

	saveRDS(chunk_stat_data,file=paste0(path_output,"chunk.",idx_chunk,".difference.rds"))

	return(chunk_stat_data)
}




