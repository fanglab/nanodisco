#!/usr/bin/env Rscript

suppressMessages(library(optparse))

# Parsing arguments
option_list <- list(
	make_option(c("-x", "--exec_type"), type="character", default="seq", help="Execution type (seq or batch)", metavar="<batch|seq>"),
	make_option(c("-f", "--first_chunk"), type="integer", default=NULL, help="Chunk index", metavar="<chunk_id>"),
	make_option(c("-l", "--last_chunk"), type="integer", default=NULL, help="Chunk index", metavar="<chunk_id>"),
	make_option(c("-p", "--nb_threads"), type="integer", default=1, help="Number of threads to use for each chunk", metavar="<integer>"),
	make_option(c("-w", "--sample_name_wga"), type="character", default=NULL, help="WGA sample name ", metavar="<name>"),
	make_option(c("-n", "--sample_name_nat"), type="character", default=NULL, help="Native sample name", metavar="<name>"),
	make_option(c("-i", "--path_input"), type="character", default=NULL, help="Path input", metavar="<path>"),
	make_option(c("-o", "--path_output"), type="character", default=NULL, help="Path output", metavar="<path>"),
	make_option(c("-r", "--genome"), type="character", default=NULL, help="Path to reference genome (.fasta)", metavar="<path>"),
	make_option(c("-a", "--IQR_factor"), type="double", default=0, help="IQR factor for outliers removal (0 to skip; smaller is harsher)", metavar="<double>"),
	make_option(c("-z", "--normalized"), type="integer", default=1, help="Either normalize or not (1|0|2)", metavar="<0|1|2>"),
	make_option(c("-s", "--chunk_size"), type="integer", default=5000, help="Size of chunk to use (default is 5000)", metavar="<integer>"),
	make_option(c("-c", "--corr_type"), type="character", default="nanopolish", help="Type of event correction (nanopolish|nanoraw)", metavar="<nanopolish>"),
	make_option(c("-b", "--sig_norm"), type="character", default="revc", help="Correct strand bias (ori|revc)", metavar="<ori|revc>"),
	make_option(c("-q", "--inSilico"), type="character", default="no", help="Analyze using in silico model (no|ONT)", metavar="<no|ONT|other>"),
	make_option(c("-d", "--inSilico_model"), type="character", default=NULL, help="Path to alternate model", metavar="<path>"),
	make_option(c("-e", "--min_coverage"), type="integer", default=5, help="Minimum number of events per position", metavar="<integer>"),
	make_option(c("-j", "--map_type"), type="character", default="noAddSupp", help="Type of filtering for mapping", metavar="<type>"),
	make_option(c("-k", "--min_read_length"), type="integer", default=0, help="Minimum mapped read length (0 default)", metavar="<integer>"),
	make_option(c("-t", "--path_t_model"), type="character", default=NULL, help="Path for nanopolish r9.4_450bps.nucleotide.6mer.template.model", metavar="<path>"),
	make_option(c("--basecall_version"), type="character", default="default", help="Basecalling version (when multiple ones available)", metavar="<basecaller:version>")
)

default_usage <- c("nanodisco difference -nj <nb_jobs> -nc <nb_chunks> -p <nb_threads> -i <path_input> -o <path_output> -w <name_WGA> -n <name_native> -r <path_genome> [-f <first_chunk> -l <last_chunk> + advanced parameters]")
opt_parser <- OptionParser(usage=default_usage, option_list=option_list)
opt <- parse_args(opt_parser)

# Load parameters and functions
path_script <- "/home/nanodisco/code/"
source(paste0(path_script,"difference_functions.R"), keep.source=TRUE)
suppressMessages(load.libraries())

setDTthreads(threads=opt$nb_threads)

if(opt$exec_type=="batch"){
	if(is.null(opt$first_chunk)|is.null(opt$last_chunk)){
		print_help(opt_parser)
		stop("Chunk index needed.", call.=FALSE)
	}else{
		list_chunks <- seq(opt$first_chunk,opt$last_chunk)
	}
	if(is.null(opt$nb_threads)){
		print_help(opt_parser)
		stop("Number of CPU needed.", call.=FALSE)
	}else{
		nb_threads <- opt$nb_threads
	}

	path_t_model <- opt$path_t_model
	t_model <- fread(path_t_model, sep="\t", header=TRUE, skip=5)

	# TODO add control
	sample_name_wga <- opt$sample_name_wga
	sample_name_nat <- opt$sample_name_nat
	path_input <- paste0(normalizePath(opt$path_input),"/") # Deal with relative path, stay usable afetr setwd()
	if(!dir.exists(opt$path_output)){
		dir.create(opt$path_output)
	}
	path_output <- paste0(normalizePath(opt$path_output),"/") # Deal with relative path, stay usable afetr setwd()
	genome <- normalizePath(opt$genome) # Deal with relative path, stay usable afetr setwd()

	IQR_factor <- opt$IQR_factor # 0 | +/- 1.5
	normalized <- opt$normalized # 0 | 1 | 2
	chunk_size <- opt$chunk_size # 5000
	corr_type <- opt$corr_type # nanopolish | nanoraw
	sig_norm <- opt$sig_norm # ori | revc
	inSilico <- opt$inSilico # no | ONT
	if(inSilico=="ONT"){ # TODO check inSilico also
		inSilico_model <- t_model
	}else{
		inSilico_model <- opt$inSilico_model # Is null
	}
	min_coverage <- opt$min_coverage # default 5
	map_type <- opt$map_type # default 5
	min_read_length <- opt$min_read_length # default 5

	basecall_version <- opt$basecall_version # default is default
}else if(opt$exec_type=="seq"){ # If launched using source; use default values
	
	# Not used
}else{
	# If input error
	stop("Unknown execution types.", call.=FALSE)
}

setwd(path_input)

check.reference(genome)

bc_version <- check.basecall.version(path_input, sample_name_nat, sample_name_wga, basecall_version)

Sys.time()
index_wga <- prepare.index(sample_name_wga, path_input, path_output, genome, chunk_size, list_chunks) # TODO try to save in hdf5 files
index_nat <- prepare.index(sample_name_nat, path_input, path_output, genome, chunk_size, list_chunks) # TODO try to save in hdf5 files

Sys.time()
processed_reads_wga <- ""
processed_reads_nat <- ""
corrected_data <- ""
prev_corrected_data <- ""
final_stat_data <- foreach(idx_chunk=list_chunks, .combine=rbind) %:% when(handle.empty.chunk(idx_chunk, index_wga, index_nat, sample_name_wga, sample_name_nat, genome, chunk_size, path_output, inSilico, bc_version)) %do% {
	print(paste0("Processing chunk #",idx_chunk))
	# Create fasta subset
	reads_wga <- prepare.input.data(index_wga, idx_chunk, path_output, sample_name_wga, processed_reads_wga, corr_type, nb_threads)
	reads_nat <- prepare.input.data(index_nat, idx_chunk, path_output, sample_name_nat, processed_reads_nat, corr_type, nb_threads)

	# Removed analyzed reads
	if(!is.null(nrow(corrected_data))){ # If first chunk, could switch to a flag
		prev_corrected_data <- subset(corrected_data, read_name %in% c(reads_wga$toKeep, reads_nat$toKeep))
	}
	rm(corrected_data)

	# Correct event calling
	corrected_data <- correct.data(corr_type, path_output, idx_chunk, sample_name_wga, sample_name_nat, genome, chunk_size, path_script, nb_threads, sig_norm, map_type)
	gc()

	if(!is.null(corrected_data) | !is.null(nrow(prev_corrected_data))){ # If still data for chunk idx_chunk
		# Read filtering
		if(min_read_length>0 & !is.null(corrected_data)){
			corrected_data <- filter.mapped.reads(corrected_data, min_read_length)
		}

		# Normalization
		if(normalized>0 & !is.null(corrected_data)){
			corrected_data <- normalize.data(corrected_data, "event_level_mean", idx_chunk, nb_threads, t_model, normalized) # Reduce memory footprint by diminishing nb_threads; nb_threads-3
		}else{ # Add norm_mean column with nanopolish means without normalization
			corrected_data$norm_mean <- corrected_data$event_level_mean 
		}

		# Combine newly processed reads with previous those of chunk
		if(!is.null(nrow(prev_corrected_data))){
			corrected_data <- rbind(prev_corrected_data,corrected_data)
			rm(prev_corrected_data)
			gc()
		}

		# Outlier removal
		if(IQR_factor>0){
			corrected_data <- remove.outlier(corrected_data, idx_chunk, chunk_size, genome, IQR_factor)
			gc()
		}

		# Compute chunk statistics
		stat_data <- compute.statistic(corrected_data, idx_chunk, chunk_size, sample_name_wga, sample_name_nat, path_output, c(reads_wga$toProcess,reads_nat$toProcess), nb_threads, genome, inSilico, inSilico_model, min_coverage, bc_version)
		gc()

		# Record processed reads
		processed_reads_wga <- c(reads_wga$toKeep, reads_wga$toProcess)
		processed_reads_nat <- c(reads_nat$toKeep, reads_nat$toProcess)

		clean.temporary.files(path_output, idx_chunk)

		return(stat_data)
	}else{ # If no data for chunk idx_chunk; e.g. Multimapped reads/repeat region/genomic SV
		print(paste0("  No data for chunk #",idx_chunk))
		
		# Create empty stat file to keep track of processed chunks
		stat_data <- create.empty.stat_data.file(genome, chunk_size, path_output, idx_chunk, inSilico, bc_version)

		clean.temporary.files(path_output, idx_chunk)

		return(stat_data)
	}
}
stifle <- file.remove(paste0(path_output,"tmp.",paste(range(list_chunks), collapse="_"),".",sample_name_wga,".fasta")) # . ok because not pattern
stifle <- file.remove(paste0(path_output,"tmp.",paste(range(list_chunks), collapse="_"),".",sample_name_nat,".fasta"))
Sys.time()

