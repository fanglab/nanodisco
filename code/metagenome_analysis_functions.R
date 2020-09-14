
source("/home/nanodisco/code/difference_functions.R")
source("/home/nanodisco/code/analysis_functions.R")

path_list_motifs <- "/home/nanodisco/reference/list_motifs.RDS"

load.libraries.metagenome <- function(){
	library(optparse)
	library(stringr)
	library(Biostrings) # For sequence composition
	library(Rsamtools)
	library(seqinr)
	library(egg)
	library(scales)
	library(plyr)
	library(dplyr)
	library(doMC)
	library(hues)
	library(tidyr)
	library(Rtsne)
	library(progress)
}

print_message <- function(message){
	# Print message to terminal
	cat(paste0("[",Sys.time(),"] ",message,".\n"))
}

#   _____ _               _      _                   _       
#  /  __ \ |             | |    (_)                 | |      
#  | /  \/ |__   ___  ___| | __  _ _ __  _ __  _   _| |_ ___ 
#  | |   | '_ \ / _ \/ __| |/ / | | '_ \| '_ \| | | | __/ __|
#  | \__/\ | | |  __/ (__|   <  | | | | | |_) | |_| | |_\__ \
#   \____/_| |_|\___|\___|_|\_\ |_|_| |_| .__/ \__,_|\__|___/
#                                       | |                  
#                                       |_|                  

check.input.profile <- function(opt){
	# Check if current differences file exist.
	if(is.null(opt$path_difference)){
		cat("Parameter -d/--path_difference is missing. Please provide the path to a current differences file.\n")
		# print_help(opt_parser)
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$path_difference)){
			cat(paste0("Current differences file doesn't exist (",opt$path_difference,").\n"))
			cat("Please check -d/--path_difference parameter and/or run nanodisco difference.\n")
			quit(save="no", status=4)
		}
	}
	# Check if output directory exist and create if not.
	if(!dir.exists(opt$path_output)){
		dir.create(opt$path_output, recursive=TRUE)
	}
	# Check if coverage files exist.
	if(is.null(opt$path_cov_wga)){
		cat("Parameter -w/--path_cov_wga is missing. Please provide the path to WGA sample coverage file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$path_cov_wga)){
			cat(paste0("WGA sample coverage file doesn't exist (",opt$path_cov_wga,").\n"))
			cat("Please check -w/--path_cov_wga parameter and/or run nanodisco coverage.\n")
			quit(save="no", status=4)
		}
	}
	if(is.null(opt$path_cov_nat)){
		cat("Parameter -n/--path_cov_nat is missing. Please provide the path to WGA sample coverage file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$path_cov_nat)){
			cat(paste0("Native sample coverage file doesn't exist (",opt$path_cov_nat,").\n"))
			cat("Please check -n/--path_cov_nat parameter and/or run nanodisco coverage.\n")
			quit(save="no", status=4)
		}
	}
	# Check if reference metagenome file exist.
	if(is.null(opt$metagenome)){
		cat("Parameter -r/--metagenome is missing. Please provide the path to a reference metagenome file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$metagenome)){
			cat(paste0("Reference metagenome file doesn't exist (",opt$metagenome,").\n"))
			cat("Please check -r/--metagenome parameter. Path to reference metagenome (.fasta or .fa).\n")
			quit(save="no", status=4)
		}
	}
	# Check if methylation profiling processing type is supplied (automated with -a/--auto or motif discovery with -m/--motifs_file)
	if(is.null(opt$list_motif) && is.null(opt$motifs_file)){
		if(!is.null(opt$auto)){
			if(opt$auto=="all"){
				if(!file.exists(path_list_motifs)){
					cat(paste0("List of common methylation motifs (",path_list_motifs,") is missing. Please reach out for help on GitHub.\n"))
					quit(save="no", status=4)
				}

				cat(paste0("Methylation profile are computed for all predefined common motifs (n=210,176) on long contigs only (>=",opt$min_contig_len," bp). This can take a while (>24h).\n"))
				
				print_message("Read list of common motifs")
				opt$list_motif <- paste0(readRDS(path_list_motifs), collapse=",")
			}else if(opt$auto=="4mer"){
				print_message("Generate list of common 4-mer motifs")
				opt$list_motif <- paste0(generate.motifs.dictionnary(4, 4, FALSE), collapse=",")
			}else if(opt$auto=="5mer"){
				print_message("Generate list of common 5-mer motifs")
				opt$list_motif <- paste0(generate.motifs.dictionnary(5, 5, FALSE), collapse=",")
			}else if(opt$auto=="6mer"){
				print_message("Generate list of common 6-mer motifs")
				opt$list_motif <- paste0(generate.motifs.dictionnary(6, 6, FALSE), collapse=",")
			}else if(opt$auto=="noBi"){
				print_message("Generate list of common non-bipartite motifs")
				opt$list_motif <- paste0(generate.motifs.dictionnary(4, 6, FALSE), collapse=",")
			}else{
				cat("Parameter value for -a/--auto isn't recognized. Please select the appropriate option (all|4mer|5mer|6mer|noBi).\n")
				quit(save="no", status=4)
			}
		}else{
			cat("Parameter -a/--auto or -m/--list_motif/--motifs_file is missing. Please choose between:\n")
			cat("\t- adding -a/--auto to compute profile from predefined common motifs\n\t- providing a comma separated list of motifs (e.g. GATC,CCWGG) with -m/--list_motif\n\t- providing the path to a file with list of motifs (one per line) with --motifs_file.\n")
			quit(save="no", status=3)
		}
	}else{
		if(!is.null(opt$auto)){
			cat("Parameter -a/--auto and -m/--list_motif/--motifs_file cannot be both supplied. Please select the appropriate option.\n")
			quit(save="no", status=4)
		}else{
			if(!is.null(opt$list_motif) && !is.null(opt$motifs_file)){
				cat("Parameter -m/--list_motif and --motifs_file cannot be both supplied. Please select one of the two options.\n")
				quit(save="no", status=4)
			}else{
				if(!is.null(opt$list_motif)){
					# Check list of motif formating
					if(!grepl(pattern="^[ACGTRYSWKMBDHVN,]+$", x=opt$list_motif)){
						cat(paste0("Unknown character found in comma separated list of motifs (",opt$list_motif,").\n"))
						cat("Please check -m/--list_motif parameter. Only the following characters are recognized: 'ACGTRYSWKMBDHVN,' (e.g. GATC,CCWGG).\n")
						quit(save="no", status=4)
					}
				}else if(!is.null(opt$motifs_file)){
					if(!file.exists(opt$motifs_file)){
						cat(paste0("List of motifs file doesn't exist (",opt$motifs_file,").\n"))
						cat("Please check --motifs_file parameter. Path to file with list of motifs (one per line).\n")
						quit(save="no", status=4)
					}else{
						motifs_from_file <- read.table(opt$motifs_file, stringsAsFactors=FALSE)
						if(!grepl(pattern="^[ACGTRYSWKMBDHVN,]+$", x=paste0(motifs_from_file$V1, collapse=","))){
							cat(paste0("Unknown character found in list of motifs from ",opt$motifs_file," file.\n"))
							cat("Please check --motifs_file parameter. Only the following characters are recognized: 'ACGTRYSWKMBDHVN,' with one motif per line.\n")
							quit(save="no", status=4)
						}else{
							opt$list_motif <- paste0(motifs_from_file$V1, collapse=",")
						}
					}
				}
			}
		}
	}

	return(opt)
}

check.input.feature <- function(opt){
	# Check if methylation profile file exist.
	if(is.null(opt$path_profile)){
		cat("Parameter -s/--path_profile is missing. Please provide the path to a methylation profile file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$path_profile)){
			cat(paste0("Methylation profile file doesn't exist (",opt$path_profile,").\n"))
			cat("Please check -s/--path_profile parameter and/or run nanodisco profile.\n")
			quit(save="no", status=4)
		}else{ # Check if methylation file generated with nanodisco profile -a
			# Load supplied methylation profile matrix
			methylation_profile <- readRDS(opt$path_profile)
			original_min_contig_len <- attr(methylation_profile, "min_contig_len")
			if(is.null(original_min_contig_len)){
				cat(paste0("Methylation profile file (",opt$path_profile,") was not generated with -a parameter.\n"))
				cat("Please check -s/--path_profile parameter and/or run nanodisco profile -a [+ parameters].\n")
				quit(save="no", status=4)
			}else{
				if(!is.null(opt$fsel_min_contig_len)){
					if(opt$fsel_min_contig_len > original_min_contig_len){
						cat(paste0("Minimum length to consider a contig for feature selection was increased.\n"))
						cat(paste0("Remove --fsel_min_contig_len parameter if you want to use the same threshold than for `nanodisco profile -a`.\n"))
					}else{
						opt$fsel_min_contig_len <- original_min_contig_len
					}
				}else{
					opt$fsel_min_contig_len <- original_min_contig_len
				}
			}
		}
	}
	# Check if output directory exist and create if not.
	if(!dir.exists(opt$path_output)){
		dir.create(opt$path_output, recursive=TRUE)
	}
	# Check if reference metagenome file exist.
	if(is.null(opt$metagenome)){
		cat("Parameter -r/--metagenome is missing. Please provide the path to a reference metagenome file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$metagenome)){
			cat(paste0("Reference metagenome file doesn't exist (",opt$metagenome,").\n"))
			cat("Please check -r/--metagenome parameter. Path to reference metagenome (.fasta or .fa).\n")
			quit(save="no", status=4)
		}
	}

	return(opt)
}

check.input.filter <- function(opt){
	# Check if current differences file exist.
	if(is.null(opt$path_difference)){
		cat("Parameter -d/--path_difference is missing. Please provide the path to a current differences file.\n")
		# print_help(opt_parser)
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$path_difference)){
			cat(paste0("Current differences file doesn't exist (",opt$path_difference,").\n"))
			cat("Please check -d/--path_difference parameter and/or run nanodisco difference.\n")
			quit(save="no", status=4)
		}
	}
	# Check if selected features file exist.
	if(is.null(opt$path_feature)){
		cat("Parameter -f/--path_feature is missing. Please provide the path to a selected features file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$path_feature)){
			cat(paste0("Selected features file doesn't exist (",opt$path_feature,").\n"))
			cat("Please check -f/--path_feature parameter and/or run `nanodisco selected_feature`\n")
			quit(save="no", status=4)
		}
	}
	# Check if output directory exist and create if not.
	if(!dir.exists(opt$path_output)){
		dir.create(opt$path_output, recursive=TRUE)
	}
	# Check if reference metagenome file exist.
	if(is.null(opt$metagenome)){
		cat("Parameter -r/--metagenome is missing. Please provide the path to a reference metagenome file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$metagenome)){
			cat(paste0("Reference metagenome file doesn't exist (",opt$metagenome,").\n"))
			cat("Please check -r/--metagenome parameter. Path to reference metagenome (.fasta or .fa).\n")
			quit(save="no", status=4)
		}
	}

	return(opt)
}

check.input.binning <- function(opt){
	# Check if reference metagenome file exist.
	if(is.null(opt$metagenome)){
		cat("Parameter -r/--metagenome is missing. Please provide the path to a reference metagenome file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$metagenome)){
			cat(paste0("Reference metagenome file doesn't exist (",opt$metagenome,").\n"))
			cat("Please check -r/--metagenome parameter. Path to reference metagenome (.fasta or .fa).\n")
			quit(save="no", status=4)
		}
	}
	# Check if methylation profile file exist.
	if(is.null(opt$path_profile)){
		cat("Parameter -s/--path_profile is missing. Please provide the path to a methylation profile file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$path_profile)){
			cat(paste0("Methylation profile file doesn't exist (",opt$path_profile,").\n"))
			cat("Please check -s/--path_profile parameter and/or run nanodisco profile.\n")
			quit(save="no", status=4)
		}
	}
	# Check if output directory exist and create if not.
	if(!dir.exists(opt$path_output)){
		dir.create(opt$path_output, recursive=TRUE)
	}

	return(opt)
}

check.input.plot.binning <- function(opt){
	# Check if reference metagenome file exist.
	if(is.null(opt$metagenome)){
		cat("Parameter -r/--metagenome is missing. Please provide the path to a reference metagenome file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$metagenome)){
			cat(paste0("Reference metagenome file doesn't exist (",opt$metagenome,").\n"))
			cat("Please check -r/--metagenome parameter. Path to reference metagenome (.fasta or .fa).\n")
			quit(save="no", status=4)
		}
	}
	# Check if methylation profile file exist.
	if(is.null(opt$path_methylation_binning)){
		cat("Parameter -u/--path_methylation_binning is missing. Please provide the path to the methylation binning results.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$path_methylation_binning)){
			cat(paste0("Methylation binning file doesn't exist (",opt$path_methylation_binning,").\n"))
			cat("Please check -u/--path_methylation_binning parameter and/or run nanodisco binning.\n")
			quit(save="no", status=4)
		}
	}
	# Check if output directory exist and create if not.
	if(!dir.exists(opt$path_output)){
		dir.create(opt$path_output, recursive=TRUE)
	}
	# Check if methylation profile file exist.
	if(!is.null(opt$path_annotation)){
		if(!file.exists(opt$path_annotation)){
			cat(paste0("Contig annotation file doesn't exist (",opt$path_annotation,").\n"))
			cat("Please check -a/--path_annotation parameter.\n")
			quit(save="no", status=4)
		}else{
			# Check if RDS file
			opt$type_annotation <- tryCatch({
				tmp <- readRDS(opt$path_annotation)
				c("is_rds")
			}, warning = function(w) {
				# This is not an RDS file

				return(NULL)
			}, error = function(e) {
				# This is not an RDS file
				# Check if txt file
				opt$type_annotation <- tryCatch({
					read.table(opt$path_annotation, stringsAsFactors=FALSE)
					c("is_txt")
				}, warning = function(w) {
					# This is not an txt file

					return(NULL)
				}, error = function(e) {
					# This is not an txt file

					return(NULL)
				})

				return(type_annotation)
			})
			if(is.null(opt$type_annotation)){
				cat(paste0("Contig annotation file content not recognized (",opt$path_annotation,").\n"))
				cat("Please check -a/--path_annotation parameter. We expect two columns .txt or .RDS file with contig_name<tab>custom_name).\n")
				quit(save="no", status=4)
			}
		}
	}
	# Check if some MGEs contigs should be highlighted.
	if(is.null(opt$list_MGE_contig) && is.null(opt$MGEs_file)){
		opt$list_MGE_contig <- NA # No MGE contigs 
	}else{
		if(!is.null(opt$list_MGE_contig) && !is.null(opt$MGEs_file)){
			cat("Parameter -c/--list_MGE_contig and --MGEs_file cannot be both supplied. Please select one of the two options.\n")
			quit(save="no", status=4)
		}else{
			if(!is.null(opt$list_MGE_contig)){
				# TODO could check if contig exist in metagenome.fasta and/or in methylation profile file
			}else if(!is.null(opt$MGEs_file)){
				if(!file.exists(opt$MGEs_file)){
					cat(paste0("List of MGE contigs file doesn't exist (",opt$MGEs_file,").\n"))
					cat("Please check --MGEs_file parameter. Path to file with of MGE contigs (one per line).\n")
					quit(save="no", status=4)
				}else{
					MGEs_from_file <- read.table(opt$MGEs_file, stringsAsFactors=FALSE)
					# TODO could check if contig exist in metagenome.fasta and/or in methylation profile file
					opt$list_MGE_contig <- paste0(MGEs_from_file$V1, collapse=",")
				}
			}
		}
	}
	# Check if zooming parameters are understood
	if(!is.null(opt$xlim)){
		if(!grepl("^[+-]?[0-9]+:[+-]?[0-9]+$", x=opt$xlim)){
			cat(paste0("Optional x-axis zooming (",opt$xlim,") is not recognized.\n"))
			cat("Please check --xlim parameter (e.g. -5:10)\n")
			quit(save="no", status=4)
		}
	}else{
		opt$xlim <- NA # No zoom performed
	}
	if(!is.null(opt$ylim)){
		if(!grepl("^[+-]?[0-9]+:[+-]?[0-9]+$", x=opt$ylim)){
			cat(paste0("Optional y-axis zooming (",opt$ylim,") is not recognized.\n"))
			cat("Please check --ylim parameter (e.g. -10:9)\n")
			quit(save="no", status=4)
		}
	}else{
		opt$ylim <- NA # No zoom performed
	}

	return(opt)
}

#   _____             _   _                                     _        _   _             
#  /  __ \           | | (_)                                   | |      | | (_)            
#  | /  \/ ___  _ __ | |_ _  __ _ ___    __ _ _ __  _ __   ___ | |_ __ _| |_ _  ___  _ __  
#  | |    / _ \| '_ \| __| |/ _` / __|  / _` | '_ \| '_ \ / _ \| __/ _` | __| |/ _ \| '_ \ 
#  | \__/\ (_) | | | | |_| | (_| \__ \ | (_| | | | | | | | (_) | || (_| | |_| | (_) | | | |
#   \____/\___/|_| |_|\__|_|\__, |___/  \__,_|_| |_|_| |_|\___/ \__\__,_|\__|_|\___/|_| |_|
#                            __/ |                                                         
#                           |___/                                                          

create.contigs.info <- function(path_metagenome, genome_id, data_type){
	sequence_metagenome <- readDNAStringSet(path_metagenome)
	if(data_type=="mock"){
		contigs_info <- colsplit(gsub("(.*)_([0-9]+):([0-9]+)","\\1 \\2 \\3",names(sequence_metagenome)), " ", names=c("genome","start","end")) %>% # TODO use width
			mutate(contig=names(sequence_metagenome), length=end-start+1, id=genome_id$id[match(genome, genome_id$contig_name)])
	}else if(data_type=="mock_denovo"){
		contigs_info <- data.frame(contig=names(sequence_metagenome), length=width(sequence_metagenome)) %>%
			mutate(genome=genome_id$ref[match(contig, genome_id$contig)], id=genome_id$id[match(contig, genome_id$contig)])
	}else if(data_type=="real"){
		contigs_info <- data.frame(
			genome=str_split(names(sequence_metagenome)," ",simplify=TRUE)[,1], #names(sequence_metagenome), #NEW
			contig=names(sequence_metagenome),
			length=width(sequence_metagenome),
			start=1, #NEW
			end=width(sequence_metagenome), #NEW
			id=genome_id$id[match(names(sequence_metagenome), genome_id$contig)]
		)
	}

	return(contigs_info)
}

prepare.metagenome.info <- function(path_metagenome){
	metagenome <- readDNAStringSet(path_metagenome)
	metagenome_info <- data.frame(
		ncbi_name=names(metagenome),
		contig_name=gsub(".*(unitig_[0-9]+),.*","\\1",names(metagenome)),
		length=width(metagenome),
		chr=str_split(names(metagenome)," ",simplify=TRUE)[,1]
	)

	return(metagenome_info)
}

process.coverage <- function(input_coverage){
	data_cov <- read.table(input_coverage, header=F)
	colnames(data_cov) <- c("chr","depth","nb_pos","chrom_size","frac")
	data_cov <- subset(data_cov, chr!="genome") # Remove global coverage "genome"

	contig_annotation <- data_cov %>%
		group_by(chr) %>%
		summarize(contig_length=unique(chrom_size), avg_cov=sum(depth*nb_pos)/contig_length, .groups="drop_last")

	return(contig_annotation)
}

prepare.contig.annotation <- function(input_coverage, contig_clustering, bin_annotation){
	contig_annotation <- process.coverage(input_coverage)
	contig_annotation <- merge(contig_annotation, contig_clustering, by="chr", all.x=TRUE)
	contig_annotation <- merge(contig_annotation, bin_annotation, by="contig_name", all.x=TRUE)
	contig_annotation <- contig_annotation %>%
		mutate(cols=as.factor(ifelse(is.na(bin),ifelse(is.na(length.PacBio),"Not Clustered","Not Binned"),as.character(bin)))) %>%
		mutate(sizes=ifelse(cols=="Not Binned", 1, 3))

	return(contig_annotation)
}

prepare.contig.clustering <- function(metagenome_info, path_clustering){
	contig_clustering <- read.table(path_clustering, header=TRUE)
	contig_clustering <- merge(metagenome_info, contig_clustering, by="contig_name", suffixes=c(".ONT",".PacBio"), all.x=TRUE)

	return(contig_clustering)
}

prepare.bin.annotation <- function(path_binning){
	bin_annotation <- foreach(fasta_bin_file=list.files(path_binning, pattern="*.fasta"), .combine=rbind) %do% {
		fasta_bin <- readDNAStringSet(paste0(path_binning,"/",fasta_bin_file))
		contig_name <- str_split(names(fasta_bin),"\\|", simplify=TRUE)[,1]
		bin <- as.integer(gsub(".*([0-9]).*","\\1",fasta_bin_file, fixed=FALSE))

		return(data.frame(contig_name=contig_name, bin=paste0("Bin_",bin)))
	}

	return(bin_annotation)
}

find.information.contig <- function(contig_name, bin_id, kraken_id){
	annotation <- merge(bin_id, kraken_id, by=c("contig"), suffixes=c("bin","kraken"))
	results <- subset(annotation, grepl(contig_name, contig))

	print(results)
}

read.kraken.labels <- function(path_labels, id_column, split_all=FALSE, keep_all=FALSE){
	kraken_annotation <- read.table(path_labels, header=FALSE, stringsAsFactors=FALSE)
	colnames(kraken_annotation) <- c("contig", "Kraken_annotation")
	kraken_annotation <- kraken_annotation %>%
		mutate(contig=str_split(contig,"\\|",simplify=TRUE)[,1])
	
	if(!is.na(id_column)){
		kraken_annotation <- kraken_annotation %>%
			mutate(id=gsub("[dpcofgs]__","",str_split(Kraken_annotation,"\\|",simplify=TRUE)[,id_column]))
	}

	if(split_all){
		kraken_annotation <- kraken_annotation %>%
			mutate(Kraken_annotation=gsub("[dpcofgs]__","",Kraken_annotation)) %>%
			separate(Kraken_annotation, into=c("superkingdom","kingdom","phylum","class","order","family","genus","species"), sep="\\|")
	}

	if(!keep_all & !split_all){
		kraken_annotation <- kraken_annotation %>%
			dplyr::select(-Kraken_annotation)
	}

	return(kraken_annotation)
}

#   _____                                     
#  /  __ \                                    
#  | /  \/ _____   _____ _ __ __ _  __ _  ___ 
#  | |    / _ \ \ / / _ \ '__/ _` |/ _` |/ _ \
#  | \__/\ (_) \ V /  __/ | | (_| | (_| |  __/
#   \____/\___/ \_/ \___|_|  \__,_|\__, |\___|
#                                   __/ |     
#                                  |___/      

compute.effective.coverage <- function(methylation_signal, contigs_info){
	effective_coverage <- methylation_signal %>%
		filter(!is.na(mean_diff)) %>% # Remove position without values
		group_by(contig) %>%
		summarize(sum_N_wga=sum(N_wga, na.rm=TRUE), sum_N_nat=sum(N_nat, na.rm=TRUE), .groups="drop_last") %>%
		mutate(contig=as.character(contig)) %>%
		inner_join(contigs_info, by="contig") %>%
		mutate(avg_cov_wga=sum_N_wga/length, avg_cov_nat=sum_N_nat/length)

	return(effective_coverage)
}

#   _____                           _        ___  ___      _   _  __     
#  |  __ \                         | |       |  \/  |     | | (_)/ _|    
#  | |  \/ ___ _ __   ___ _ __ __ _| |_ ___  | .  . | ___ | |_ _| |_ ___ 
#  | | __ / _ \ '_ \ / _ \ '__/ _` | __/ _ \ | |\/| |/ _ \| __| |  _/ __|
#  | |_\ \  __/ | | |  __/ | | (_| | ||  __/ | |  | | (_) | |_| | | \__ \
#   \____/\___|_| |_|\___|_|  \__,_|\__\___| \_|  |_/\___/ \__|_|_| |___/
#                                                                        
#                                                                        

generate.ungapped.motifs.dictionnary <- function(ungapped_motif_min_length, ungapped_motif_max_length, nucleotides){
	list_ungapped_motifs <- foreach(length_motifs=seq(ungapped_motif_min_length, ungapped_motif_max_length), .combine=c) %do% {
		ungapped_motifs_base <- vector("list", length_motifs)
		for(i in seq(1,length_motifs)) {
			ungapped_motifs_base[[i]] <- nucleotides
		}

		ungapped_motifs_df <- expand.grid(ungapped_motifs_base)
		subset_ungapped_motifs <- apply(ungapped_motifs_df, 1, paste, collapse="")

		return(subset_ungapped_motifs)
	}

	return(list_ungapped_motifs)
}

generate.gapped.motifs.dictionnary <- function(parameters_bipartite, nucleotides){
	gapped_motif_left_min_length <- parameters_bipartite[["left_min_len"]]
	gapped_motif_left_max_length <- parameters_bipartite[["left_max_len"]]
	gapped_motif_gap_min_length <- parameters_bipartite[["gap_min_len"]]
	gapped_motif_gap_max_length <- parameters_bipartite[["gap_max_len"]]
	gapped_motif_right_min_length <- parameters_bipartite[["right_min_len"]]
	gapped_motif_right_max_length <- parameters_bipartite[["right_max_len"]]

	list_leftside_motifs <- generate.ungapped.motifs.dictionnary(gapped_motif_left_min_length, gapped_motif_left_max_length, nucleotides)
	list_gaps <- generate.ungapped.motifs.dictionnary(gapped_motif_gap_min_length, gapped_motif_gap_max_length, c("N"))
	list_rightside_motifs <- generate.ungapped.motifs.dictionnary(gapped_motif_right_min_length, gapped_motif_right_max_length, nucleotides)

	gapped_motifs_df <- expand.grid(list_leftside_motifs, list_gaps, list_rightside_motifs)
	list_gapped_motifs <- apply(gapped_motifs_df, 1, paste, collapse="") # Long exec time, TODO parSapply?

	return(list_gapped_motifs)
}

generate.motifs.dictionnary <- function(ungapped_motif_min_length, ungapped_motif_max_length, generate_bipartite=FALSE, parameters_bipartite=NA){
	# Global parameters
	nucleotides <- c("A","C","G","T")

	# Ungapped motifs
	list_ungapped_motifs <- generate.ungapped.motifs.dictionnary(ungapped_motif_min_length, ungapped_motif_max_length, nucleotides)

	if(generate_bipartite){
		# Gapped motifs
		list_gapped_motifs <- generate.gapped.motifs.dictionnary(parameters_bipartite, nucleotides)

		return(c(list_ungapped_motifs, list_gapped_motifs))
	}else{
		return(list_ungapped_motifs)
	}
}

decompose.degenerated.motifs <- function(list_motifs, iupac_nc, partial=TRUE){
	custom_iupac_nc <- iupac_nc %>%
		add_row(code="n",pattern="[ACGT]",choice="n") %>%
		mutate(choice_vector=str_split(choice,""))

	list_decomposed_motifs <- foreach(motif=list_motifs, .combine=c) %do% {
		length_motif <- nchar(motif)
		decomposed_motif <- as.data.frame(str_split(motif,"",simplify=TRUE))
		colnames(decomposed_motif) <- paste0("pos",seq(1,length_motif))

		if(partial %in% c(1,2)){
			is_N <- grepl("N",t(decomposed_motif))
			N_stretch <- rle(is_N)
			N_stretch_idx <- which(N_stretch$values==TRUE & N_stretch$lengths>=3)
			N_stretch_len_cumsum <- cumsum(N_stretch$lengths)
			N_stretch_starts <- N_stretch_len_cumsum[N_stretch_idx-1] + 1
			N_stretch_ends <- N_stretch_len_cumsum[N_stretch_idx]
			N_to_protect <- is_N
			if(partial==1){
				N_to_protect[c(N_stretch_starts,N_stretch_ends)] <- FALSE # not helpful in reality
			}
			decomposed_motif[,N_to_protect] <- "n"
		}

		motifs_base <- vector("list", length_motif)
		for(i in seq(1,length_motif)) {
			motifs_base[[i]] <- unlist(custom_iupac_nc$choice_vector[match(as.character(decomposed_motif[,i]),custom_iupac_nc$code)])
		}

		motifs_df <- expand.grid(motifs_base)
		decomposed_motifs <- apply(motifs_df, 1, paste, collapse="")

		if(partial){
			decomposed_motifs <- toupper(decomposed_motifs)
		}

		return(decomposed_motifs)
	}

	return(list_decomposed_motifs)
}

select.random.motifs <- function(list_motifs, nb_motifs, per_length=FALSE, resample=NA){
	df_motifs <- as.data.frame(list_motifs, stringsAsFactors=FALSE)
	colnames(df_motifs) <- "motif" # Replace default name "list_motifs"
	df_motifs$len_motif <- nchar(df_motifs$motif)

	if(per_length){
		list_random_motifs <- foreach(motif_len=unique(df_motifs$len_motif), .combine=c) %do% {
			subset_random_motifs <- sample_n(subset(df_motifs, len_motif==motif_len), size=nb_motifs, replace=FALSE)$motif

			return(subset_random_motifs)
		}
	}else if(any(!is.na(resample))){
		df_resample <- as.data.frame(resample, stringsAsFactors=FALSE)
		colnames(df_resample) <- "motif" # Replace default name "list_motifs"
		df_resample$len_motif <- nchar(df_resample$motif)
		df_resample <- df_resample %>%
			group_by(len_motif) %>%
			summarize(n=n(), .groups="drop_last")

		list_random_motifs <- foreach(idx=seq(1,nrow(df_resample)), .combine=c) %do% {
			subset_random_motifs <- sample_n(subset(df_motifs, len_motif==df_resample$len_motif[idx]), size=df_resample$n[idx], replace=FALSE)$motif

			return(subset_random_motifs)
		}
	}else{
		list_random_motifs <- sample_n(df_motifs, size=nb_motifs, replace=FALSE)$motif
	}

	return(list_random_motifs)
}

complete.decompose.motifs <- function(list_motifs, iupac_nc){
	custom_iupac_nc <- iupac_nc %>%
		mutate(choice_vector=str_split(choice,""))

	df_decomposed_motifs <- foreach(motif=list_motifs, .combine=rbind) %do% {
		length_motif <- nchar(motif)
		decomposed_motif <- as.data.frame(str_split(motif,"",simplify=TRUE))
		colnames(decomposed_motif) <- paste0("pos",seq(1,length_motif))

		motifs_base <- vector("list", length_motif)
		for(i in seq(1,length_motif)) {
			motifs_base[[i]] <- unlist(custom_iupac_nc$choice_vector[match(as.character(decomposed_motif[,i]),custom_iupac_nc$code)])
		}

		motifs_df <- expand.grid(motifs_base)
		decomposed_motifs <- apply(motifs_df, 1, paste, collapse="")

		return(data.frame(motifs=motif, decomposed_motifs=decomposed_motifs, n=length(decomposed_motifs)))
	}

	return(df_decomposed_motifs)
}

estimate.motifs.overlaps <- function(list_motifs_of_interest, reference_decomposed_motifs){
	results <- foreach(motif_of_interest=list_motifs_of_interest, .combine=rbind) %do% {
		matches <- grepl(convert.motif.grep(motif_of_interest), reference_decomposed_motifs$decomposed_motifs)
		reference_decomposed_motifs$matched <- ifelse(matches,"match","no_match")

		tmp_results <- reference_decomposed_motifs %>%
			group_by(motifs, n, matched) %>%
			summarize(n_matches=n(), .groups="drop_last") %>%
			spread(matched, n_matches, fill=0)
		
		# Add potential missing columns
		if(! c("match") %in% colnames(tmp_results)){
			tmp_results$match <- 0
		}
		if(! c("no_match") %in% colnames(tmp_results)){
			tmp_results$no_match <- 0
		}

		tmp_results <- tmp_results %>%
			dplyr::select(motifs, n, match, no_match) %>% # Reorder columns
			mutate(ratio=match/n)
		tmp_results$motif_of_interest <- motif_of_interest

		return(tmp_results)
	}

	return(results)
}

select.unmodified.random.motifs <- function(list_motifs_to_pick_from, list_expected_motifs, max_ratio,  nb_motifs, per_length=TRUE){
	list_final_random_motifs <- NULL
	flag_while <- TRUE
	flag_resample <- NA
	while(flag_while){
		list_random_motifs <- select.random.motifs(list_motifs_to_pick_from, nb_motifs, per_length, flag_resample)

		full_decomposed_random_motifs <- complete.decompose.motifs(list_random_motifs, iupac_nc)

		overlaps_expected <- estimate.motifs.overlaps(list_expected_motifs, full_decomposed_random_motifs)

		list_modified_motifs <- as.character(subset(overlaps_expected, ratio>max_ratio)$motifs)

		if(length(list_modified_motifs)>0){
			list_motifs_to_pick_from <- list_motifs_to_pick_from[!list_motifs_to_pick_from %in% list_modified_motifs] # Remove bad pick from list
			nb_motifs <- length(list_modified_motifs) # Will be ignored because non-NA flag_resample
			flag_resample <- list_modified_motifs # Indicate length of motifs to pick
			per_length <- FALSE # Disable to use resampling
			list_final_random_motifs <- c(list_final_random_motifs, list_random_motifs[!list_random_motifs %in% list_modified_motifs])
		}else{
			flag_while <- FALSE
		}
	}

	return(list_final_random_motifs)
}

# Generating list of motifs; Not used
# Ungapped
# ungapped_motif_min_length <- 4
# ungapped_motif_max_length <- 6
# Gapped
# parameters_bipartite <- list(
# 	left_min_len=3, left_max_len=4, # Have 2 and 5
# 	gap_min_len=5, gap_max_len=6, # Have 4 and 7
# 	right_min_len=3, right_max_len=4 # Have 2
# )

# list_motifs <- generate.motifs.dictionnary(ungapped_motif_min_length, ungapped_motif_max_length, TRUE, parameters_bipartite)
# saveRDS(list_motifs, file="list_motifs.RDS")


#   _____                   _    ___  ___      _   _  __     
#  /  __ \                 | |   |  \/  |     | | (_)/ _|    
#  | /  \/ ___  _   _ _ __ | |_  | .  . | ___ | |_ _| |_ ___ 
#  | |    / _ \| | | | '_ \| __| | |\/| |/ _ \| __| |  _/ __|
#  | \__/\ (_) | |_| | | | | |_  | |  | | (_) | |_| | | \__ \
#   \____/\___/ \__,_|_| |_|\__| \_|  |_/\___/ \__|_|_| |___/
#                                                            
#                                                            

prepare.splitted.contifs.annotation <- function(contigs_info, target_length=100000, keep_all=FALSE){
	contigs_info <- contigs_info %>%
		mutate(nb_split=round(length/target_length, 0))
	contigs_toSplit <- contigs_info %>%
		filter(nb_split>1)

	splitted_contigs <- foreach(contig_idx=seq(1,nrow(contigs_toSplit)), .combine=rbind) %do% {
		contig_toSplit <- contigs_toSplit[contig_idx,]
		contig_start <- contig_toSplit$start
		contig_end <- contig_toSplit$end
		nb_split <- contig_toSplit$nb_split

		# contig_breaks <- round(seq(contig_start, contig_end, length.out=nb_split+1), 0)
		contig_breaks <- round(seq(0, contig_end - (contig_start - 1), length.out=nb_split+1), 0)
		starts <- contig_breaks[seq(1,nb_split)]
		ends <- c(contig_breaks[seq(2,nb_split)] - 1, contig_breaks[nb_split+1])

		splitted_contig <- data.frame(
			genome=contig_toSplit$genome,
			start=starts,
			end=ends,
			contig=paste0(contig_toSplit$contig,"_",starts,":",ends), # was genome
			length=ends-(starts-1),
			id=contig_toSplit$id,
			parent_contig=str_split(contig_toSplit$contig," ",simplify=TRUE)[,1],
			split_id=seq(1,nb_split)
		)
		
		return(splitted_contig)
	}

	if(keep_all){
		contigs_toConserve <- contigs_info %>%
			filter(nb_split<=1) %>%
			dplyr::select(-nb_split)
		contigs_toConserve$split_id <- 1
		contigs_toConserve$parent_contig <- contigs_toConserve$contig
		splitted_contigs <- rbind(contigs_toConserve, splitted_contigs)
	}

	return(splitted_contigs)
}

split.motifs.long.contigs <- function(motifs, splitted_contigs_info){
	splitted_contigs_info <- subset(splitted_contigs_info, contig!=parent_contig) # Handle motifs & signal

	motifs_toSplit <- subset(motifs, contig_name %in% splitted_contigs_info$parent_contig)
	motifs_toConserve <- subset(motifs, !contig_name %in% splitted_contigs_info$parent_contig)
	motifs_splitted <- foreach(splitted_contigs_idx=seq(1,nrow(splitted_contigs_info)), .combine=rbind) %do% {
		splitted_contig_info <- splitted_contigs_info[splitted_contigs_idx,]

		subset_tmp <- subset(motifs_toSplit, contig_name==as.character(splitted_contig_info$parent_contig) & contig_pos_motif>=splitted_contig_info$start & contig_pos_motif<=splitted_contig_info$end)

		if(nrow(subset_tmp)>0){
			subset_tmp$contig_name <- splitted_contig_info$contig
		}

		return(subset_tmp)
	}
	motifs_splitted <- rbind(motifs_toConserve, motifs_splitted)

	return(motifs_splitted)
}

count.metagenome.motifs <- function(path_metagenome, genome_id, list_motifs_to_process, target_length, iupac_nc, nbCPU){
	left_signal <- -1 # Keep overlapping motifs, cannot be exactly same position because one by one
	right_signal <- -1 # Keep overlapping motifs
	error_margin <- -1 # Keep overlapping motifs
	arbitrary_mod_pos <- 1

	contigs_info <- create.contigs.info(path_metagenome, genome_id, "mock")
	if(!is.na(target_length)){
		contigs_info <- prepare.splitted.contifs.annotation(contigs_info, target_length, TRUE)
	}

	registerDoMC(nbCPU) # Silent error if too high
	number_motifs <- foreach(current_motif=list_motifs_to_process, .combine=rbind.fill, .multicombine=TRUE, .maxcombine=100) %dopar% { # TODO look multiple motifs
		motif_length <- nchar(current_motif)
		motif_summary <- data.frame(motif=current_motif, mod_pos=arbitrary_mod_pos, stringsAsFactors=FALSE)
		motifs <- find.isolated.motifs(path_metagenome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, 1, FALSE) # Long

		if(!is.na(target_length)){
			motifs <- split.motifs.long.contigs(motifs, contigs_info)
		}

		number_motif <- motifs %>%
			group_by(contig_name, motif) %>%
			summarize(n=n(), .groups="drop_last")

		return(number_motif)
	}
	registerDoSEQ()

	annotated_number_motifs <- merge(number_motifs, contigs_info, by.x=c("contig_name"), by.y=c("contig"))

	return(annotated_number_motifs)
}

heatmap.motifs.distribution <- function(distribution_motifs, contigs_info, base_name){
	if("n" %in% colnames(distribution_motifs)){
		gp <- ggplot(distribution_motifs) +
			geom_tile(aes(x=motif, y=reorder(contig_name, length), fill=n)) +
			scale_fill_distiller(palette="Spectral", trans='log10') +
			facet_grid(id~., space="free_y", scales="free_y") +
			theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			theme(axis.text.x=element_text(angle=45, hjust=1)) +
			labs(x="Motifs", y="Contigs", fill="Number of\nOccurrences")
	}else if("average_nb_occurrence" %in% colnames(distribution_motifs)){
		gp <- ggplot(merge(distribution_motifs, contigs_info)) +
			geom_tile(aes(x=motif, y=reorder(contig, length), fill=average_nb_occurrence)) +
			scale_fill_distiller(palette="Spectral", trans='log10') +
			facet_grid(id~., space="free_y", scales="free_y") +
			theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			theme(axis.text.x=element_text(angle=45, hjust=1)) +
			labs(x="Motifs", y="Contigs", fill="Average\nnumber of\nOccurrences")
	}else if("nb_occurrence" %in% colnames(distribution_motifs)){
		gp <- ggplot(merge(distribution_motifs, contigs_info)) +
			geom_tile(aes(x=motif, y=reorder(contig, length), fill=nb_occurrence)) +
			scale_fill_distiller(palette="Spectral", trans='log10') +
			facet_grid(id~., space="free_y", scales="free_y") +
			theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			theme(axis.text.x=element_text(angle=45, hjust=1)) +
			labs(x="Motifs", y="Contigs", fill="Average\nnumber of\nOccurrences")
	}
	pdf(paste0("Distribution_motifs.",base_name,".v1.pdf"), height=30, width=16)
	print(gp)
	dev.off()
}

table.motifs.distribution <- function(motifs_info, distribution_motifs, id_genome, length_threshold){
	res <- distribution_motifs %>%
		filter(id==id_genome & motif %in% motifs_info$motif[which(motifs_info$id==id_genome)]) %>%
		spread(motif,n) %>%
		filter(length>=length_threshold) %>%
		arrange(length)

	return(res)
}

#   _____                     ___  ___      _   _  __     
#  /  ___|                    |  \/  |     | | (_)/ _|    
#  \ `--.  ___ ___  _ __ ___  | .  . | ___ | |_ _| |_ ___ 
#   `--. \/ __/ _ \| '__/ _ \ | |\/| |/ _ \| __| |  _/ __|
#  /\__/ / (_| (_) | | |  __/ | |  | | (_) | |_| | | \__ \
#  \____/ \___\___/|_|  \___| \_|  |_/\___/ \__|_|_| |___/
#                                                         
#                                                         

split.signal.long.contigs <- function(contigs_methylation_signal, contigs_info, target_length){
	splitted_contigs_info <- prepare.splitted.contifs.annotation(contigs_info, target_length, FALSE)

	contigs_methylation_signal_toSplit <- subset(contigs_methylation_signal, contig %in% splitted_contigs_info$parent_contig)
	contigs_methylation_signal_toConserve <- subset(contigs_methylation_signal, !contig %in% splitted_contigs_info$parent_contig)
	contigs_methylation_signal_splitted <- foreach(splitted_contigs_idx=seq(1,nrow(splitted_contigs_info)), .combine=rbind) %do% {
		splitted_contig_info <- splitted_contigs_info[splitted_contigs_idx,]

		subset_tmp <- subset(contigs_methylation_signal_toSplit, contig==as.character(splitted_contig_info$parent_contig) & contig_pos_motif>=splitted_contig_info$start & contig_pos_motif<=splitted_contig_info$end)

		if(nrow(subset_tmp)>0){
			subset_tmp$contig <- splitted_contig_info$contig
		}

		return(subset_tmp)
	}
	contigs_methylation_signal_splitted <- rbind(contigs_methylation_signal_toConserve, contigs_methylation_signal_splitted)

	return(contigs_methylation_signal_splitted)
}

make.progress.bar <- function(list_motifs_to_process, nbCPU){
	pb <- progress_bar$new(format=" Motifs processed (:detail): [:bar] :percent eta: :eta (elapsed: :elapsedfull)", total=length(list_motifs_to_process)/nbCPU, show_after=0)

	return(pb)
}

initialize.progress.bar <- function(pb, list_motifs_to_process){
	detail <- paste0(0,"/",length(list_motifs_to_process))
	pb$tick(0, tokens=list(detail=detail))
}

progress.tracker <- function(pb, current_motif, list_motifs_to_process, nbCPU){
	motif_idx <- match(current_motif,list_motifs_to_process)
	if(motif_idx %% nbCPU == 0){
		detail <- paste0(motif_idx,"/",length(list_motifs_to_process))
		pb$tick(tokens=list(detail=detail))
	}
}

terminate.progress.bar <- function(pb, list_motifs_to_process, nbCPU){
	nb_motifs <- length(list_motifs_to_process)
	detail <- paste0(nb_motifs,"/",nb_motifs)
	pb$tick(nb_motifs/nbCPU, tokens=list(detail=detail))
}

score.metagenome.motifs <- function(genome_id, path_metagenome, data_type, list_motifs_to_process, score_type, methylation_signal, gr_methylated_motifs_signal, min_cov, signature_state, target_length, iupac_nc, nbCPU){
	left_signal <- -1 # Keep overlapping motifs, cannot be exactly same position because one by one
	right_signal <- -1 # Keep overlapping motifs
	error_margin <- -1 # Keep overlapping motifs
	arbitrary_mod_pos <- 1
	expected_signal_left <- -6
	expected_signal_right <- -1
	signal_margin <- 0

	# Prepare contigs information
	contigs_info <- create.contigs.info(path_metagenome, genome_id, data_type)

	print_message("    Initialize methylation feature computation")

	# Prepare methylation GRanges
	gr_methylation <- GRanges(
		seqnames=methylation_signal$contig,
		ranges=IRanges(methylation_signal$position, methylation_signal$position),
		strand=ifelse(methylation_signal$dir=="fwd","+","-")
	)

	# Mark modified position
	methylation_signal$signal_stat <- factor("no_signal", levels=c("no_signal", "as_signal"))
	if(length(gr_methylated_motifs_signal)>1){ # If NA then = 1
		overlaps_tmp <- findOverlaps(gr_methylation, gr_methylated_motifs_signal, type="any", select="all")
		methylation_signal$signal_stat[overlaps_tmp@from] <- as.factor("as_signal")
	}

	pb <- make.progress.bar(list_motifs_to_process, nbCPU)
	print_message("    Processing motifs")
	initialize.progress.bar(pb, list_motifs_to_process)
	registerDoMC(nbCPU) # Silent error if too high
	scored_motifs <- foreach(current_motif=list_motifs_to_process, .combine=rbind.fill, .multicombine=TRUE, .maxcombine=100) %dopar%{
		motif_length <- nchar(current_motif)
		motif_summary <- data.frame(motif=current_motif, mod_pos=arbitrary_mod_pos, stringsAsFactors=FALSE)
		motifs <- find.isolated.motifs(path_metagenome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, 1, FALSE) # Long

		if(nrow(motifs)==0){
			# No matching motif in all contigs
			empty_scored_motif <- data.frame(contig="None", motif=current_motif, distance_motif=0, signal_ratio=NA, dist_score=0.0, nb_occurrence=0)
			empty_scored_motif$nb_occurrence <- as.integer(empty_scored_motif$nb_occurrence)

			return(empty_scored_motif)
		}

		expected_motifs_signal <- motifs %>%
			mutate(strand=as.factor(ifelse(dir=="fwd","+","-"))) %>%
			mutate(left_side=contig_pos_motif + expected_signal_left - signal_margin) %>%
			mutate(right_side=contig_pos_motif + expected_signal_right + signal_margin) %>%
			mutate(left_side=ifelse(strand=="+", left_side - (arbitrary_mod_pos - 1), left_side - (motif_length - arbitrary_mod_pos))) %>% # Add full motif window with strand correction
			mutate(right_side=ifelse(strand=="+", right_side + (motif_length - arbitrary_mod_pos), right_side + (arbitrary_mod_pos - 1)))

		gr_expected_motifs_signal <- GRanges(
			seqnames=expected_motifs_signal$contig_name,
			ranges=IRanges(expected_motifs_signal$left_side, expected_motifs_signal$right_side),
			strand=expected_motifs_signal$strand
		)

		overlaps_motifs_methylation <- findOverlaps(gr_expected_motifs_signal, gr_methylation, type="any", select="all")
		matched_expected_motifs_signal <- subset(expected_motifs_signal[overlaps_motifs_methylation@from,], select=-c(contig_name,dir,strand)) # TODO could remove more
		matched_methylation_signal <- methylation_signal[overlaps_motifs_methylation@to,]
		contigs_methylation_signal <- cbind(matched_expected_motifs_signal, matched_methylation_signal)

		if(!is.na(target_length)){ # Split long contigs if asked for
			contigs_methylation_signal <- split.signal.long.contigs(contigs_methylation_signal, contigs_info, target_length)
		}

		# Preparing data and confidence filtering
		if(all(c("N_wga","N_nat") %in% colnames(contigs_methylation_signal))){ # natVSwga dataset
			scored_motif <- contigs_methylation_signal %>% 
				mutate(distance_motif=ifelse(dir=="fwd",position-contig_pos_motif,(-(position-contig_pos_motif)) - 7)) %>% # Relative distance to mod_pos with strand correction
				filter(N_wga>=min_cov & N_nat>=min_cov)
		}else if(c("N_val") %in% colnames(contigs_methylation_signal)){ # in silico dataset
			scored_motif <- contigs_methylation_signal %>% 
				mutate(distance_motif=ifelse(dir=="fwd",position-contig_pos_motif,(-(position-contig_pos_motif)) - 7)) %>% # Relative distance to mod_pos with strand correction
				filter(N_val>=min_cov)
		}

		if(signature_state=="complete"){ # Keep point with complete signature
			scored_motif <- scored_motif %>%
				group_by(contig, contig_pos_motif, dir, motif) %>%
				filter(n()==(expected_signal_right - expected_signal_left) + 1 + signal_margin*2 + (motif_length - 1)) # Filtered out incomplete occurrences
		}
		scored_motif <- scored_motif %>%
			group_by(contig, motif, distance_motif, signal_stat) %>%
			summarize(means_diff=list(mean_diff), .groups="drop_last")

		# Compute signal ratio
		if(length(gr_methylated_motifs_signal)>1){ # If NA then = 1
			list_signal_stat <- unique(scored_motif$signal_stat)
			if(all(c("no_signal","as_signal") %in% list_signal_stat)){
				scored_motif <- scored_motif %>%
					spread(signal_stat, means_diff) %>%
					mutate(signal_ratio=length(unlist(as_signal))/(length(unlist(as_signal)) + length(unlist(no_signal)))) %>%
					mutate(means_diff=list(c(unlist(as_signal),unlist(no_signal)))) %>%
					dplyr::select(-c(as_signal, no_signal))
			}else if(!c("no_signal") %in% list_signal_stat){
				scored_motif <- scored_motif %>%
					mutate(signal_ratio=1) %>%
					dplyr::select(-c(signal_stat))
			}else if(!c("as_signal") %in% list_signal_stat){
				scored_motif <- scored_motif %>%
					mutate(signal_ratio=0) %>%
					dplyr::select(-c(signal_stat))
			}
		}else{
			scored_motif <- scored_motif %>%
				mutate(signal_ratio=NA) %>% # Keep for score_type !dist
				dplyr::select(-c(signal_stat))
		}

		# Score relative position
		if(grepl("abs",score_type)){
			scored_motif <- scored_motif %>%
				group_by(contig, motif, distance_motif) %>%
				mutate(dist_score=mean(abs(unlist(means_diff)), na.rm=TRUE), nb_occurrence=length(unlist(means_diff))) %>%
				dplyr::select(-c(means_diff))
		}else if(grepl("real",score_type)){
			scored_motif <- scored_motif %>%
				group_by(contig, motif, distance_motif) %>%
				mutate(dist_score=mean(unlist(means_diff), na.rm=TRUE), nb_occurrence=length(unlist(means_diff))) %>%
				dplyr::select(-c(means_diff))
		}

		# Combine score per motif
		if(!grepl("dist",score_type)){
			# TODO maybe abs if some case???
			scored_motif <- scored_motif %>% # nb_occurrence accurate only if complete signature
				group_by(contig, motif) %>%
				top_n(n=6, wt=dist_score) %>% # Keep only 6 max position to be comparable between motifs
				arrange(contig, motif, desc(dist_score)) %>%
				group_by(contig, motif) %>%
				filter(row_number()<=6) %>% # Remove potential ties; Rare
				summarize(score=sum(dist_score), average_nb_occurrence=mean(nb_occurrence), average_signal_ratio=mean(signal_ratio), .groups="drop_last") # average_signal_ratio could recompute
		}

		progress.tracker(pb, current_motif, list_motifs_to_process, nbCPU)

		return(scored_motif)
	}
	registerDoSEQ()
	terminate.progress.bar(pb, list_motifs_to_process, nbCPU)

	nb_processed_motifs <- unique(as.character(scored_motifs$motif))
	if(length(nb_processed_motifs)!=length(list_motifs_to_process)){
		print_message("Possible error during scoring.") # Silent error appended when out of memory
	}

	if("None" %in% levels(scored_motifs$contig)){
		scored_motifs <- subset(scored_motifs, contig!="None") # Remove potential entry from missing motifs (no occurrence across all contigs)
		scored_motifs$contig <- droplevels(scored_motifs$contig)
		scored_motifs$motif <- droplevels(scored_motifs$motif)
	}

	return(scored_motifs)
}

mark.methylated.position <- function(path_metagenome, motifs_info, iupac_nc, nbCPU){
	left_signal <- -1 # Keep overlapping motifs, cannot be exactly same position because one by one
	right_signal <- -1 # Keep overlapping motifs
	error_margin <- -1 # Keep overlapping motifs
	expected_signal_left <- -6
	expected_signal_right <- -1
	signal_margin <- 0
	nbCPU <- 10

	gr_methylated_motifs_signal <- foreach(genome_id=unique(motifs_info$id), .combine=c) %do% {
		methyated_motifs <- subset(motifs_info, id==genome_id)

		motifs_metagenome <- find.isolated.motifs(path_metagenome, methyated_motifs, iupac_nc, left_signal, right_signal, error_margin, nbCPU, FALSE) # Long
		methylated_motifs_contigs <- motifs_metagenome[motifs_metagenome$contig_name %in% subset(contigs_info, id==genome_id)$contig,]

		methylated_motifs_signal_contigs <- methylated_motifs_contigs %>%
			mutate(strand=as.factor(ifelse(dir=="fwd","+","-"))) %>%
			mutate(left_side=contig_pos_motif + expected_signal_left - signal_margin) %>%
			mutate(right_side=contig_pos_motif + expected_signal_right + signal_margin)

		gr_methylated_motifs_signal_contigs <- GRanges(
			seqnames=methylated_motifs_signal_contigs$contig_name,
			ranges=IRanges(methylated_motifs_signal_contigs$left_side, methylated_motifs_signal_contigs$right_side),
			strand=methylated_motifs_signal_contigs$strand
		)

		return(gr_methylated_motifs_signal_contigs)
	}

	return(gr_methylated_motifs_signal)
}

annotate.score <- function(scored_motifs, contigs_info){
	annotated_scored_motifs <- merge(scored_motifs, contigs_info)

	return(annotated_scored_motifs)
}

plot.motifs.dist.score <- function(dist_scored_motifs, contigs_info, occ_threshold, score_threshold){
	contigs_info$contig <- str_split(contigs_info$contig," ",simplify=TRUE)[,1] # NEW
	annotated_scored_motifs <- annotate.score(dist_scored_motifs, contigs_info)

	gp <- ggplot(annotated_scored_motifs)
	if(all(is.na(annotated_scored_motifs$signal_ratio))){ # Don't know for non-mock
		gp <- gp +
			geom_point(aes(x=dist_score, y=log10(nb_occurrence)), pch=46)
	}else{
		gp <- gp +
			geom_point(aes(x=dist_score, y=log10(nb_occurrence), col=signal_ratio), pch=46)
	}
	if("genome" %in% colnames(annotated_scored_motifs)){ # TODO could replace by any small number id
		gp <- gp +
			facet_wrap(~genome)
	}
	gp <- gp +
		geom_hline(yintercept=log10(occ_threshold)) +
		geom_vline(xintercept=score_threshold) +
		labs(pch="Type of motif")

	print(gp)
}

plot.motifs.score <- function(scored_motifs, contigs_info, occ_threshold){
	annotated_scored_motifs <- annotate.score(scored_motifs, contigs_info)
	
	gp <- ggplot(annotated_scored_motifs) +
		geom_point(aes(x=score, y=log10(average_nb_occurrence), col=average_signal_ratio), pch=46) + # TODO add color if should be modified
		geom_hline(yintercept=log10(occ_threshold)) +
		facet_wrap(~genome) +
		labs(pch="Type of motif")
	print(gp)
}

convert.scored.motifs.to.matrix <- function(scored_motifs, occ_threshold, col_spread=NA, col_noSpread=NA, fill_missing=NA){
	if(is.na(col_spread) & is.na(col_noSpread)){
		if("dist_score" %in% colnames(scored_motifs)){ # With distance_motif
			matrix_scored_motifs <- scored_motifs %>%
				filter(nb_occurrence>=occ_threshold) %>% # Entirely remove some contigs
				mutate(motif=paste0(motif,"_",distance_motif)) %>% # Integrate relative position in motif name
				dplyr::select(-distance_motif, -nb_occurrence, -signal_ratio) %>%
				spread(motif, dist_score, fill=fill_missing)
		}else{ # Without distance_motif
			matrix_scored_motifs <- scored_motifs %>%
				filter(average_nb_occurrence>=occ_threshold) %>% # Could remove some contigs
				dplyr::select(-average_nb_occurrence, -average_signal_ratio) %>%
				spread(motif, score, fill=fill_missing)
		}
	}else{
		matrix_scored_motifs <- scored_motifs %>%
			filter(nb_occurrence>=occ_threshold) %>% # Entirely remove some contigs
			mutate(motif=paste0(motif,"_",distance_motif)) %>% # Integrate relative position in motif name
			dplyr::select(-one_of(c("distance_motif", "nb_occurrence", col_noSpread))) %>%
			spread_("motif", col_spread, fill=fill_missing)
	}

	return(matrix_scored_motifs)
}

#   _____ _           _            _             
#  /  __ \ |         | |          (_)            
#  | /  \/ |_   _ ___| |_ ___ _ __ _ _ __   __ _ 
#  | |   | | | | / __| __/ _ \ '__| | '_ \ / _` |
#  | \__/\ | |_| \__ \ ||  __/ |  | | | | | (_| |
#   \____/_|\__,_|___/\__\___|_|  |_|_| |_|\__, |
#                                           __/ |
#                                          |___/ 

replace.na.by.func <- function(x, type){
	if(type=="min"){
		x[is.na(x)] <- min(x, na.rm=TRUE)
	}else if(type=="mean"){
		x[is.na(x)] <- mean(x, na.rm=TRUE)
	}else if(type=="median"){
		x[is.na(x)] <- median(x, na.rm=TRUE)
	}else if(type=="pseudoZero"){
		x[is.na(x)] <- runif(length(x[is.na(x)]), -0.2, 0.2)
	}

	return(x)
}

prepare.contig.coverage <- function(input_coverage_A, input_coverage_B){
	contig_annotation_A <- process.coverage(input_coverage_A)
	contig_annotation_B <- process.coverage(input_coverage_B)

	contig_coverage <- merge(contig_annotation_A, contig_annotation_B, by=c("chr","contig_length"), suffixes=c(".dataset_A",".dataset_B")) %>%
		mutate(diff=abs(avg_cov.dataset_A-avg_cov.dataset_B)) %>%
		mutate(ratio=abs((avg_cov.dataset_A + 0.001)/(avg_cov.dataset_B + 0.001))) %>%
		arrange(desc(diff))

	return(contig_coverage)
}

tsne.motifs.score.dev <- function(path_metagenome, scored_motifs, data_type, occ_threshold, length_threshold, handle_NA, genome_id, grouping_var, contig_weight_unit, contig_coverage=NA, max_relative_weight=0.05, tsne_perplexity=30, tsne_max_iter=2500, rdm_seed=42, tsne_seed=101){
	print_message("    Prepare methylation profile matrix")
	# Convert score data.frame to matrix for t-SNE
	# Handle missing scores
	if(handle_NA==0){ # TODO is possible within convert.*
		matrix_scored_motifs <- convert.scored.motifs.to.matrix(scored_motifs, occ_threshold, NA, NA, handle_NA)
	}else{
		matrix_scored_motifs <- convert.scored.motifs.to.matrix(scored_motifs, occ_threshold)
		tmp_matrix_scored_motifs <- matrix_scored_motifs[,!(names(matrix_scored_motifs) %in% c("contig"))]
		if(handle_NA=="pseudoZero"){set.seed(rdm_seed)} # Keep tSNE reproducible
		tmp_matrix_scored_motifs <- apply(tmp_matrix_scored_motifs, 2, replace.na.by.func, type=handle_NA) # Replace NA per function results TODO is possible within convert.*
		matrix_scored_motifs[,!(names(matrix_scored_motifs) %in% c("contig"))] <- tmp_matrix_scored_motifs
	}

	# Prepare contigs information
	contigs_info <- create.contigs.info(path_metagenome, genome_id, data_type)
	if(contig_weight_unit %in% c(-1,-2) & data_type=="mock"){
		splitted_contigs_info <- colsplit(gsub("(.*)_([0-9]+):([0-9]+)","\\1 \\2 \\3",matrix_scored_motifs$contig), " ", names=c("genome","start","end")) %>%
			mutate(contig=as.character(matrix_scored_motifs$contig)) %>%
			filter(as.character(genome)!=as.character(start)) %>% # Keep splitted contig only
			mutate(start=as.numeric(as.character(start)), end=as.numeric(as.character(end))) %>%
			mutate(length=end-start + 1) %>%
			dplyr::select(genome, contig, length, start, end)
		# splitted_contigs_info <- merge(splitted_contigs_info, subset(contigs_info, select=c(contig, id)), by="genome", all.x=TRUE)
		splitted_contigs_info <- merge(splitted_contigs_info, subset(contigs_info, select=c(contig, id)), by.x="genome", by.y="contig", all.x=TRUE)
		# splitted_contigs_info <- merge(splitted_contigs_info, subset(contigs_info, select=c(contig, id)), by="contig", all.x=TRUE)
		# splitted_contigs_info <- splitted_contigs_info[!duplicated(splitted_contigs_info),] # Not needed.
		splitted_contigs_info <- subset(splitted_contigs_info, !is.na(id))

		unsplitted_contigs_info <- contigs_info %>%
			mutate(contig=str_split(contig, " ", simplify=TRUE)[,1]) %>%
			dplyr::select(colnames(splitted_contigs_info)) %>%
			filter(!contig %in% splitted_contigs_info$contig) # was genome
		contigs_info <- rbind(splitted_contigs_info, unsplitted_contigs_info) %>%
			mutate(contig=as.factor(contig))

		# contigs_info <- contigs_info[!duplicated(contigs_info),] # NEW not sure
	}

	# Filter out small contigs
	if(!is.na(length_threshold)){
		if(data_type %in% c("mock","mock_denovo")){
			contigs_info <- contigs_info %>% #NEW
				filter(length>=length_threshold) %>% # Remove contig below length threshold
				filter(contig %in% matrix_scored_motifs$contig) # Remove contig without scored motifs (e.g. nb occ < threshold)
		}else if(data_type=="real"){
			contigs_info <- contigs_info %>% #NEW
				filter(length>=length_threshold) %>% # Remove contig below length threshold
				filter(genome %in% matrix_scored_motifs$contig) # Remove contig without scored motifs (e.g. nb occ < threshold)
		}
		if(any(contig_weight_unit %in% c(-1,-2) | !is.na(contig_weight_unit), na.rm=TRUE)){ # Contigs were splitted, used original length in threshold filtering
			if(data_type %in% c("mock","mock_denovo")){
				matrix_scored_motifs <- subset(matrix_scored_motifs, contig %in% contigs_info$contig)
			}else if(data_type=="real"){
				matrix_scored_motifs <- subset(matrix_scored_motifs, gsub("(.*)_[0-9]+:[0-9]+","\\1",contig) %in% str_split(contigs_info$contig, " ", simplify=TRUE)[,1][contigs_info$length>=length_threshold])
			}
		}else{
			matrix_scored_motifs <- subset(matrix_scored_motifs, contig %in% str_split(contigs_info$contig, " ", simplify=TRUE)[,1][contigs_info$length>=length_threshold])
		}
		matrix_scored_motifs$contig <- droplevels(matrix_scored_motifs$contig)
	}

	# Weight contigs if requested
	if(!is.na(contig_weight_unit) & !contig_weight_unit %in% c(-1,-2)){ # Don't do with splitted contigs
		nb_remaining_contigs <- nrow(contigs_info)
		max_weight <- min(25, max(1, round(nb_remaining_contigs * max_relative_weight))) # Add ceiling to weights

		contigs_info <- contigs_info %>%
			mutate(weight=ceiling(length/contig_weight_unit)) %>%
			mutate(weight=ifelse(weight > max_weight, max_weight, weight))

		if(!is.null(nrow(contig_coverage))){ # Reduce weight of poorly covered contigs
			if(data_type %in% c("mock","mock_denovo")){
				contigs_info <- merge(contigs_info, subset(contig_coverage, select=-c(contig_length)), by.x=c("contig"), by.y=c("chr"))
			}else if(data_type=="real"){
				contigs_info <- merge(contigs_info, subset(contig_coverage, select=-c(contig_length)), by.x=c("genome"), by.y=c("chr"))
			}
			contigs_info <- contigs_info %>%
				mutate(weight=ifelse(avg_cov.dataset_A<10 | avg_cov.dataset_B<10,1,weight))
		}

		element_to_repeat <- match(rep(str_split(contigs_info$contig, " ", simplify=TRUE)[,1], times=contigs_info$weight), matrix_scored_motifs$contig)
		matrix_scored_motifs <- matrix_scored_motifs[element_to_repeat[!is.na(element_to_repeat)],] # Some contigs are missing after filter(nb_occurrence>=5)
	}

	tsne_matrix <- as.matrix(subset(matrix_scored_motifs, select=-c(contig)))
	print_message("    Dimentionality reduction")
	set.seed(tsne_seed)
	tsne_data <- as.data.frame(Rtsne(tsne_matrix, check_duplicates=FALSE, perplexity=tsne_perplexity, max_iter=tsne_max_iter)$Y)

	tsne_data <- tsne_data %>%
		mutate(contig=matrix_scored_motifs$contig) %>% # Link back contig names
		plyr::rename(c('V1'='tSNE_1','V2'='tSNE_2')) %>% # Change colnames
		mutate(contig_length=contigs_info$length[match(contig, str_split(contigs_info$contig, " ", simplify=TRUE)[,1])]) %>%
		mutate(id=paste0(contigs_info$id[match(contig, str_split(contigs_info$contig, " ", simplify=TRUE)[,1])]))

	# Summarize weighted contigs for t-SNE
	if(!is.na(contig_weight_unit) & contig_weight_unit!=-1){ # Don't do with splitted -1 contigs 
		if(contig_weight_unit!=-2){
			print_message("    Finding weigthed contigs centroids")
		}else if(contig_weight_unit==-2){
			print_message("    Finding splitted contigs centroids")
			contigs_info_original <- create.contigs.info(path_metagenome, genome_id, data_type)
			unsplitted_tsne_data <- subset(tsne_data, contig %in% contigs_info_original$contig)
			splitted_tsne_data <- subset(tsne_data, !contig %in% contigs_info_original$contig)

			contigs_info_split <- colsplit(gsub("(.*)_([0-9]+):([0-9]+)","\\1 \\2 \\3",splitted_tsne_data$contig), " ", names=c("genome","start","end"))
			splitted_tsne_data <- cbind(splitted_tsne_data, contigs_info_split)
			splitted_tsne_data$contig <- as.character(splitted_tsne_data$contig) # Generate NA if unexpected contig
			stifle <- foreach(row_idx=seq(1,nrow(splitted_tsne_data))) %do% {
				matching_contigs_info_original <- subset(contigs_info_original, contig==as.character(splitted_tsne_data$genome[row_idx]))
				splitted_tsne_data$contig[row_idx] <- matching_contigs_info_original$contig # Is already in genome
				splitted_tsne_data$contig_length[row_idx] <- matching_contigs_info_original$length

				return(NA)
			}
			splitted_tsne_data$contig <- as.factor(splitted_tsne_data$contig)
			tsne_data <- rbind(subset(splitted_tsne_data, select=-c(genome, start, end)), unsplitted_tsne_data)
		}
		tsne_data <- tsne_data %>%
			group_by(contig) %>%
			mutate(centroid_V1=mean(tSNE_1),centroid_V2=mean(tSNE_2)) %>%
			distinct(contig, .keep_all=TRUE) %>%
			mutate(tSNE_1=centroid_V1,tSNE_2=centroid_V2) %>%
			dplyr::select(-c(centroid_V1, centroid_V2))
	}

	tsne_data <- tsne_data %>%
		group_by(.dots=grouping_var) %>%
		dplyr::mutate(n_contig=n()) %>%
		ungroup() %>%
		mutate(id=paste0(id," (n=",n_contig,")"))

	return(tsne_data)
}

tsne.motifs.score <- function(path_metagenome, scored_motifs, data_type, occ_threshold, length_threshold, handle_NA, genome_id, grouping_var, contig_weight_unit, contig_coverage=NA, max_relative_weight=0.05, tsne_perplexity=30, tsne_max_iter=2500, rdm_seed=42, tsne_seed=101){
	print_message("    Prepare methylation profile matrix")
	# Convert score data.frame to matrix for t-SNE
	# Handle missing scores
	if(handle_NA==0){ # TODO is possible within convert.*
		matrix_scored_motifs <- convert.scored.motifs.to.matrix(scored_motifs, occ_threshold, NA, NA, handle_NA)
	}else{
		matrix_scored_motifs <- convert.scored.motifs.to.matrix(scored_motifs, occ_threshold)
		tmp_matrix_scored_motifs <- matrix_scored_motifs[,!(names(matrix_scored_motifs) %in% c("contig"))]
		if(handle_NA=="pseudoZero"){set.seed(rdm_seed)} # Keep tSNE reproducible
		tmp_matrix_scored_motifs <- apply(tmp_matrix_scored_motifs, 2, replace.na.by.func, type=handle_NA) # Replace NA per function results TODO is possible within convert.*
		matrix_scored_motifs[,!(names(matrix_scored_motifs) %in% c("contig"))] <- tmp_matrix_scored_motifs
	}

	# Prepare contigs information
	contigs_info <- create.contigs.info(path_metagenome, genome_id, data_type)
	if(contig_weight_unit %in% c(-1,-2) & data_type=="mock"){
		splitted_contigs_info <- colsplit(gsub("(.*)_([0-9]+):([0-9]+)","\\1 \\2 \\3",matrix_scored_motifs$contig), " ", names=c("genome","start","end")) %>%
			mutate(contig=as.character(matrix_scored_motifs$contig)) %>%
			filter(as.character(genome)!=as.character(start)) %>% # Keep splitted contig only
			mutate(start=as.numeric(as.character(start)), end=as.numeric(as.character(end))) %>%
			mutate(length=end-start + 1) %>%
			dplyr::select(genome, contig, length, start, end)
		splitted_contigs_info <- merge(splitted_contigs_info, subset(contigs_info, select=c(contig, id)), by.x="genome", by.y="contig", all.x=TRUE)
		splitted_contigs_info <- subset(splitted_contigs_info, !is.na(id))

		unsplitted_contigs_info <- contigs_info %>%
			mutate(contig=str_split(contig, " ", simplify=TRUE)[,1]) %>%
			dplyr::select(colnames(splitted_contigs_info)) %>%
			filter(!contig %in% splitted_contigs_info$contig) # was genome
		contigs_info <- rbind(splitted_contigs_info, unsplitted_contigs_info) %>%
			mutate(contig=as.factor(contig))
	}

	# Filter out small contigs
	if(!is.na(length_threshold)){
		if(data_type %in% c("mock","mock_denovo")){
			contigs_info <- contigs_info %>% #NEW
				filter(length>=length_threshold) %>% # Remove contig below length threshold
				filter(contig %in% matrix_scored_motifs$contig) # Remove contig without scored motifs (e.g. nb occ < threshold)
		}else if(data_type=="real"){
			contigs_info <- contigs_info %>% #NEW
				filter(length>=length_threshold) %>% # Remove contig below length threshold
				filter(genome %in% matrix_scored_motifs$contig) # Remove contig without scored motifs (e.g. nb occ < threshold)
		}
		if(any(contig_weight_unit %in% c(-1,-2) | !is.na(contig_weight_unit), na.rm=TRUE)){ # Contigs were splitted, used original length in threshold filtering
			if(data_type %in% c("mock","mock_denovo")){
				matrix_scored_motifs <- subset(matrix_scored_motifs, contig %in% contigs_info$contig)
			}else if(data_type=="real"){
				matrix_scored_motifs <- subset(matrix_scored_motifs, gsub("(.*)_[0-9]+:[0-9]+","\\1",contig) %in% str_split(contigs_info$contig, " ", simplify=TRUE)[,1][contigs_info$length>=length_threshold])
			}
		}else{
			matrix_scored_motifs <- subset(matrix_scored_motifs, contig %in% str_split(contigs_info$contig, " ", simplify=TRUE)[,1][contigs_info$length>=length_threshold])
		}
		matrix_scored_motifs$contig <- droplevels(matrix_scored_motifs$contig)
	}

	# Weight contigs if requested
	if(!is.na(contig_weight_unit) & !contig_weight_unit %in% c(-1,-2)){ # Don't do with splitted contigs
		nb_remaining_contigs <- nrow(contigs_info)
		max_weight <- min(25, max(1, round(nb_remaining_contigs * max_relative_weight))) # Add ceiling to weights

		contigs_info <- contigs_info %>%
			mutate(weight=ceiling(length/contig_weight_unit)) %>%
			mutate(weight=ifelse(weight > max_weight, max_weight, weight))

		if(!is.null(nrow(contig_coverage))){ # Reduce weight of poorly covered contigs
			if(data_type %in% c("mock","mock_denovo")){
				contigs_info <- merge(contigs_info, subset(contig_coverage, select=-c(contig_length)), by.x=c("contig"), by.y=c("chr"))
			}else if(data_type=="real"){
				contigs_info <- merge(contigs_info, subset(contig_coverage, select=-c(contig_length)), by.x=c("genome"), by.y=c("chr"))
			}
			contigs_info <- contigs_info %>%
				mutate(weight=ifelse(avg_cov.dataset_A<10 | avg_cov.dataset_B<10,1,weight))
		}

		element_to_repeat <- match(rep(str_split(contigs_info$contig, " ", simplify=TRUE)[,1], times=contigs_info$weight), matrix_scored_motifs$contig)
		matrix_scored_motifs <- matrix_scored_motifs[element_to_repeat[!is.na(element_to_repeat)],] # Some contigs are missing after filter(nb_occurrence>=5)
	}

	tsne_matrix <- as.matrix(subset(matrix_scored_motifs, select=-c(contig)))
	print_message("    Dimentionality reduction")
	set.seed(tsne_seed)
	tsne_data <- as.data.frame(Rtsne(tsne_matrix, check_duplicates=FALSE, perplexity=tsne_perplexity, max_iter=tsne_max_iter)$Y)

	tsne_data <- tsne_data %>%
		mutate(contig=matrix_scored_motifs$contig) %>% # Link back contig names
		plyr::rename(c('V1'='tSNE_1','V2'='tSNE_2')) %>% # Change colnames
		mutate(contig_length=contigs_info$length[match(contig, str_split(contigs_info$contig, " ", simplify=TRUE)[,1])]) %>%
		mutate(id=paste0(contigs_info$id[match(contig, str_split(contigs_info$contig, " ", simplify=TRUE)[,1])]))

	# Summarize weighted contigs for t-SNE
	if(!is.na(contig_weight_unit) & contig_weight_unit!=-1){ # Don't do with splitted -1 contigs 
		if(contig_weight_unit!=-2){
			print_message("    Finding weigthed contigs centroids")
		}else if(contig_weight_unit==-2){
			print_message("    Finding splitted contigs centroids")
			contigs_info_original <- create.contigs.info(path_metagenome, genome_id, data_type)
			unsplitted_tsne_data <- subset(tsne_data, contig %in% contigs_info_original$contig)
			splitted_tsne_data <- subset(tsne_data, !contig %in% contigs_info_original$contig)

			contigs_info_split <- colsplit(gsub("(.*)_([0-9]+):([0-9]+)","\\1 \\2 \\3",splitted_tsne_data$contig), " ", names=c("genome","start","end"))
			splitted_tsne_data <- cbind(splitted_tsne_data, contigs_info_split)
			splitted_tsne_data$contig <- as.character(splitted_tsne_data$contig) # Generate NA if unexpected contig
			stifle <- foreach(row_idx=seq(1,nrow(splitted_tsne_data))) %do% {
				matching_contigs_info_original <- subset(contigs_info_original, contig==as.character(splitted_tsne_data$genome[row_idx]))
				splitted_tsne_data$contig[row_idx] <- matching_contigs_info_original$contig # Is already in genome
				splitted_tsne_data$contig_length[row_idx] <- matching_contigs_info_original$length

				return(NA)
			}
			splitted_tsne_data$contig <- as.factor(splitted_tsne_data$contig)
			tsne_data <- rbind(subset(splitted_tsne_data, select=-c(genome, start, end)), unsplitted_tsne_data)
		}
		tsne_data <- tsne_data %>%
			group_by(contig) %>%
			mutate(centroid_V1=mean(tSNE_1),centroid_V2=mean(tSNE_2)) %>%
			distinct(contig, .keep_all=TRUE) %>%
			mutate(tSNE_1=centroid_V1,tSNE_2=centroid_V2) %>%
			dplyr::select(-c(centroid_V1, centroid_V2))
	}

	return(tsne_data)
}

big_steps_log10_trans <- function(x){
	trans_new("big_steps_log10", function(x) 1/x, function(x) 1/x)
}

plot.tsne.motifs.score.dev <- function(tsne_data, base_name, is_weighted=FALSE, list_MGEs=NA, list_color=NA){
	if(any(is.na(list_color))){
		contig_groups <- sort(unique(tsne_data$id))
		color_palette <- as.vector(iwanthue(length(contig_groups)))
		names(color_palette) <- contig_groups
	}else{
		contig_groups <- sort(unique(tsne_data$id))
		color_palette <- list_color
		names(color_palette) <- contig_groups
	}

	if(any(is.na(list_MGEs))){
		gp <- ggplot(tsne_data) +
			geom_point(aes(x=tSNE_1 ,y=tSNE_2, col=id, size=contig_length), shape=1) +
			scale_colour_manual(values=color_palette) +
			scale_size_continuous(breaks=c(10000,100000,1000000), labels=c(10000,100000,1000000), range=c(0.1,9), limits=c(min(10000,tsne_data$contig_length),max(1000000,tsne_data$contig_length))) + #, limits=c(min(tsne_data$contig_length),max(tsne_data$contig_length))
			labs(x="t-SNE 1", y="t-SNE 2", colour="Genome of Origin", size="Contigs Length") +
			guides(color=guide_legend(order=1, override.aes=list(shape=15, size=3)), size=guide_legend(order=0)) + # Decreasing order
			theme_bw() +
			theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
			theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
	}else{
		gp <- ggplot() +
			geom_point(data=subset(tsne_data,! contig %in% list_MGEs), aes(x=tSNE_1, y=tSNE_2, col=id, size=contig_length, shape="Genome")) +
			# geom_point(data=subset(tsne_data, contig %in% list_MGEs), aes(x=tSNE_1, y=tSNE_2, col=id, shape="MGEs"), size=3) +
			geom_point(data=subset(tsne_data, contig %in% list_MGEs), aes(x=tSNE_1, y=tSNE_2, shape="MGEs"), size=3) +
			scale_colour_manual(values=color_palette) +
			scale_size_continuous(breaks=c(10000,100000,1000000), labels=c(10000,100000,1000000), range=c(0.1,9), limits=c(min(10000,tsne_data$contig_length),max(1000000,tsne_data$contig_length))) + #, limits=c(min(tsne_data$contig_length),max(tsne_data$contig_length))
			scale_shape_manual(name="Type", values=c("Genome"=1,"MGEs"=2), labels=c("Genome","MGEs")) +
			labs(x="t-SNE 1", y="t-SNE 2", colour="Bin of Origin", size="Contigs Length") +
			guides(color=guide_legend(order=1, override.aes=list(shape=15, size=3))) +
			guides(size=guide_legend(order=0, override.aes=list(shape=1))) + # Decreasing order
			guides(shape=guide_legend(override.aes=list(size=2))) +
			theme_bw() +
			theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
			theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
	}

	# Look for large gap in clustering
	gaps <- tsne_data %>%
		arrange(tSNE_1) %>%
		mutate(diff1=tSNE_1 - lag(tSNE_1, default=0)) %>%
		arrange(tSNE_2) %>%
		mutate(diff2=tSNE_2 - lag(tSNE_2, default=0)) %>%
		filter(diff1>40 | diff2>40)

	if(nrow(gaps)>0){
		gp <- gp +
			theme(legend.position="none")

		pdf(paste0("Contigs_methylation_tsne_",base_name,"_noZoom_v2.pdf"), width=5, height=5)
		print(gp)
		dev.off()

		gp <- gp +
			theme(legend.position="right")

		if(gaps$diff1>40 & gaps$tSNE_1<0){
			gp <- gp +
				coord_cartesian(xlim=c(gaps$tSNE_1, max(tsne_data$tSNE_1)))
		}else if(gaps$diff1>40 & gaps$tSNE_1>0){
			gp <- gp +
				coord_cartesian(xlim=c(min(tsne_data$tSNE_1), gaps$tSNE_1 - gaps$diff1))
		}
		if(gaps$diff2>40 & gaps$tSNE_2<0){
			gp <- gp +
				coord_cartesian(xlim=c(gaps$tSNE_2, max(tsne_data$tSNE_2)))
		}else if(gaps$diff2>40 & gaps$tSNE_2>0){
			gp <- gp +
				coord_cartesian(xlim=c(min(tsne_data$tSNE_2), gaps$tSNE_2 - gaps$diff2))
		}

		if(length(unique(tsne_data$id))>20){ # If too many colors plot separate legend
			gp_with_legend <- gp
			pdf(file=NULL)
			gp_legend <- gglegend(gp_with_legend)
			dev.off()
			gp <- gp +
				theme(legend.position="none")

			pdf(paste0("Contigs_methylation_tsne_",base_name,"_Zoom_v2.pdf"), width=5, height=5)
			print(gp)
			grid.newpage()
			grid.draw(gp_legend)
			dev.off()
		}else{
			pdf(paste0("Contigs_methylation_tsne_",base_name,"_Zoom_v2.pdf"), width=7, height=5)
			print(gp)
			dev.off()
		}
	}else{
		if(length(unique(tsne_data$id))>20){ # If too many colors plot separate legend
			gp_with_legend <- gp
			pdf(file=NULL)
			gp_legend <- gglegend(gp_with_legend)
			dev.off()
			gp <- gp +
				theme(legend.position="none")

			pdf(paste0("Contigs_methylation_tsne_",base_name,"_v2.pdf"), width=5, height=5)
			print(gp)
			grid.newpage()
			grid.draw(gp_legend)
			dev.off()
		}else{
			pdf(paste0("Contigs_methylation_tsne_",base_name,"_v2.pdf"), width=7, height=5)
			print(gp)
			dev.off()
		}
	}

	if(is_weighted){
		gp <- ggplot(tsne_data) +
			geom_point(aes(x=tSNE_1, y=tSNE_2, col=contig_length)) +
			facet_wrap(~id) +
			labs(x="t-SNE 1", y="t-SNE 2", colour="Contig of Origin") +
			theme(legend.position="none")
		pdf(paste0("Contigs_methylation_tsne_",base_name,"_detailed_v2.pdf"), height=9, width=9)
		print(gp)
		dev.off()
	}
}

plot.tsne.motifs.score.custom <- function(tsne_data, base_name, list_MGEs=NA, list_color=NA, my_xlim=NA, my_ylim=NA){
	if(any(is.na(list_color))){
		contig_groups <- sort(unique(tsne_data$id))
		color_palette <- as.vector(iwanthue(length(contig_groups)))
		names(color_palette) <- contig_groups
	}else{
		contig_groups <- sort(unique(tsne_data$id))
		color_palette <- list_color
		names(color_palette) <- contig_groups
	}

	if(any(is.na(list_MGEs))){
		gp <- ggplot(tsne_data) +
			geom_point(aes(x=tSNE_1 ,y=tSNE_2, col=id, size=contig_length), shape=1) +
			scale_colour_manual(values=color_palette) +
			scale_size_continuous(breaks=c(10000,100000,1000000), labels=c(10000,100000,1000000), range=c(0.1,9), limits=c(min(10000,tsne_data$contig_length),max(1000000,tsne_data$contig_length))) + #, limits=c(min(tsne_data$contig_length),max(tsne_data$contig_length))
			labs(x="t-SNE 1", y="t-SNE 2", colour="Genome of Origin", size="Contigs Length") +
			guides(color=guide_legend(order=1, override.aes=list(shape=15, size=3)), size=guide_legend(order=0)) + # Decreasing order
			theme_bw() +
			theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
			theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
	}else{
		gp <- ggplot() +
			geom_point(data=subset(tsne_data,! contig %in% list_MGEs), aes(x=tSNE_1, y=tSNE_2, col=id, size=contig_length, shape="Genome")) +
			# geom_point(data=subset(tsne_data, contig %in% list_MGEs), aes(x=tSNE_1, y=tSNE_2, col=id, shape="MGEs"), size=3) +
			geom_point(data=subset(tsne_data, contig %in% list_MGEs), aes(x=tSNE_1, y=tSNE_2, shape="MGEs"), size=3) +
			scale_colour_manual(values=color_palette) +
			scale_size_continuous(breaks=c(10000,100000,1000000), labels=c(10000,100000,1000000), range=c(0.1,9), limits=c(min(10000,tsne_data$contig_length),max(1000000,tsne_data$contig_length))) + #, limits=c(min(tsne_data$contig_length),max(tsne_data$contig_length))
			scale_shape_manual(name="Type", values=c("Genome"=1,"MGEs"=2), labels=c("Genome","MGEs")) +
			labs(x="t-SNE 1", y="t-SNE 2", colour="Bin of Origin", size="Contigs Length") +
			guides(color=guide_legend(order=1, override.aes=list(shape=15, size=3))) +
			guides(size=guide_legend(order=0, override.aes=list(shape=1))) + # Decreasing order
			guides(shape=guide_legend(override.aes=list(size=2))) +
			theme_bw() +
			theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
			theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
	}

	if(!any(is.na(my_xlim)) & any(is.na(my_ylim))){
		gp <- gp +
			coord_cartesian(xlim=my_xlim)
	}else if(any(is.na(my_xlim)) & !any(is.na(my_ylim))){
		gp <- gp +
			coord_cartesian(ylim=my_ylim)
	}else if(!any(is.na(my_xlim)) & !any(is.na(my_ylim))){
		gp <- gp +
			coord_cartesian(xlim=my_xlim, ylim=my_ylim)
	}

	if(length(unique(tsne_data$id))>20){ # If too many colors plot separate legend
		gp_with_legend <- gp
		pdf(file=NULL)
		gp_legend <- gglegend(gp_with_legend)
		dev.off()
		gp <- gp +
			theme(legend.position="none")

		pdf(paste0("Contigs_methylation_tsne_",base_name,"_v2.pdf"), width=5, height=5)
		print(gp)
		grid.newpage()
		grid.draw(gp_legend)
		dev.off()
	}else{
		pdf(paste0("Contigs_methylation_tsne_",base_name,"_v2.pdf"), width=7, height=5)
		print(gp)
		dev.off()
	}
}

find.tsne.clusters <- function(tsne_binning_de_novo){
	db <- dbscan(subset(tsne_binning_de_novo, select=c(tSNE_1,tSNE_2)), eps=5, minPts=3)
	tsne_binning_de_novo$db_clust <- db$cluster
	tsne_binning_de_novo$id <- paste0("dbscan Bin ", db$cluster)

	return(tsne_binning_de_novo)
}

listing.binned.contigs <- function(tsne_data_dbclust, background_bins){
	list_contigs_tmp <- subset(tsne_data_dbclust, !db_clust %in% background_bins) %>%
		dplyr::select(contig,db_clust) %>%
		group_by(db_clust) %>%
		summarize(contigs=list(as.character(contig)), .groups="drop_last") %>%
		as.list()
	list_contigs <- list_contigs_tmp$contigs
	names(list_contigs) <- list_contigs_tmp$db_clust
	
	return(list_contigs)
}

draw.heatmap.clustering <- function(scored_motifs, list_contigs_subset, base_name, score_threshold, occ_threshold, noRow=FALSE, noCol=FALSE, plot_width=20, plot_height=8){
	# Default color palette
	color_palette <- colorRampPalette(colors=c("darkred", "white", "darkgreen"))

	# Create data.frame with bin annotation 
	bin_annotation <- as.data.frame(as.vector(unlist(list_contigs_subset)), stringsAsFactors=FALSE)
	colnames(bin_annotation) <- "contig"
	bin_annotation$bin <- rep(names(list_contigs_subset),sapply(list_contigs_subset, length))

	# Filter features and convert to matrix-like structure
	scored_motifs <- scored_motifs %>%
		filter(contig %in% unlist(list_contigs_subset)) %>%
		filter(nb_occurrence>=occ_threshold) %>%
		mutate(features=paste0(motif,"_",distance_motif)) %>%
		group_by(features) %>%
		filter(max(abs(dist_score))>=score_threshold) %>%
		ungroup() %>%
		dplyr::select(c(contig, features, dist_score)) %>%
		spread(features, dist_score, fill=0)

	# Annotated bins with names
	scored_motifs_annotated <- merge(scored_motifs, bin_annotation, all.x=TRUE) %>%
		arrange(bin) %>%
		mutate(contig=paste0(contig," (",bin,")"))

	# Convert to matrix
	scored_motifs_annotated_matrix <- as.matrix(subset(scored_motifs_annotated, select=c(-contig, -bin)))
	rownames(scored_motifs_annotated_matrix) <- scored_motifs_annotated$contig
	data_hc_rows <- hclust(dist(scored_motifs_annotated_matrix))
	data_hc_cols <- hclust(dist(t(scored_motifs_annotated_matrix)))

	data_hc_rows_groups <- gsub(".*\\((.*)\\)", "\\1", data_hc_rows$labels)
	# data_hc_cols_color <- gsub(".*\\((.*)\\)", "\1", data_hc_cols$labels)
	color_groups <- as.vector(iwanthue(length(unique(data_hc_rows_groups))))
	names(color_groups) <- sort(unique(data_hc_rows_groups))
	data_hc_rows_color <- as.vector(color_groups[match(data_hc_rows_groups,names(color_groups))])

	hc_features_height_top <- 1
	scale_height <- 1
	hc_features_height_bottom <- scale_height
	hc_features_height <- hc_features_height_top + hc_features_height_bottom
	hc_contigs_width <- 2
	hc_groups_width <- 0.5

	hc_contigs_height <- plot_height - hc_features_height
	hc_features_width <- plot_width - hc_contigs_width - hc_groups_width

	pdf(paste0("Heatmap_clustering_",base_name,"_v1.pdf"), width=plot_width, height=plot_height)
	if(noRow & noCol){
		heatmap.2(
			scored_motifs_annotated_matrix, dendrogram="none", Rowv=FALSE, Colv=FALSE,
			RowSideColors=data_hc_rows_color, #ColSideColors=data_hc_cols_color, 
			scale="none", trace="none",
			col=rev(color_palette(100)), breaks=seq(-5,5, length=101),
			keysize=0.5, key.title=NA, key.xlab="Methylation feature values", density.info="none",
			lmat=rbind(c(0,0,4),c(5,0,4),c(3,1,2)), lwid=c(hc_contigs_width,hc_groups_width,hc_features_width), lhei=c(hc_features_height_top,scale_height,hc_contigs_height), margins=c(7,8.5)
		)
	}else if(noRow){
		heatmap.2(
			scored_motifs_annotated_matrix, dendrogram="column", Rowv=FALSE, Colv=as.dendrogram(data_hc_cols),
			RowSideColors=data_hc_rows_color, #ColSideColors=data_hc_cols_color, 
			scale="none", trace="none",
			col=rev(color_palette(100)), breaks=seq(-5,5, length=101),
			keysize=0.5, key.title=NA, key.xlab="Methylation feature values", density.info="none",
			lmat=rbind(c(0,0,4),c(5,0,4),c(3,1,2)), lwid=c(hc_contigs_width,hc_groups_width,hc_features_width), lhei=c(hc_features_height_top,scale_height,hc_contigs_height), margins=c(7,8.5)
		)
	}else if(noCol){
		heatmap.2(
			scored_motifs_annotated_matrix, dendrogram="row", Rowv=as.dendrogram(data_hc_rows), Colv=FALSE,
			RowSideColors=data_hc_rows_color, #ColSideColors=data_hc_cols_color, 
			scale="none", trace="none",
			col=rev(color_palette(100)), breaks=seq(-5,5, length=101),
			keysize=0.5, key.title=NA, key.xlab="Methylation feature values", density.info="none",
			lmat=rbind(c(0,0,4),c(5,0,4),c(3,1,2)), lwid=c(hc_contigs_width,hc_groups_width,hc_features_width), lhei=c(hc_features_height_top,scale_height,hc_contigs_height), margins=c(7,8.5)
		)
	}else{
		heatmap.2(
			scored_motifs_annotated_matrix, Rowv=as.dendrogram(data_hc_rows), Colv=as.dendrogram(data_hc_cols),
			RowSideColors=data_hc_rows_color, #ColSideColors=data_hc_cols_color, 
			scale="none", trace="none",
			col=rev(color_palette(100)), breaks=seq(-5,5, length=101),
			keysize=0.5, key.title=NA, key.xlab="Methylation feature values", density.info="none",
			lmat=rbind(c(0,0,4),c(5,0,4),c(3,1,2)), lwid=c(hc_contigs_width,hc_groups_width,hc_features_width), lhei=c(hc_features_height_top,scale_height,hc_contigs_height), margins=c(7,8.5)
		)
	}
	par(lend=1)
	legend("topright", legend=names(color_groups), col=color_groups, lty=1, lwd=10)
	dev.off()
}

plot.tsne.other.color.scheme <- function(main_tsne, color_tsne, plot_name, my_xlim=NA, my_ylim=NA){
	main_tsne_colored <- merge(subset(main_tsne, select=c(-id)), subset(color_tsne, select=c(contig, id)), by=c("contig"), all.x=TRUE)
	main_tsne_colored$id <- ifelse(is.na(main_tsne_colored$id),"none",main_tsne_colored$id) # Add color for contigs missing from annotation

	if(any(is.na(my_xlim)) | any(is.na(my_ylim))){
		plot.tsne.motifs.score(main_tsne_colored, plot_name)
	}else{
		plot.tsne.motifs.score.custom(main_tsne_colored, plot_name, NA, NA, my_xlim, my_ylim)
	}
}

generate.color.palette <- function(annotation){
	contig_groups <- levels(annotation$id)
	if("Not binned" %in% contig_groups){
		color_palette <-  c(as.vector(iwanthue(length(contig_groups) - 1)), "#000000")
		names(color_palette) <- c(contig_groups[!contig_groups %in% "Not binned"], "Not binned")
	}else if("Unknown" %in% contig_groups){
		color_palette <-  c(as.vector(iwanthue(length(contig_groups) - 1)), "#000000")
		names(color_palette) <- c(contig_groups[!contig_groups %in% "Unknown"], "Unknown")
	}else if("No annotation" %in% contig_groups){
		color_palette <-  c("#000000")
		names(color_palette) <- c("No annotation")
	}else{
		color_palette <- c(as.vector(iwanthue(length(contig_groups))), "#000000")
		names(color_palette) <- c(contig_groups, "Not binned")
	}

	return(color_palette)
}

add.contigs.annotation <- function(tsne_data, annotation){
	tsne_data$id <- NULL
	tsne_data$db_clust <- NULL
	annotation$contig <- str_split(annotation$contig," ",simplify=TRUE)[,1]

	tsne_data <- merge(tsne_data, annotation, by=c("contig"), all.x=TRUE)
	if("Not binned" %in% levels(tsne_data$id)){
		tsne_data$id <- factor(tsne_data$id, levels=levels(tsne_data$id))
	}else if("No annotation" %in% levels(tsne_data$id)){
		tsne_data$id <- factor(tsne_data$id, levels=levels(tsne_data$id))
	}else if("Unknown" %in% levels(tsne_data$id)){
		tsne_data$id <- factor(tsne_data$id, levels=levels(tsne_data$id))
	}else{
		tsne_data$id <- factor(tsne_data$id, levels=c(levels(tsne_data$id), "Not binned"))
		tsne_data$id[is.na(tsne_data$id)] <- as.factor("Not binned")
	}

	return(tsne_data)
}

plot.tsne.motifs.score <- function(tsne_data, annotation, base_name, color_palette, list_MGEs=NA, my_xlim=NA, my_ylim=NA, min_len=25000, path_output="./"){
	size_default <- 7
	size_legend <- 6
	default_text <- theme(
		text=element_text(size=size_default),
		axis.title=element_text(size=size_default, face="bold"),
		axis.text=element_text(size=size_default),
		strip.text=element_text(size=size_default),
		legend.title=element_text(size=size_default),
		legend.text=element_text(size=size_legend),
		legend.margin=ggplot2::margin(0,0,0,0, unit="in"),
		legend.box.margin=ggplot2::margin(0,0,0,0, unit="in"),
		legend.box.spacing=ggplot2::margin(0.02,0.02,0.02,0.02, unit="in"),
		legend.key.width=unit(0.1, "in"),
		legend.key.height=unit(0.1, "in")
	)

	tsne_data <- add.contigs.annotation(tsne_data, annotation)

	tsne_data <- tsne_data %>%
		mutate(contig_length=ifelse(contig_length>1000000,1000001,contig_length))

	if(any(is.na(list_MGEs))){
		gp <- ggplot(tsne_data) +
			geom_point(aes(x=tSNE_1 ,y=tSNE_2, col=id, size=contig_length), shape=1) +
			scale_colour_manual(values=color_palette) +
			scale_size_continuous(breaks=c(min_len,100000,1000001), labels=c(min_len,100000,"> 1000000"), range=c(0.1,9), limits=c(min_len,1000001)) +
			labs(x="t-SNE 1", y="t-SNE 2", colour="Genome of Origin", size="Contig length") +
			guides(color=guide_legend(order=1, override.aes=list(shape=15, size=3)), size=guide_legend(order=2)) +
			theme_bw() +
			theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
			theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			theme(legend.text.align=0) +
			default_text +
			labs(title=paste0("Methylation binning of ",base_name), color=NULL) + # "Bin of origin"
			guides(size=guide_legend(keyheight=0.20, default.unit="in", override.aes=list(shape=1)))
	}else{
		gp <- ggplot() +
			geom_point(data=subset(tsne_data,! contig %in% list_MGEs), aes(x=tSNE_1, y=tSNE_2, col=id, size=contig_length, shape="Genome")) +
			geom_point(data=subset(tsne_data, contig %in% list_MGEs), aes(x=tSNE_1, y=tSNE_2, shape="MGEs"), size=3) +
			scale_colour_manual(values=color_palette) +
			scale_size_continuous(breaks=c(min_len,100000,1000001), labels=c(min_len,100000,"> 1000000"), range=c(0.1,9), limits=c(min_len,1000001)) +
			# scale_size_continuous(breaks=c(min_len,100000,1000000), labels=c(min_len,100000,1000000), range=c(0.1,9), limits=c(min(min_len,tsne_data$contig_length),max(1000000,tsne_data$contig_length))) +
			scale_shape_manual(name="Contig type", values=c("Genome"=1,"MGEs"=2), labels=c("Genome","MGEs")) +
			labs(x="t-SNE 1", y="t-SNE 2", colour="Bin of Origin", size="Contig length") +
			guides(color=guide_legend(order=1, override.aes=list(shape=15, size=3))) +
			guides(size=guide_legend(order=2, override.aes=list(shape=1))) +
			guides(shape=guide_legend(order=3, override.aes=list(size=2))) +
			theme_bw() +
			theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
			theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			theme(legend.text.align=0) +
			default_text +
			labs(title=paste0("Methylation binning of ",base_name), color=NULL) + # "Bin of origin"
			guides(size=guide_legend(keyheight=0.20, default.unit="in", override.aes=list(shape=1)))
	}

	if(!any(is.na(my_xlim)) & any(is.na(my_ylim))){
		gp <- gp +
			coord_cartesian(xlim=my_xlim)
	}else if(any(is.na(my_xlim)) & !any(is.na(my_ylim))){
		gp <- gp +
			coord_cartesian(ylim=my_ylim)
	}else if(!any(is.na(my_xlim)) & !any(is.na(my_ylim))){
		gp <- gp +
			coord_cartesian(xlim=my_xlim, ylim=my_ylim)
	}

	if(length(unique(tsne_data$id))>20){ # If too many colors plot separate legend
		gp_with_legend <- gp
		pdf(file=NULL)
		gp_legend <- gglegend(gp_with_legend)
		dev.off()
		gp <- gp +
			theme(legend.position="none")

		pdf(paste0(path_output,"Contigs_methylation_tsne_",base_name,".pdf"), width=4.25, height=3.41)
		print(gp)
		grid.newpage()
		grid.draw(gp_legend)
		dev.off()
	}else{
		pdf(paste0(path_output,"Contigs_methylation_tsne_",base_name,".pdf"), width=4.25, height=3.41)
		print(gp)
		dev.off()
	}

	if(length(unique(tsne_data$id))>20){
		
		return(list(main=gp, legend=gp_legend))
	}else{

		return(list(main=gp))
	}
}

#  ______ _ _ _             ___  ___      _   _  __     
#  |  ___(_) | |            |  \/  |     | | (_)/ _|    
#  | |_   _| | |_ ___ _ __  | .  . | ___ | |_ _| |_ ___ 
#  |  _| | | | __/ _ \ '__| | |\/| |/ _ \| __| |  _/ __|
#  | |   | | | ||  __/ |    | |  | | (_) | |_| | | \__ \
#  \_|   |_|_|\__\___|_|    \_|  |_/\___/ \__|_|_| |___/
#                                                       
#                                                       

explore.filtering.metrics <- function(scored_motifs, contigs_info, base_name, signal_ratio_threshold, length_threshold, occ_threshold, score_threshold, figure_type){
	if("dist_score" %in% colnames(scored_motifs)){ # With distance_motif
		filtered_annotated_scored_motifs <- annotate.score(scored_motifs, contigs_info) %>%
			filter(length>=length_threshold) %>%
			mutate(to_keep=ifelse(signal_ratio>signal_ratio_threshold,"Above ratio threshold","Below ratio threshold")) %>%
			mutate(is_kept=ifelse(nb_occurrence>=occ_threshold & (dist_score<=-score_threshold | dist_score>=score_threshold),16,46))

		gp <- ggplot(filtered_annotated_scored_motifs) +
			geom_point(aes(x=log10(nb_occurrence), y=dist_score, col=signal_ratio, pch=is_kept)) +
			scale_shape_identity() +
			geom_vline(xintercept=log10(occ_threshold)) +
			geom_hline(yintercept=c(-score_threshold, score_threshold)) +
			facet_wrap(~to_keep, ncol=1) +
			labs(x="Number of motif occurrences in contig (log10)", y="Score relative position to motif", colour="Ratio of base\nwith expected signal")

		gp_detail <- ggplot(filtered_annotated_scored_motifs) +
			geom_point(aes(x=dist_score, y=log10(nb_occurrence), col=signal_ratio, pch=is_kept)) +
			scale_shape_identity() +
			geom_hline(yintercept=log10(occ_threshold)) +
			geom_vline(xintercept=c(-score_threshold, score_threshold)) +
			facet_wrap(~id) +
			labs(x="Score relative position to motif", y="Number of motif occurrences in contig (log10)", colour="Ratio of base\nwith expected signal")
	}else{ # Without distances
		filtered_annotated_scored_motifs <- annotate.score(scored_motifs, contigs_info) %>%
			filter(length>=length_threshold) %>%
			mutate(to_keep=ifelse(average_signal_ratio>signal_ratio_threshold,"Above ratio threshold","Below ratio threshold")) %>%
			mutate(is_kept=ifelse(average_nb_occurrence>=occ_threshold & score>=score_threshold,16,46))

		gp <- ggplot(filtered_annotated_scored_motifs) +
			geom_point(aes(x=log10(average_nb_occurrence), y=score, col=average_signal_ratio, pch=is_kept)) +
			scale_shape_identity() +
			geom_vline(xintercept=log10(occ_threshold)) +
			geom_hline(yintercept=score_threshold) +
			facet_wrap(~to_keep, ncol=1) +
			labs(x="Number of motif occurrences in contig (log10)", y="Motif score", colour="Ratio of base\nwith expected signal")

		gp_detail <- ggplot(filtered_annotated_scored_motifs) +
			geom_point(aes(x=score, y=log10(average_nb_occurrence), col=average_signal_ratio, pch=is_kept)) +
			scale_shape_identity() +
			geom_hline(yintercept=log10(occ_threshold)) +
			geom_vline(xintercept=score_threshold) +
			facet_wrap(~id) +
			labs(x="Motif score", y="Number of motif occurrences in contig (log10)", colour="Ratio of base\nwith expected signal")
	}
	if(figure_type=="pdf"){
		pdf(paste0("Filtering_features.",base_name,".mL",length_threshold,"_mR",signal_ratio_threshold,"_mO",occ_threshold,"_mS",score_threshold,"_v1.pdf"), width=10, height=5)
		print(gp)
		dev.off()

		pdf(paste0("Filtering_features_per_genome.",base_name,".mL",length_threshold,"_mR",signal_ratio_threshold,"_mO",occ_threshold,"_mS",score_threshold,"_v1.pdf"), width=15, height=15)
		print(gp_detail)
		dev.off()
	}else if(figure_type=="png"){
		png(paste0("Filtering_features.",base_name,".mL",length_threshold,"_mR",signal_ratio_threshold,"_mO",occ_threshold,"_mS",score_threshold,"_v1.png"), width=4096, height=2048)
		print(gp)
		dev.off()

		png(paste0("Filtering_features_per_genome.",base_name,".mL",length_threshold,"_mR",signal_ratio_threshold,"_mO",occ_threshold,"_mS",score_threshold,"_v1.png"), width=6144, height=6144)
		print(gp_detail)
		dev.off()
	}
}

select.features <- function(scored_motifs, contigs_info, length_threshold, occ_threshold, score_threshold){
	contigs_info$contig <- str_split(contigs_info$contig," ",simplify=TRUE)[,1]
	annotated_scored_motifs <- annotate.score(scored_motifs, contigs_info) # Add length but could be avoided
	rm(scored_motifs)
	gc()
	
	if("dist_score" %in% colnames(annotated_scored_motifs)){ # With distance_motif
		selected_features <- annotated_scored_motifs %>%
			filter(length >= length_threshold) %>%
			filter(nb_occurrence >= occ_threshold) %>%
			filter(dist_score <= -score_threshold | dist_score >= score_threshold) %>%
			mutate(feature_name=paste0(motif,"_",distance_motif)) %>% # Integrate relative position in motif name
			group_by(feature_name) %>%
			mutate(contigs_origin=paste(contig, collapse="|")) %>%
			dplyr::select(feature_name, motif, contigs_origin) %>%
			distinct(feature_name, .keep_all=TRUE)
	}else{ # Without distance_motif
		selected_features <- annotated_scored_motifs %>%
			filter(average_nb_occurrence>=occ_threshold) %>%
			filter(length>=length_threshold) %>%
			filter(score>=score_threshold) %>%
			dplyr::select(motif)
	}

	return(selected_features)
}

filter.features <- function(scored_motifs, selected_features){
	if("dist_score" %in% colnames(scored_motifs)){ # With distance_motif
		filtered_scored_motifs <- scored_motifs %>%
			mutate(feature_name=paste0(motif,"_",distance_motif)) %>%
			filter(feature_name %in% selected_features$feature_name) %>%
			dplyr::select(-feature_name)
	}else{ # Without distance_motif
		filtered_scored_motifs <- scored_motifs %>%
			mutate(feature_name=motif) %>%
			filter(feature_name %in% selected_features$motif) %>%
			dplyr::select(-feature_name)
	}

	return(filtered_scored_motifs)
}

select.contigs.motifs.filtering <- function(path_metagenome, min_length, name_output){
	sequence_metagenome <- readDNAStringSet(path_metagenome)
	motifs_filtering_sequence_metagenome <- sequence_metagenome[width(sequence_metagenome)>=min_length]
	writeXStringSet(motifs_filtering_sequence_metagenome, name_output)
}

score.contigs.alignments <- function(path_metagenome, bamFile, genome_bin_id, base_name, data_type, nbCPU){
	whats <- c("qname","qwidth") # c("qname","mapq","qwidth","cigar")
	tags <- c("AS") # c("NM", "AS")
	sequence_metagenome <- readDNAStringSet(path_metagenome)
	list_contigs_name <- str_split(names(sequence_metagenome)," ",simplify=TRUE)[,1]

	registerDoMC(nbCPU)
	alignment_information <- foreach(contig_name=list_contigs_name, .combine=rbind) %do%{
		contig_length <- width(subset(sequence_metagenome, str_split(names(sequence_metagenome)," ",simplify=TRUE)[,1]==contig_name))
		range <- GRanges(seqnames=contig_name, ranges=IRanges(1,contig_length))
		# param <- ScanBamParam(which=range, what=c("qname","mapq","qwidth","cigar"), tag=c("NM", "AS"))
		param <- ScanBamParam(which=range, what=whats, tag=tags)

		bam <- scanBam(bamFile, param=param) # Can fail if not correctly indexed

		alignment_information_contig <- as.data.frame(bam[[1]], stringsAsFactors=FALSE) # TODO fail if no reads mapped to chunk

		if(nrow(alignment_information_contig)==0){
			list_columns <- c(whats, tags)
			alignment_information_contig <- data.frame(matrix(NA, nrow=1, ncol=length(list_columns)))
			colnames(alignment_information_contig) <- list_columns
		}

		alignment_information_contig$contig <- as.factor(contig_name)

		return(alignment_information_contig)
	}
	registerDoSEQ()

	# Add binning information
	if(data_type=="real"){
		alignment_information <- alignment_information %>%
			mutate(bin=genome_bin_id$id[match(contig, str_split(genome_bin_id$contig," ", simplify=TRUE)[,1])])
	}else if(data_type=="mock"){ # NEW
		alignment_information <- alignment_information %>%
			mutate(bin=genome_bin_id$id[match(gsub("(.*)_([0-9]+):([0-9]+)","\\1",contig), str_split(genome_bin_id$contig_name," ", simplify=TRUE)[,1])])
	}

	summary_alignment_information <- alignment_information %>%
		group_by(bin, contig) %>%
		summarize(score=median(AS/qwidth), .groups="drop_last")

	if(nrow(summary_alignment_information)<20){
		gp <- ggplot(alignment_information) +
			geom_density(aes(x=AS/qwidth, col=bin, group=contig)) + # is tag.AS if more than one tag
			geom_vline(data=summary_alignment_information, aes(xintercept=score, col=bin, group=contig)) +
			labs(title="Alignment quality") +
			labs(x="Distribution of reads score (AS/read length)", y="Density")
		pdf(paste0("Alignment_quality_",base_name,"_v1.pdf"), width=10)
		print(gp)
		dev.off()
	}else{
		gp <- ggplot(summary_alignment_information) +
			geom_density(aes(x=score, col=bin)) + # is tag.AS if more than one tag
			facet_wrap(~bin) +
			coord_cartesian(ylim=c(0,8)) +
			labs(title="Alignment quality") +
			labs(x="Distribution of reads score (AS/read length)", y="Density")
		pdf(paste0("Alignment_quality_",base_name,"_v1.pdf"), width=10)
		print(gp)
		dev.off()			
	}
}

filtering.overlapping.motifs <- function(selected_dist_real_features, list_lengths, scope, nbCPU){
	summary_selected_dist_real_features <- selected_dist_real_features %>%
		group_by(motif) %>%
		summarize(n=n(), .groups="drop_last") %>%
		mutate(motif=as.character(motif), len=nchar(motif))

	if(scope=="all"){
		list_motifs_to_filter <- summary_selected_dist_real_features
	}else if(scope=="bipartite"){
		list_motifs_to_filter <- subset(summary_selected_dist_real_features, len>8)
	}

	registerDoMC(nbCPU)
	base_Nmer_motifs <- foreach(len_motifs=list_lengths) %do% {
		summary_overlapping_motifs <- NA
		if(scope=="all"){
			base_kmer_motifs <- subset(list_motifs_to_filter, len==len_motifs)$motif
		}else if(scope=="bipartite"){
			base_kmer_motifs <- subset(summary_selected_dist_real_features, len==len_motifs)$motif
		}
		print(paste0("Processing ",len_motifs,"-mer (n=",length(base_kmer_motifs),") motifs."))

		list_motifs_to_filter <- subset(list_motifs_to_filter,! motif %in% base_kmer_motifs) # Remove current motifs
		while(!is.null(summary_overlapping_motifs)){
			summary_overlapping_motifs <- foreach(base_kmer_motif=base_kmer_motifs, .combine=rbind) %dopar% {
				number_of_matching_motifs <- length(list_motifs_to_filter$motif[grep(base_kmer_motif, list_motifs_to_filter$motif)])

				return(data.frame(motif=base_kmer_motif, n_matches=number_of_matching_motifs))
			}
			summary_overlapping_motifs <- summary_overlapping_motifs %>%
				arrange(desc(n_matches)) %>%
				filter(n_matches>0)
			number_motifs_to_process <- nrow(summary_overlapping_motifs)

			if(number_motifs_to_process==0){
				print(paste0("No matches remaining."))
				summary_overlapping_motifs <- NULL
			}else{
				# TODO look for correlation before removing?
				matches_to_remove <- grepl(summary_overlapping_motifs$motif[1], list_motifs_to_filter$motif)
				print(paste0("Removing matches for ",summary_overlapping_motifs$motif[1]," (n=",sum(matches_to_remove, na.rm=TRUE),"). Up to ",number_motifs_to_process-1," motifs to process."))
				list_motifs_to_filter <- list_motifs_to_filter[!matches_to_remove,]
			}
		}

		return(base_kmer_motifs)
	}
	registerDoSEQ()
	names(base_Nmer_motifs) <- list_lengths

	remaining_motifs <- c(as.vector(unlist(base_Nmer_motifs)), list_motifs_to_filter$motif)

	selected_dist_real_features <- subset(selected_dist_real_features, motif %in% remaining_motifs)

	return(selected_dist_real_features)
}

convert.motif.grep <- function(motif){
	grep_iupac_nc <- data.frame(
		code=c("A","C","G","T","R","Y","S","W","K","M","B","D","H","V","N"),
		grep_pattern=c("A","C","G","T","[AG]","[CT]","[CG]","[AT]","[GT]","[AC]","[CGT]","[AGT]","[ACT]","[ACG]","[ACGTN]"),
		stringsAsFactors=FALSE
	)
	splitted_motif <- str_split(motif, "", simplify=TRUE)
	converted_splitted_motif <- grep_iupac_nc$grep_pattern[match(splitted_motif, grep_iupac_nc$code)]
	grep_motif <- paste(converted_splitted_motif, collapse="")

	return(grep_motif)
}

plot.associated.motifs.scores <- function(selected_contigSubset_dist_real_features, scored_contigSubset_dist_real_motifs, motif_of_interest, base_name){
	selected_feature_annotation <- subset(selected_contigSubset_dist_real_features, grepl(convert.motif.grep(motif_of_interest), motif)) %>%
		mutate(distance_motif=as.numeric(str_split(feature_name, "_", simplify=TRUE)[,2])) %>%
		group_by(feature_name, motif, distance_motif) %>%
		expand(contig=as.vector(str_split(contigs_origin, "\\|", simplify=TRUE)))

	subset_scored_motifs <- subset(scored_contigSubset_dist_real_motifs, grepl(convert.motif.grep(motif_of_interest), motif))
	gp <- ggplot(subset_scored_motifs) +
		geom_tile(aes(x=distance_motif, y=contig, fill=abs(dist_score))) +
		scale_fill_distiller(palette="Spectral") +
		geom_point(data=subset(subset_scored_motifs, abs(dist_score)>1.5), aes(x=distance_motif, y=contig), col=1) +
		geom_point(data=selected_feature_annotation, aes(x=distance_motif, y=contig), col=2) +
		facet_wrap(~motif) +
		labs(x="Motif distance", y="Contigs") +
		labs(fill="Methylation\nscore")
	pdf(paste0("Associated_motifs_scores_",motif_of_interest,"_",base_name,".pdf"), width=15, height=15)
	print(gp)
	dev.off()
}

summarize.motifs.filtering <- function(list_motifs_to_process, selected_dist_real_features, list_expected_motifs){
	# Look for selected motifs matching expected motifs
	df_associated_motifs <- foreach(expected_motif=list_expected_motifs, .combine=rbind) %do% {
		matching_motifs <- list_motifs_to_process[grep(convert.motif.grep(expected_motif), list_motifs_to_process)]
		
		if(length(matching_motifs)>0){
			results <- data.frame(expected_motif=expected_motif, matching_motifs=matching_motifs, stringsAsFactors=FALSE)
		}else{
			results <- data.frame(expected_motif=expected_motif, matching_motifs=NA, stringsAsFactors=FALSE)
		}

		return(results)
	}

	# How many selected motifs matched expected motifs?
	print("Worst expected motifs:")
	df_associated_motifs %>%
		group_by(expected_motif) %>%
		mutate(score=ifelse(is.na(matching_motifs), 0, 1)) %>%
		summarize(n=sum(score), .groups="drop_last") %>%
		arrange(n) %>%
		print(n=20)

	# Which selected motifs seems FP?
	explained_motifs <- as.data.frame(list_motifs_to_process[list_motifs_to_process %in% df_associated_motifs$matching_motifs], stringsAsFactors=FALSE)
	colnames(explained_motifs) <- c("motif")
	summary_explained_motifs <- explained_motifs %>% 
		mutate(len=nchar(motif)) %>%
		arrange(len) %>% group_by(len) %>%
		summarize(n=n(), perc=n/nrow(explained_motifs), .groups="drop_last")
	print("Observed motifs explained by length:")
	summary_explained_motifs  %>% print()

	unexplained_motifs <- as.data.frame(list_motifs_to_process[!list_motifs_to_process %in% df_associated_motifs$matching_motifs], stringsAsFactors=FALSE)
	colnames(unexplained_motifs) <- c("motif")
	summary_unexplained_motifs <- unexplained_motifs %>% 
		mutate(len=nchar(motif)) %>%
		arrange(len) %>%
		group_by(len) %>%
		summarize(n=n(), perc=n/nrow(unexplained_motifs), .groups="drop_last")
	print("Observed motifs unexplained by length:")
	summary_unexplained_motifs %>% print()

	# Which contig produced FP
	explained_features_per_contig <- selected_dist_real_features %>%
		filter(motif %in% list_motifs_to_process[list_motifs_to_process %in% df_associated_motifs$matching_motifs]) %>%
		group_by(contigs_origin) %>%
		summarize(n=n(), .groups="drop_last") %>%
		mutate(perc=n/sum(n)) %>%
		arrange(desc(n))
	print("Explained features per contig:")
	explained_features_per_contig %>% print()
	unexplained_features_per_contig <- selected_dist_real_features %>%
		filter(motif %in% list_motifs_to_process[!list_motifs_to_process %in% df_associated_motifs$matching_motifs]) %>%
		group_by(contigs_origin) %>%
		summarize(n=n(), .groups="drop_last") %>%
		mutate(perc=n/sum(n)) %>%
		arrange(desc(n))
	print("Unexplained features per contig:")
	unexplained_features_per_contig %>% print()

	explained_motifs_per_contig <- selected_dist_real_features %>%
		filter(motif %in% list_motifs_to_process[list_motifs_to_process %in% df_associated_motifs$matching_motifs]) %>%
		group_by(motif, contigs_origin) %>%
		summarize(n_features=n(), .groups="drop_last") %>%
		group_by(contigs_origin) %>%
		summarize(n=n(), .groups="drop_last") %>%
		mutate(perc=n/sum(n)) %>%
		arrange(desc(n))
	print("Explained motifs per contig:")
	explained_motifs_per_contig %>% print()
	unexplained_motifs_per_contig <- selected_dist_real_features %>%
		filter(motif %in% list_motifs_to_process[!list_motifs_to_process %in% df_associated_motifs$matching_motifs]) %>%
		group_by(motif, contigs_origin) %>%
		summarize(n_features=n(), .groups="drop_last") %>%
		group_by(contigs_origin) %>%
		summarize(n=n(), .groups="drop_last") %>%
		mutate(perc=n/sum(n)) %>%
		arrange(desc(n))
	print("Unexplained motifs per contig:")
	unexplained_motifs_per_contig %>% print()

	results <- list(
		associated_motifs=df_associated_motifs,
		explained_motifs=explained_motifs,
		summary_explained_motifs=summary_explained_motifs,
		unexplained_motifs=unexplained_motifs,
		summary_unexplained_motifs=summary_unexplained_motifs,
		explained_features_per_contig=explained_features_per_contig,
		unexplained_features_per_contig=unexplained_features_per_contig,
		explained_motifs_per_contig=explained_motifs_per_contig,
		unexplained_motifs_per_contig=unexplained_motifs_per_contig
	)

	return(results)
}

repetitive.grepl <- function(list_motif_patterns, list_references){
	res <- foreach(motif_pattern=list_motif_patterns, .combine=c) %do% {
		return(any(grepl(motif_pattern, list_references)))
	}

	return(res)
}

ignore.features.from.bin <- function(selected_features, bin_name, genome_bin_id){
	genome_bin_id$contig_short_name <- str_split(genome_bin_id$contig, " ", simplify=TRUE)[,1]
	contigs_to_ignore <- subset(genome_bin_id, id %in% bin_name | contig %in% bin_name)$contig_short_name
	contigs_to_keep <- subset(genome_bin_id, !contig_short_name %in% contigs_to_ignore)$contig_short_name

	list_contigs_origin <- unique(selected_features$contigs_origin)
	list_contigs_origin_to_keep <- list_contigs_origin[repetitive.grepl(list_contigs_origin, contigs_to_keep)]

	noBinX_selected_features <- subset(selected_features, contigs_origin %in% list_contigs_origin_to_keep)

	return(noBinX_selected_features)
}

show.motif.scores.distribution <- function(motif_of_interest, scored_dist_real_motifs, genome_bin_id, score_threshold_range, minOcc){
	subset_scored_dist_real_motifs <- subset(scored_dist_real_motifs, grepl(convert.motif.grep(motif_of_interest),motif))
	subset_scored_dist_real_motifs$bin <- genome_bin_id$id[match(subset_scored_dist_real_motifs$contig, str_split(genome_bin_id$contig, " ", simplify=TRUE)[,1])]

	gp <- ggplot(subset_scored_dist_real_motifs) +
		geom_point(aes(x=nb_occurrence, y=dist_score, col=bin)) +
		geom_hline(yintercept=score_threshold_range) +
		geom_vline(xintercept=minOcc) +
		coord_cartesian(ylim=c(-20,20)) +
		labs(title=paste0(motif_of_interest), subtitle=paste0("Score range from ",score_threshold_range[1]," to ",score_threshold_range[2],"; Minimum Occ. ",minOcc)) +
		labs(x="Number Occurrences features per contigs", y="Features Score", col="Bin")
	print(gp)
}

show.global.scores.distribution <- function(scored_dist_real_motifs, motifs_filtering_summary, genome_bin_id, base_name, score_threshold_range, minOcc){
	subset_scored_dist_real_motifs <- subset(scored_dist_real_motifs, motif %in% motifs_filtering_summary$unexplained_motifs$motif)
	subset_scored_dist_real_motifs$bin <- genome_bin_id$id[match(subset_scored_dist_real_motifs$contig, str_split(genome_bin_id$contig, " ", simplify=TRUE)[,1])]

	gp_short <- ggplot(subset(subset_scored_dist_real_motifs, nchar(as.character(motif))<8)) +
		geom_point(aes(x=log10(nb_occurrence), y=dist_score, col=bin)) +
		geom_hline(yintercept=score_threshold_range) +
		geom_vline(xintercept=log10(minOcc)) +
		facet_wrap(~bin) +
		coord_cartesian(ylim=c(-10,10)) +
		labs(title="Unexplained short", subtitle=paste0("Score range from ",score_threshold_range[1]," to ",score_threshold_range[2],"; Minimum Occ. ",paste0(minOcc, collapse=","))) +
		labs(x="Number Occurrences features per contigs", y="Features Score", col="Bin")
	gp_long <- ggplot(subset(subset_scored_dist_real_motifs, nchar(as.character(motif))>8)) +
		geom_point(aes(x=log10(nb_occurrence), y=dist_score, col=bin)) +
		geom_hline(yintercept=score_threshold_range) +
		geom_vline(xintercept=log10(minOcc)) +
		facet_wrap(~bin) +
		coord_cartesian(ylim=c(-10,10)) +
		labs(title="Unexplained long", subtitle=paste0("Score range from ",score_threshold_range[1]," to ",score_threshold_range[2],"; Minimum Occ. ",paste0(minOcc, collapse=","))) +
		labs(x="Number Occurrences features per contigs", y="Features Score", col="Bin")
	pdf(paste0("Motifs_score_distribution_unexplained_",base_name,"_v1.pdf"), width=15, height=15)
	print(gp_short)
	print(gp_long)
	dev.off()

	subset_scored_dist_real_motifs <- subset(scored_dist_real_motifs, motif %in% motifs_filtering_summary$explained_motifs$motif)
	subset_scored_dist_real_motifs$bin <- genome_bin_id$id[match(subset_scored_dist_real_motifs$contig, str_split(genome_bin_id$contig, " ", simplify=TRUE)[,1])]

	gp_short <- ggplot(subset(subset_scored_dist_real_motifs, nchar(as.character(motif))<8)) +
		geom_point(aes(x=log10(nb_occurrence), y=dist_score, col=bin)) +
		geom_hline(yintercept=score_threshold_range) +
		geom_vline(xintercept=log10(minOcc)) +
		facet_wrap(~bin) +
		coord_cartesian(ylim=c(-10,10)) +
		labs(title="Explained short", subtitle=paste0("Score range from ",score_threshold_range[1]," to ",score_threshold_range[2],"; Minimum Occ. ",paste0(minOcc, collapse=","))) +
		labs(x="Number Occurrences features per contigs", y="Features Score", col="Bin")
	gp_long <- ggplot(subset(subset_scored_dist_real_motifs, nchar(as.character(motif))>8)) +
		geom_point(aes(x=log10(nb_occurrence), y=dist_score, col=bin)) +
		geom_hline(yintercept=score_threshold_range) +
		geom_vline(xintercept=log10(minOcc)) +
		facet_wrap(~bin) +
		coord_cartesian(ylim=c(-10,10)) +
		labs(title="Explained long", subtitle=paste0("Score range from ",score_threshold_range[1]," to ",score_threshold_range[2],"; Minimum Occ. ",paste0(minOcc, collapse=","))) +
		labs(x="Number Occurrences features per contigs", y="Features Score", col="Bin")
	pdf(paste0("Motifs_score_distribution_explained_",base_name,"_v1.pdf"), width=15, height=15)
	print(gp_short)
	print(gp_long)
	dev.off()
}

filtering.sig.motifs <- function(selected_dist_real_features, type, values){
	summary_selected_dist_real_features <- selected_dist_real_features %>%
		group_by(motif) %>%
		summarize(n=n(), .groups="drop_last") %>%
		mutate(motif=as.character(motif), len=nchar(motif))
	if(type=="only"){
		summary_selected_dist_real_features <- summary_selected_dist_real_features %>%
			filter(n %in% values)
	}else if(type=="at_least"){
		summary_selected_dist_real_features <- summary_selected_dist_real_features %>%
			filter(n>=values)
	}
	selected_dist_real_features <- subset(selected_dist_real_features, motif %in% c(summary_selected_dist_real_features$motif))

	return(selected_dist_real_features)
}

filtering.sig.bipartite.motifs <- function(selected_dist_real_features, type, values){
	summary_selected_dist_real_features <- selected_dist_real_features %>%
		group_by(motif) %>%
		summarize(n=n(), .groups="drop_last") %>%
		mutate(motif=as.character(motif), len=nchar(motif))
	selected_dist_real_features_bipartite <- summary_selected_dist_real_features %>%
		filter(len>8)
	selected_dist_real_features_monopartite <- summary_selected_dist_real_features %>%
		filter(len<8)
	if(type=="only"){
		selected_dist_real_features_bipartite <- selected_dist_real_features_bipartite %>%
			filter(n %in% values)
	}else if(type=="at_least"){
		selected_dist_real_features_bipartite <- selected_dist_real_features_bipartite %>%
			filter(n>=values)
	}
	selected_dist_real_features <- subset(selected_dist_real_features, motif %in% c(selected_dist_real_features_monopartite$motif, selected_dist_real_features_bipartite$motif))

	return(selected_dist_real_features)
}

filtering.length.motifs <- function(selected_dist_real_features, type, values){
	summary_selected_dist_real_features <- selected_dist_real_features %>%
		group_by(motif) %>%
		summarize(n=n(), .groups="drop_last") %>%
		mutate(motif=as.character(motif), len=nchar(motif))
	if(type=="max"){
		selected_dist_real_features_lenght <- summary_selected_dist_real_features %>%
			filter(len<values)
	}
	selected_dist_real_features <- subset(selected_dist_real_features, motif %in% selected_dist_real_features_lenght$motif)

	return(selected_dist_real_features)
}

filtering.sigOnly.features <- function(scored_dist_real_filtered_motifs, selected_dist_real_features, type){
	scored_dist_real_filtered_motifs <- scored_dist_real_filtered_motifs %>%
		mutate(feature_name=paste0(motif,"_",distance_motif))
	if(type=="all"){
		scored_dist_real_filtered_motifs <- scored_dist_real_filtered_motifs %>%
			filter(feature_name %in% selected_dist_real_features$feature_name) %>%
			dplyr::select(-c(feature_name))
	}else if(type=="bipartite"){
		scored_dist_real_filtered_motifs <- scored_dist_real_filtered_motifs %>%
			filter(nchar(as.character(motif))<8 | feature_name %in% selected_dist_real_features$feature_name) %>%
			dplyr::select(-c(feature_name))
	}

	return(scored_dist_real_filtered_motifs)
}

show.motifs.scores.distribution <- function(scored_dist_real_motifs, list_motifs, genome_id, name_output, score_threshold_range, minOcc){
	pdf(name_output, width=7, height=7)
	foreach(motif_of_interest=list_motifs) %do% {
		show.motif.scores.distribution(motif_of_interest, scored_dist_real_motifs, genome_id, score_threshold_range, minOcc)
	}
	dev.off()
}

check.motifs.selection <- function(selected_dist_real_features, list_expected_motifs){
	list_selected_motifs <- unique(selected_dist_real_features$motif)

	recovery <- foreach(motif_of_interest=list_expected_motifs, .combine=rbind) %do% {
		matching_selected_motifs <- list_selected_motifs[grepl(convert.motif.grep(motif_of_interest),list_selected_motifs)]

		return(data.frame(motif=motif_of_interest, n_matches=length(matching_selected_motifs)))
	}

	accuracy <- selected_dist_real_features
	accuracy$n_matches <- 0
	stiffle <- foreach(motif_of_interest=list_expected_motifs, .combine=rbind) %do% {
		is_matches <- grepl(convert.motif.grep(motif_of_interest),accuracy$motif)
		accuracy$n_matches[is_matches] <- accuracy$n_matches[is_matches] + 1

		return(NA)
	}

	recovery <- recovery %>%
		mutate(is_found=ifelse(n_matches>0, "Yes", "No")) %>%
		group_by(is_found) %>%
		summarize(score=n(), .groups="drop_last") %>%
		spread(is_found, score) %>%
		mutate(pseudo_comp=(Yes/(Yes+No))*100)
	print(paste0(round(recovery$pseudo_comp,2),"% of expected motifs covered by selected motifs (",recovery$Yes,"/",recovery$Yes+recovery$No,")."))

	sum_feat_accuracy <- accuracy %>%
		mutate(is_explained=ifelse(n_matches>0, "Yes", "No")) %>%
		group_by(is_explained) %>%
		summarize(score=n(), .groups="drop_last") %>%
		spread(is_explained, score) %>%
		mutate(pseudo_acc=(Yes/(Yes+No))*100)
	sum_mot_accuracy <- accuracy %>%
		group_by(motif, n_matches) %>%
		summarize(n=n(), .groups="drop_last") %>%
		mutate(is_explained=ifelse(n_matches>0, "Yes", "No")) %>%
		group_by(is_explained) %>%
		summarize(score=n(), .groups="drop_last") %>%
		spread(is_explained, score) %>%
		mutate(pseudo_acc=(Yes/(Yes+No))*100)
	print(paste0(round(sum_feat_accuracy$pseudo_acc,2),"% of selected features were expected (",sum_feat_accuracy$Yes,"/",sum_feat_accuracy$Yes+sum_feat_accuracy$No,")."))
	print(paste0(round(sum_mot_accuracy$pseudo_acc,2),"% of selected motifs were expected (",sum_mot_accuracy$Yes,"/",sum_mot_accuracy$Yes+sum_mot_accuracy$No,")."))
}

#  ___  ___      _   _  __      ______     _            _   _             
#  |  \/  |     | | (_)/ _|     |  _  \   | |          | | (_)            
#  | .  . | ___ | |_ _| |_ ___  | | | |___| |_ ___  ___| |_ _  ___  _ __  
#  | |\/| |/ _ \| __| |  _/ __| | | | / _ \ __/ _ \/ __| __| |/ _ \| '_ \ 
#  | |  | | (_) | |_| | | \__ \ | |/ /  __/ ||  __/ (__| |_| | (_) | | | |
#  \_|  |_/\___/ \__|_|_| |___/ |___/ \___|\__\___|\___|\__|_|\___/|_| |_|
#                                                                         
#                                                                         

contigs.motif.detection <- function(methylation_signal, list_contig_name, path_metagenome, detection_params, iupac_nc, contigs_group_name=NA, discovered_motifs_contig=NULL){
	methylation_signal_contig <- subset(methylation_signal, contig %in% list_contig_name)

	if(length(list_contig_name)==1){
		output_dir_name <- paste0("./motifs_detection/",list_contig_name)
		if(!dir.exists(output_dir_name)){
			dir.create(output_dir_name)
		}
	}else{
		if(!is.na(contigs_group_name)){
			output_dir_name <- paste0("./motifs_detection/",contigs_group_name)
			if(!dir.exists(output_dir_name)){
				dir.create(output_dir_name)
			}
		}else{
			stop("Missing contigs group name.")
		}
	}
	discovered_motifs_contig <- wrapper.motif.detection(methylation_signal_contig, path_metagenome, output_dir_name, detection_params, iupac_nc, discovered_motifs_contig)

	return(discovered_motifs_contig)
}

contigs.refine.motif <- function(motif_of_interest, discovered_motifs_contig, methylation_signal, list_contig_name, path_metagenome, detection_params, iupac_nc, contigs_group_name=NA){
	methylation_signal_contig <- subset(methylation_signal, contig %in% list_contig_name)

	if(length(list_contig_name)==1){
		output_dir_name <- paste0("./motifs_detection/",list_contig_name)
	}else{
		if(!is.na(contigs_group_name)){
			output_dir_name <- paste0("./motifs_detection/",contigs_group_name)
		}else{
			stop("Missing contigs group name.")
		}
	}
	refine.motif(motif_of_interest, discovered_motifs_contig, methylation_signal_contig, path_metagenome, detection_params$nbCPU, detection_params$seq_params, iupac_nc, output_dir_name, "tmp")
}

confirm.bin.identity <- function(discovered_motifs, methylation_signal, list_contigs, path_metagenome, detection_params, iupac_nc, output_name){
	stifle <- foreach(motif=discovered_motifs) %do% {
		contigs.refine.motif(motif, discovered_motifs[!discovered_motifs %in% motif], methylation_signal, list_contigs, path_metagenome, detection_params, iupac_nc, output_name)
	
		return(NA)
	}
}

find.signature.center.bin <- function(methylation_signal, list_contigs, discovered_motifs, path_metagenome, nbCPU, seq_params, iupac_nc, signal_margin, signal_len, base_name){
	methylation_signal_contig <- subset(methylation_signal, contig %in% list_contigs)
	methylation_signal_contig_summary <- data.frame(motif=discovered_motifs, mod_pos=rep(NA,length(discovered_motifs)), stringsAsFactors=FALSE)

	motif_center_summary <- find.signature.center(methylation_signal_contig, methylation_signal_contig_summary, path_metagenome, nbCPU, seq_params, iupac_nc, signal_margin, signal_len, base_name)

	return(motif_center_summary)
}

classify.detected.motifs.bin <- function(methylation_signal, list_contigs, bin_id, motif_center_summary, model, path_metagenome, min_cov, keepIsolated, iupac_nc, nbCPU){
	methylation_signal_contig <- subset(methylation_signal, contig %in% list_contigs)

	classification_results <- classify.detected.motifs(methylation_signal_contig, bin_id, motif_center_summary, model, path_metagenome, min_cov, keepIsolated, iupac_nc, nbCPU)

	return(classification_results)
}

#  ___  ____                                  _     _        ______     _            _   _             
#  |  \/  (_)                                | |   | |       |  _  \   | |          | | (_)            
#  | .  . |_ ___  __ _ ___ ___  ___ _ __ ___ | |__ | |_   _  | | | |___| |_ ___  ___| |_ _  ___  _ __  
#  | |\/| | / __|/ _` / __/ __|/ _ \ '_ ` _ \| '_ \| | | | | | | | / _ \ __/ _ \/ __| __| |/ _ \| '_ \ 
#  | |  | | \__ \ (_| \__ \__ \  __/ | | | | | |_) | | |_| | | |/ /  __/ ||  __/ (__| |_| | (_) | | | |
#  \_|  |_/_|___/\__,_|___/___/\___|_| |_| |_|_.__/|_|\__, | |___/ \___|\__\___|\___|\__|_|\___/|_| |_|
#                                                      __/ |                                           
#                                                     |___/                                            

score.contig.motif <- function(motif_to_score, discovered_motifs, contig_name, methylation_signal, path_metagenome, window_size, scoring_type="top"){
	expected_signal_left <- -6
	expected_signal_right <- -1
	signal_margin <- 4
	min_cov <- 0

	methylation_signal_contig <- subset(methylation_signal, contig==contig_name)

	# Detect misassemblies
	tagged_motif <- tag.motifs(motif_to_score, discovered_motifs, methylation_signal_contig, path_metagenome, iupac_nc, expected_signal_left, expected_signal_right, signal_margin, min_cov)

	# Score occurrences
	if(scoring_type=="top_win"){
		min_distance <- min(tagged_motif$distance, na.rm=TRUE)
		max_distance <- max(tagged_motif$distance, na.rm=TRUE)
		scored_motif <- tagged_motif %>%
			group_by(contig, pos_motif, dir, motif) %>%
			mutate(cov_wga=mean(N_wga, na.rm=TRUE), cov_nat=mean(N_nat, na.rm=TRUE)) %>%
			mutate(score=rollapplyr(abs(mean_diff), 6, mean, partial=TRUE, fill=NA, align="center", na.rm=TRUE)) %>%
			filter(distance %in% seq(min_distance + 2, max_distance - 3)) %>%
			filter(score==max(score, na.rm=TRUE)) %>%
			filter(row_number()==1) # Remove potential/rare ties; Keep first distance only
	}else if(scoring_type=="top"){
		scored_motif <- tagged_motif %>%
			group_by(contig, pos_motif, dir, motif) %>%
			mutate(cov_wga=mean(N_wga, na.rm=TRUE), cov_nat=mean(N_nat, na.rm=TRUE)) %>%
			top_n(n=6, wt=abs(mean_diff)) %>% # Keep only 6 max position to be comparable between motifs
			arrange(contig, pos_motif, dir, motif, desc(abs(mean_diff))) %>%
			group_by(contig, pos_motif, dir, motif) %>%
			filter(row_number()<=6) %>% # Remove potential/rare ties
			summarize(score=mean(abs(mean_diff), na.rm=TRUE), cov_wga=unique(cov_wga), cov_nat=unique(cov_nat), .groups="drop_last")
	}

	contig_ranges <- range(scored_motif$pos_motif)

	contig_length <- contig_ranges[2] - (contig_ranges[1] - 1)
	if(is.na(window_size)){
		nb_occurrences <- nrow(scored_motif)
		nb_avg_motifs_per_window <- 100

		motif_span <- min(nb_avg_motifs_per_window/nb_occurrences, 0.1) # Fix maximum at 10% contig's occurrences
	}else{
		motif_span <- window_size/contig_length
	}
	fit_score <- loess(score ~ pos_motif, data=scored_motif, span=motif_span)
	scored_motif$fitted_score <- predict(fit_score)
	fit_cov_wga <- loess(cov_wga ~ pos_motif, data=scored_motif, span=motif_span)
	scored_motif$fitted_cov_wga <- predict(fit_cov_wga)
	fit_cov_nat <- loess(cov_nat ~ pos_motif, data=scored_motif, span=motif_span)
	scored_motif$fitted_cov_nat <- predict(fit_cov_nat)

	return(scored_motif)
}

compute.contig.composition <- function(contig_name, path_metagenome, window_size, step_size, nbCPU){
	g_seqs <- read.fasta(path_metagenome)
	g_seq <- g_seqs[attr(g_seqs,"name")==contig_name]

	len_contig <- length(g_seq[[1]])
	starts_windows <- seq(1, len_contig - window_size, step_size)
	relative_position <- (window_size/2) - 1
	registerDoMC(nbCPU)
	contig_content <- foreach(start_window=starts_windows, .combine=rbind) %dopar% {
		nucleotide_n <- seqinr::count(g_seq[[1]][seq(start_window,start_window + (window_size - 1))],1)
		result <- data.frame((nucleotide_n * 100)/window_size)
		result$start <- start_window + relative_position

		return(result)
	}
	registerDoSEQ()

	colnames(contig_content) <- c("Base","Frequency","Position")
	contig_content$Contig_Name <- as.factor(paste0(gsub(">","",attr(g_seq[[1]], "Annot"))," (",len_contig," bp)"))
	contig_content$Base <- mapvalues(contig_content$Base, c("a","c","g","t"), c("A","C","G","T"))

	contig_content <- dcast(contig_content, Position + Contig_Name ~ Base, value.var="Frequency")

	return(contig_content)
}

detection.contig.misassemblies <- function(motif_to_score, discovered_motifs, contig_name, base_name, methylation_signal, path_metagenome, scoring_type, window_size_motif, window_size_content, step_size, nbCPU){
	scored_motif <- score.contig.motif(motif_to_score, discovered_motifs, contig_name, methylation_signal, path_metagenome, window_size_motif, scoring_type)

	gp_scores <- ggplot(scored_motif) +
		geom_point(aes(x=pos_motif, y=score)) +
		geom_line(aes(x=pos_motif, y=fitted_score), col=2, size=1.5) +
		scale_x_continuous(expand=c(0,0)) +
		coord_cartesian(ylim=c(0,max(scored_motif$score))) +
		labs(title="Motif Scores", subtitle=paste0(motif_to_score," motif; Loess with ",window_size_motif," bp windows.")) +
		labs(x="Contig position", y="Motif Score (abs top 6)")

	gp_coverages <- ggplot(scored_motif) +
		geom_line(aes(x=pos_motif, y=fitted_cov_nat, col="NAT"), size=1.5) +
		geom_line(aes(x=pos_motif, y=fitted_cov_wga, col="WGA"), size=1.5) +
		scale_x_continuous(expand=c(0,0)) +
		coord_cartesian(ylim=c(0,max(scored_motif$fitted_cov_nat, scored_motif$fitted_cov_wga))) +
		scale_colour_manual(name="Dataset", values=c("NAT"="#0059D3","WGA"="#D35900")) +
		labs(title="Contig coverage", subtitle=paste0(motif_to_score," motif; Loess with ",window_size_motif," bp windows.")) +
		labs(x="Contig position", y="Coverage")

	contig_content <- compute.contig.composition(contig_name, path_metagenome, window_size_content, step_size, nbCPU)

	gp_content <- ggplot(contig_content) +
		geom_line(aes(x=Position, y=G, col="G"), size=0.1) +
		geom_line(aes(x=Position, y=C, col="C"), size=0.1) +
		geom_line(aes(x=Position, y=A, col="A"), size=0.1) +
		geom_line(aes(x=Position, y=T, col="T"), size=0.1) +
		scale_x_continuous(expand=c(0,0)) +
		scale_colour_manual(name="Base", values=c("A"="#0059D3","C"="#D35900","G"="#EE1010","T"="#01AB25")) +
		labs(title="Contig content", subtitle=paste0("Sliding window of ",window_size_content," bp with ",step_size," steps size.")) +
		labs(x="Contig position", y="Base %") +
		guides(colour=guide_legend(override.aes=list(size=2)))

	pdf(paste0("Misassembly_detection_",base_name,"_",scoring_type,"_v1.pdf"), width=20, height=10) # TODO tweak
	ggarrange(gp_scores, gp_coverages, gp_content, nrow=3, newpage=FALSE, top=paste0("Detection of misassemblies in ",contig_name," contig"))
	dev.off()
}

# Can generate warning with rare motifs
detection.contig.misassemblies.binMotifs <- function(discovered_motifs, contig_name, base_name, methylation_signal, path_metagenome, scoring_type, window_size_motif, window_size_content, step_size){
	stifle <- foreach(motif_to_score=discovered_motifs) %do% {
		detection.contig.misassemblies(motif_to_score, discovered_motifs[!grepl(motif_to_score, discovered_motifs)], contig_name, paste0(contig_name,"_",motif_to_score,"_",base_name), methylation_signal, path_metagenome, scoring_type, window_size_motif, window_size_content, step_size)

		return(NA)
	}
}

detection.contig.misassemblies.simple <- function(discovered_motifs, list_contig_name, base_name, methylation_signal, path_metagenome, scoring_type, window_size_motif, additionnal_motifs=NA, rename_contigs=NA, name_motifs_legend=NA){
	discovered_motifs <- unique(discovered_motifs) # Remove duplicate if used multiple lists

	registerDoMC(nbCPU)
	scored_motif <- foreach(contig_name=list_contig_name, .combine=rbind) %do% {
		contig_scored_motif <- foreach(motif_to_score=discovered_motifs, .combine=rbind) %dopar% {
			if(any(is.na(additionnal_motifs))){
				background_motifs <- discovered_motifs[!grepl(motif_to_score, discovered_motifs)]
				subset_scored_motif <- score.contig.motif(motif_to_score, background_motifs, contig_name, methylation_signal, path_metagenome, window_size_motif, scoring_type)
			}else{
				background_motifs <- unique(additionnal_motifs, discovered_motifs)
				background_motifs <- background_motifs[!grepl(motif_to_score, background_motifs)] # Remove exact or more precise motifs
				background_motifs <- background_motifs[!sapply(background_motifs,function(x) grepl(x,motif_to_score))] # Remove "master motifs" (less precise)
				subset_scored_motif <- score.contig.motif(motif_to_score, background_motifs, contig_name, methylation_signal, path_metagenome, window_size_motif, scoring_type)
			}
			subset_scored_motif$motif <- as.character(subset_scored_motif$motif)

			return(subset_scored_motif)
		}

		return(contig_scored_motif)
	}
	registerDoSEQ()
	scored_motif$motif <- as.factor(scored_motif$motif)

	if(!any(is.na(rename_contigs))){
		scored_motif$contig <- droplevels(scored_motif$contig)
		scored_motif$contig <- mapvalues(scored_motif$contig, from=list_contig_name, to=rename_contigs)
	}

	gp_scores <- ggplot(scored_motif) +
		geom_line(aes(x=pos_motif, y=fitted_score, col=motif)) +
		geom_hline(aes(yintercept=1, lty="dashed")) +
		facet_wrap(~contig, ncol=1, scales="free_x") +
		scale_x_continuous(expand=c(0,0)) +
		scale_linetype_identity(name="", labels=c("Baseline"), guide="legend") +
		# coord_cartesian(ylim=ifelse(scoring_type=="top",c(0.75,3),c(0,3))) + # Should cover most cases
		labs(title="Evaluate contig misassemblies") +
		labs(subtitle=paste0("Loess with ",ifelse(is.na(window_size_motif),paste0("adaptative windows."),paste0(window_size_motif," bp windows.")))) +
		labs(x="Contig position", y="Smoothed pA differences") +
		labs(col="Motifs") +
		theme_bw() +
		theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
		theme(strip.background=element_rect(fill="white"), strip.text=element_text(size=12)) +
		guides(colour=guide_legend(override.aes=list(shape=15, size=4))) +
		theme(legend.text=element_text(size=10))

	if(scoring_type=="top"){
		gp_scores <- gp_scores +
			coord_cartesian(ylim=c(0.75,3)) # Should cover most cases
	}else if(scoring_type=="top_win"){
		gp_scores <- gp_scores +
			coord_cartesian(ylim=c(0.5,3)) # Should cover most cases
	}

	if(!is.na(name_motifs_legend)){
		gp_scores <- gp_scores +
			labs(col=name_motifs_legend)
	}

	# Add columns when many motifs
	if(length(discovered_motifs)>8){
		gp_scores <- gp_scores +
			guides(col=guide_legend(ncol=ceiling(length(discovered_motifs)/8), override.aes=list(shape=15, size=4)))
	}

	pdf(paste0("Misassembly_detection_simple_",base_name,"_",scoring_type,"_v1.pdf"), width=10, height=0.5 + 2.5*length(list_contig_name))
	print(gp_scores)
	dev.off()
}

detection.contig.misassemblies.complex <- function(list_discovered_motifs, list_contig_name, base_name, methylation_signal, path_metagenome, scoring_type, window_size_motif, additionnal_motifs=NA, rename_contigs=NA){
	list_gp_scores <- foreach(idx_discovered_motifs=seq(1,length(list_discovered_motifs))) %do% {
		elem_discovered_motifs <- list_discovered_motifs[idx_discovered_motifs]
		discovered_motifs <- elem_discovered_motifs[[1]]
		discovered_motifs <- unique(discovered_motifs) # Remove duplicate if used multiple lists

		registerDoMC(nbCPU)
		scored_motif <- foreach(contig_name=list_contig_name, .combine=rbind) %do% {
			contig_scored_motif <- foreach(motif_to_score=discovered_motifs, .combine=rbind) %dopar% {
				if(any(is.na(additionnal_motifs))){
					background_motifs <- discovered_motifs[!grepl(motif_to_score, discovered_motifs)]
					subset_scored_motif <- score.contig.motif(motif_to_score, background_motifs, contig_name, methylation_signal, path_metagenome, window_size_motif, scoring_type)
				}else{
					background_motifs <- unique(additionnal_motifs, discovered_motifs)
					background_motifs <- background_motifs[!grepl(motif_to_score, background_motifs)] # Remove exact or more precise motifs
					background_motifs <- background_motifs[!sapply(background_motifs,function(x) grepl(x,motif_to_score))] # Remove "master motifs" (less precise)
					subset_scored_motif <- score.contig.motif(motif_to_score, background_motifs, contig_name, methylation_signal, path_metagenome, window_size_motif, scoring_type)
				}
				subset_scored_motif$motif <- as.character(subset_scored_motif$motif)

				return(subset_scored_motif)
			}

			return(contig_scored_motif)
		}
		registerDoSEQ()
		scored_motif$motif <- as.factor(scored_motif$motif)

		if(!any(is.na(rename_contigs))){
			scored_motif$contig <- droplevels(scored_motif$contig)
			scored_motif$contig <- mapvalues(scored_motif$contig, from=list_contig_name, to=rename_contigs)
		}

		gp_scores <- ggplot(scored_motif) +
			geom_line(aes(x=pos_motif, y=fitted_score, col=motif)) +
			# geom_hline(aes(yintercept=1, lty="dashed")) +
			# facet_wrap(~contig, ncol=1, scales="free_x") +
			scale_x_continuous(expand=c(0,0)) +
			# coord_cartesian(ylim=ifelse(scoring_type=="top",c(0.75,3),c(0,3))) + # Should cover most cases
			labs(x=NULL, y=NULL) +
			labs(col=paste0("Methylation motifs\nfrom ",names(elem_discovered_motifs))) +
			theme_bw() +
			theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
			guides(colour=guide_legend(override.aes=list(shape=15, size=4))) +
			theme(legend.text=element_text(size=10))

		if(scoring_type=="top"){
			gp_scores <- gp_scores +
				coord_cartesian(ylim=c(0.75,3)) # Should cover most cases
		}else if(scoring_type=="top_win"){
			gp_scores <- gp_scores +
				coord_cartesian(ylim=c(0.5,3)) # Should cover most cases
		}

		if(idx_discovered_motifs!=1){
			gp_scores <- gp_scores +
				geom_hline(aes(yintercept=1), lty="dashed") +
				guides(linetype=FALSE)
		}

		if(idx_discovered_motifs==1){
			gp_scores <- gp_scores +
				geom_hline(aes(yintercept=1, lty="dashed")) +
				scale_linetype_identity(name="", labels=c("Baseline"), guide="legend")
		}

		if(idx_discovered_motifs!=length(list_discovered_motifs)){
			gp_scores <- gp_scores +
				theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
		}
		if(idx_discovered_motifs==ceiling(length(list_discovered_motifs)/2)){
			gp_scores <- gp_scores +
				labs(y="Smoothed pA differences") +
				theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
		}

		if(idx_discovered_motifs==length(list_discovered_motifs)){
			gp_scores <- gp_scores +
				labs(x="Contig position")
		}

		# Add columns when many motifs
		if(length(discovered_motifs)>8){
			gp_scores <- gp_scores +
				guides(col=guide_legend(ncol=ceiling(length(discovered_motifs)/8), override.aes=list(shape=15, size=4)))
		}

		return(gp_scores)
	}

	pdf(paste0("Misassembly_detection_complex_",base_name,"_",scoring_type,"_v1.pdf"), width=10, height=0.5 + 3*length(list_discovered_motifs))
	ggarrange(plots=list_gp_scores, newpage=FALSE, ncol=1)
	dev.off()
}

#   _____                                        _____                                 _ _   _             
#  /  ___|                                      /  __ \                               (_) | (_)            
#  \ `--.  ___  __ _ _   _  ___ _ __   ___ ___  | /  \/ ___  _ __ ___  _ __   ___  ___ _| |_ _  ___  _ __  
#   `--. \/ _ \/ _` | | | |/ _ \ '_ \ / __/ _ \ | |    / _ \| '_ ` _ \| '_ \ / _ \/ __| | __| |/ _ \| '_ \ 
#  /\__/ /  __/ (_| | |_| |  __/ | | | (_|  __/ | \__/\ (_) | | | | | | |_) | (_) \__ \ | |_| | (_) | | | |
#  \____/ \___|\__, |\__,_|\___|_| |_|\___\___|  \____/\___/|_| |_| |_| .__/ \___/|___/_|\__|_|\___/|_| |_|
#                 | |                                                 | |                                  
#                 |_|                                                 |_|                                  

# Do I need to combine with reverse complement? Alneberg justify by biderectionnal sequencing.
# rowsum(kmer_frequency)

normalize.kmer.frequency <- function(x){

	return(log(x/sum(x)))
}

prepare.kmer.frequency <- function(kmer_frequency, pseudocount){
	kmer_frequency <- kmer_frequency + pseudocount # Add pseudocount
	kmer_frequency <- apply(kmer_frequency, 1, normalize.kmer.frequency)

	return(t(kmer_frequency))
}

tsne.kmer.frequency <- function(kmer_frequency, sequence_metagenome, contigs_info, grouping_var, contig_weight_unit=NA, seed=101){
	kmer_frequency <- prepare.kmer.frequency(kmer_frequency, 1)

	# Weight contigs if requested
	if(!is.na(contig_weight_unit)){
		contigs_info <- contigs_info %>%
			mutate(weight=ceiling(length/contig_weight_unit))

		element_to_repeat <- match(rep(contigs_info$contig, times=contigs_info$weight), names(sequence_metagenome))
		kmer_frequency <- kmer_frequency[element_to_repeat[!is.na(element_to_repeat)],] # Some contigs are missing after filter(nb_occurrence>=5)
		contig_name_to_add <- names(sequence_metagenome)[element_to_repeat[!is.na(element_to_repeat)]]
	}else{
		contig_name_to_add <- names(sequence_metagenome)
	}

	set.seed(seed)
	tsne_kmer_frequency_data <- as.data.frame(Rtsne(kmer_frequency, check_duplicates=FALSE)$Y)

	tsne_kmer_frequency_data <- tsne_kmer_frequency_data %>%
		mutate(contig=as.factor(contig_name_to_add)) %>%
		plyr::rename(c('V1'='tSNE_1','V2'='tSNE_2')) %>%
		mutate(genome=as.factor(str_split(contig, " ", simplify=TRUE)[,1])) %>%
		mutate(contig_length=contigs_info$length[match(contig, contigs_info$contig)]) %>%
		mutate(id=paste0(contigs_info$id[match(contig, contigs_info$contig)]))
	
	if(!is.na(contig_weight_unit)){
		tsne_kmer_frequency_data <- tsne_kmer_frequency_data %>%
			group_by(contig) %>%
			mutate(centroid_V1=mean(tSNE_1),centroid_V2=mean(tSNE_2)) %>%
			distinct(contig, .keep_all=TRUE) %>%
			mutate(tSNE_1=centroid_V1,tSNE_2=centroid_V2) %>%
			dplyr::select(-c(centroid_V1, centroid_V2))
	}

	tsne_kmer_frequency_data <- tsne_kmer_frequency_data %>%
		group_by_(.dots=grouping_var) %>%
		mutate(n_contig=n()) %>%
		ungroup() %>%
		mutate(id=paste0(id," (n=",n_contig,")"))

	return(tsne_kmer_frequency_data)
}

plot.tsne.kmer.frequency <- function(tsne_kmer_frequency_data, base_name, is_weighted=FALSE){
	gp <- ggplot(tsne_kmer_frequency_data) +
		geom_point(aes(x=tSNE_1 ,y=tSNE_2, col=id, size=contig_length), shape=1) +
		scale_size_continuous(breaks=c(10000,100000,1000000), labels=c(10000,100000,1000000), range=c(0.1,9)) +
		labs(x="t-SNE 1", y="t-SNE 2", colour="Genome of Origin", size="Contigs Length") +
		guides(color=guide_legend(order=1), size=guide_legend(order=0)) # Decreasing order

	if(length(unique(tsne_kmer_frequency_data$id))>20){ # If too many colors plot separate legend
		gp_with_legend <- gp
		pdf(file=NULL)
		gp_legend <- gglegend(gp_with_legend)
		dev.off()
		gp <- gp +
			theme(legend.position="none")
		pdf(paste0("Contigs_composition_tsne_",base_name,"_v2.pdf"), width=9)
		print(gp)
		grid.newpage()
		grid.draw(gp_legend)
		dev.off()
	}else{ # Plot with legend
		pdf(paste0("Contigs_composition_tsne_",base_name,"_v2.pdf"), width=9)
		print(gp)
		dev.off()
	}

	if(is_weighted){
		gp <- ggplot(tsne_kmer_frequency_data) +
			geom_point(aes(x=tSNE_1, y=tSNE_2, col=contig_length)) +
			facet_wrap(~id) +
			labs(x="t-SNE 1", y="t-SNE 2", colour="Contig of Origin") +
			theme(legend.position="none")
		pdf(paste0("Contigs_composition_tsne_",base_name,"_detailed_v2.pdf"), height=9, width=9)
		print(gp)
		dev.off()
	}
}


