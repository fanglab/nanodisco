#!/usr/bin/env bash

# Minimum input command
# nanodisco difference -nj 1 -nc 1 -p 1 -f 1 -l 2 -i . -o . -w wga -n nat -r ref.fa
# Example command using provide dataset with get_data_bacteria.
# nanodisco difference -nj 1 -nc 1 -p 2 -x batch -f 281 -l 290 -i ./ -o ./processed_example -w EC_WGA -n EC_NAT -r /home/nanodisco/reference/Ecoli_K12_MG1655_ATCC47076.fasta -b revc -a 1.5 -z 2 -e 5 -j noAddSupp -k 0 -t /usr/bin/r9.4_450bps.nucleotide.6mer.template.model

###
# Wrapper for computing current differences
###
# 1. Control input validity 
# 2. Divide queried chunks across available ressources
# 3. Computing current differences

## Describe parameters
nb_jobs_def="Specify how many jobs to run in parallel (affect CPU and memory usage)."
nb_chunks_def="Specify how many chunks to process in a row (affect memory usage)."
nb_threads_def="Specify how many threads to use per job."
first_chunk_def="Specify which is the 1st chunk you want to process."
last_chunk_def="Specify which is the last chunk you want to process."
path_input_def="Specify the path of input data (folder with .fasta, .bam, and .bam.bai)."
path_output_def="Specify the path of output results (current_difference file and logs)."
wga_name_def="Specify the name of WGA sample (same than for nanodisco preprocess -s <sample_name>)."
nat_name_def="Specify the name of native sample (same than for nanodisco preprocess -s <sample_name>)."
genome_def="Specify the path of reference genome (.fasta)."
sig_norm_def="Specify the type of correction for strand bias (ori|revc, default is revc)."
IQR_factor_def="Specify the IQR factor for outliers removal (0.0 to skip; default is 1.5; smaller is harsher)."
normalized_def="Specify the type of normalization (0:nanopolish only, 1: np+lm, 2: np+rlm; default is 2)."
min_coverage_def="Specify the minimum number of events per position (default is 5)."
map_type_def="Specify the type of filtering for mapping (default is noAddSupp, see realign_events.sh for other types)."
min_read_length_def="Specify the minimum mapped read length (default is 0)."
path_ont_model_def="Specify the path to the normalization model (default is r9.4_450bps.nucleotide.6mer.template.model from nanopolish)."
exec_type_def="Specify which type of exectution is required (batch or seq; default is batch; seq is for developpement only)."

## Default parameters
# Best parameters for individual bacteria processing
sig_norm="revc"                # type of correction for strand bias
IQR_factor=1.5                 # value k of Tukey’s outlier removal method
normalized=2                   # type of normalization (1: lm, 2: rlm)
min_coverage=5                 # minimum coverage to compute current difference and statistics
map_type="noAddSupp"           # type of alignment filtering
min_read_length=0              # minimum mapped read length
path_ont_model='/usr/bin/r9.4_450bps.nucleotide.6mer.template.model' # Path to 6mer model for event signal normalization
exec_type="batch"              # seq is for developpement only

# Read options
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -nj|--nb_jobs)
      if [[ $2 =~ ^[\-0-9]+$ ]] && (( $2 > 0)); then
        nb_jobs="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-nj/--nb_jobs positive integer required." >&2
        echo -e "\t"$nb_jobs_def
        exit 5
      fi
      ;;
    -nc|--nb_chunks)
      if [[ $2 =~ ^[\-0-9]+$ ]] && (( $2 > 0)); then
        nb_chunks="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-nc/--nb_chunks positive integer required." >&2
        echo -e "\t"$nb_chunks_def
        exit 5
      fi
      ;;
    -p|--nb_threads)
      if [[ $2 =~ ^[\-0-9]+$ ]] && (( $2 > 0)); then
        nb_threads="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-p/--nb_threads positive integer required." >&2
        echo -e "\t"$nb_threads_def
        exit 5
      fi
      ;;
    -x|--exec_type)
      if [[ $2 =~ (batch|seq) ]]; then
        exec_type="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-x/--exec_type can be batch or seq." >&2
        echo -e "\t"$exec_type_def
        exit 5
      fi
      ;;
    -f|--first_chunk)
      if [[ $2 =~ ^[\-0-9]+$ ]] && (( $2 > 0)); then
        first_chunk="$2"
        shift # pass argument
        shift # pass value
      # else
        # echo "-f/--first_chunk positive integer required." >&2
        # echo -e "\t"$first_chunk_def
        # exit 5
      fi
      ;;
    -l|--last_chunk)
      if [[ $2 =~ ^[\-0-9]+$ ]] && (( $2 > 0)); then
        last_chunk="$2"
        shift # pass argument
        shift # pass value
      # else
        # echo "-l/--last_chunk positive integer required." >&2
        # echo -e "\t"$last_chunk_def
        # exit 5
      fi
      ;;
    -i|--path_input)
      if [[ -d $2 ]]; then
        path_input="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-i/--path_input directory doesn't exist." >&2
        echo -e "\t"$path_input_def
        exit 5
      fi
      ;;
    -o|--path_output)
      if [[ -d $2 ]]; then
        echo "-o/--path_output directory already exist. Some data could be overwrite." >&2
      fi
      path_output="$2"
      shift # pass argument
      shift # pass value
      ;;
    -w|--wga_name)
      wga_name="$2"
      shift # pass argument
      shift # pass value
      ;;
    -n|--nat_name)
      nat_name="$2"
      shift # pass argument
      shift # pass value
      ;;
    -r|--genome)
      if [[ -f $2 ]]; then
        genome="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-r/--genome file doesn't exist." >&2
        echo -e "\t"$genome_def >&2
        exit 5
      fi
      ;;
    -b|--sig_norm)
      if [[ $2 =~ (revc|ori) ]]; then
        sig_norm="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-b/--sig_norm can be revc or ori only." >&2
        echo -e "\t"$sig_norm_def >&2
        exit 5
      fi
      ;;
    -a|--IQR_factor)
      if [[ $2 =~ ^[\-0-9]+\.?[0-9]*$ ]] && (( $(awk 'BEGIN {print ("'$2'" >= 0)}') )); then
        IQR_factor="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-a/--IQR_factor 0.0 or positive float required." >&2
        echo -e "\t"$IQR_factor_def >&2
        exit 5
      fi
      ;;
    -z|--normalized)
      if [[ $2 =~ ^[012]$ ]]; then
        normalized="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-z/--normalized can be 0, 1, or 2." >&2
        echo -e "\t"$normalized_def >&2
        exit 5
      fi
      ;;
    -e|--min_coverage)
      if [[ $2 =~ ^[\-0-9]+$ ]] && (( $2 > 0)); then
        min_coverage="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-e/--min_coverage positive integer required." >&2
        echo -e "\t"$min_coverage_def >&2
        exit 5
      fi
      ;;
    -j|--map_type)
      if [[ $2 =~ (all|noMulti|noAddSupp|noSupp|noMultiAddSupp|none) ]]; then
        map_type="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-j/--map_type isn't recognized." >&2
        echo -e "\t"$map_type_def >&2
        exit 5
      fi
      ;;
    -k|--min_read_length)
      if [[ $2 =~ ^[\-0-9]+$ ]] && (( $2 >= 0)); then
        min_read_length="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-k/--min_read_length positive integer required." >&2
        echo -e "\t"$min_read_length_def >&2
        exit 5
      fi
      ;;
    -t|--path_ont_model)
      if [[ -f $2 ]]; then
        path_ont_model="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-t/--path_ont_model file doesn't exist." >&2
        echo -e "\t"$path_ont_model_def >&2
        exit 5
      fi
      ;;
    -h|--help) # Print help
      echo "Usage: nanodisco difference -nj <nb_jobs> -nc <nb_chunks> -p <nb_threads> [-f <first_chunk> -l <last_chunk>] -i <path_input> -o <path_output> -w <name_WGA> -n <name_native> -r <path_genome> [+ advanced parameters]" >&2
      echo -e "\tAdditional information can be found in our GitHub repository." >&2
      exit 3
      ;;
    *) # unknown option
      echo "Option $1 isn't recognized." >&2
      exit 3
      ;;
  esac
done

## Check that all parameters are set
flag_mp=0 # Flag for missing parameters
 
# Parallelization parameters
if [ -z "$nb_jobs" ]; then
  flag_mp=1 && echo -e "-nj/--nb_jobs is missing. \n\t$nb_jobs_def" >&2
fi
if [ -z "$nb_chunks" ]; then
  flag_mp=1 && echo -e "-nc/--nb_chunks is missing. \n\t$nb_chunks_def" >&2
fi
if [ -z "$nb_threads" ]; then
  flag_mp=1 && echo -e "-p/--nb_threads is missing. \n\t$nb_threads_def" >&2
fi

# Datasets processing parameters
if [ -z "$genome" ]; then
  flag_mp=1 && echo -e "-r/--genome is missing. \n\t$genome_def" >&2
fi
if [ -z "$first_chunk" -a -z "$last_chunk" ]; then
  first_chunk=1
  last_chunk=$(nanodisco chunk_info -r "$genome" | sed 's/.*: //')
else
  if [ -z "$first_chunk" ]; then
    flag_mp=1 && echo -e "-f/--first_chunk is missing. \n\t$first_chunk_def" >&2
  fi
  if [ -z "$last_chunk" ]; then
    flag_mp=1 && echo -e "-l/--last_chunk is missing. \n\t$last_chunk_def" >&2
  fi
fi
if [ -z "$path_input" ]; then
  flag_mp=1 && echo -e "-i/--path_input is missing. \n\t$path_input_def" >&2
fi
if [ -z "$path_output" ]; then
  flag_mp=1 && echo -e "-o/--path_output is missing. \n\t$path_output_def" >&2
fi
if [ -z "$wga_name" ]; then
  flag_mp=1 && echo -e "-w/--wga_name is missing. \n\t$wga_name_def" >&2
fi
if [ -z "$nat_name" ]; then
  flag_mp=1 && echo -e "-n/--nat_name is missing. \n\t$nat_name_def" >&2
fi

# Exit if missing parameters
if [ $flag_mp -eq 1 ]; then
  echo -e "\nPlease provide the parameters mentioned above and run again." >&2
  exit 4
fi

# Check if input data exist (nanodisco preprocessing)
missing_input_data="Please check -i/--path_input and/or run nanodisco preprocessing."
if [[ ! -f "${path_input%/}/$wga_name.fasta" ]]; then
  echo "${path_input%/}/$wga_name.fasta file doesn't exist. $missing_input_data" >&2
  exit 6
fi
if [[ ! -f "${path_input%/}/$wga_name.sorted.bam" ]]; then
  echo "${path_input%/}/$wga_name.sorted.bam file doesn't exist. $missing_input_data" >&2
  exit 6
fi
if [[ ! -f "${path_input%/}/$wga_name.sorted.bam.bai" ]]; then
  echo "${path_input%/}/$wga_name.sorted.bam.bai file doesn't exist. $missing_input_data" >&2
  exit 6
fi
if [[ ! -f "${path_input%/}/$nat_name.fasta" ]]; then
  echo "${path_input%/}/$nat_name.fasta file doesn't exist. $missing_input_data" >&2
  exit 6
fi
if [[ ! -f "${path_input%/}/$nat_name.sorted.bam" ]]; then
  echo "${path_input%/}/$nat_name.sorted.bam file doesn't exist. $missing_input_data" >&2
  exit 6
fi
if [[ ! -f "${path_input%/}/$nat_name.sorted.bam.bai" ]]; then
  echo "${path_input%/}/$nat_name.sorted.bam.bai file doesn't exist. $missing_input_data" >&2
  exit 6
fi

# Create a temporary function with requested parameters
eval "compute_difference() { chunks=\$1; Rscript --vanilla /home/nanodisco/code/difference_chunks.R -p $nb_threads -x $exec_type -f \${chunks%_*} -l \${chunks#*_} -i $path_input -o $path_output -w $wga_name -n $nat_name -r $genome -b $sig_norm -a $IQR_factor -z $normalized -e $min_coverage -j $map_type -k $min_read_length -t $path_ont_model; }"
export -f compute_difference # Export the function to be accessible by parallel

# Process genomic chunks; can be speed up by increasing nb_chunks but will impact memory usage
# ~6 min with 25 threads available (nb_jobs=5, nb_chunks=1 and nb_threads=5), chunk by chunk (each 5000 bp)
# Analysis logs are conserved in ./difference_log/1 (folder are named <first_chunk>_<last_chunk>)
if [ $nb_jobs -eq 1 ]; then
  compute_difference "${first_chunk}_${last_chunk}"
else
  for chunk in $(seq $first_chunk $nb_chunks $last_chunk); do
    chunk_to=$((chunk+$nb_chunks-1))
    chunk_to=$((chunk_to>last_chunk ? last_chunk : chunk_to))
    echo "${chunk}_${chunk_to}"
  done | parallel --no-notice --progress --results $path_output/difference_log -P $nb_jobs compute_difference {} 1> /dev/null
fi

