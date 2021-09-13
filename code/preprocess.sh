#!/bin/bash

# Example command using provide dataset with get_data_bacteria.
# nanodisco preprocess -p 40 -f /home/nanodisco/dataset/EC_WGA -s EC_WGA -o /home/nanodisco/analysis/test -r /home/nanodisco/reference/Ecoli_K12_MG1655_ATCC47076.fasta
# nanodisco preprocess -p 40 -f /home/nanodisco/dataset/EC_NAT -s EC_NAT -o /home/nanodisco/analysis/test -r /home/nanodisco/reference/Ecoli_K12_MG1655_ATCC47076.fasta

###
# Pre-processing a new dataset
###
# 1. Fasta extraction
# 2. Reference genome preparation
# 3. Read mapping

function print_message (){
  echo "["$(date +"%Y-%m-%d %T")"] $1."
}

function check_error (){
  if [[ "$1" -ne 0 ]]; then
    cat "$2_log"
  fi
  rm "$2_log"
}

## Describe parameters
nb_threads_def="Specify how many threads to use."
path_fast5_def="Specify name for the sample (e.g. EC_NAT for a native E. coli sample)"
name_sample_def="Specify the name of the sample."
path_output_def="Specify the path of output results (.fasta, sorted.bam, and .sorted.bam.bai)."
path_reference_genome_def="Specify the path of reference genome (.fasta)."

# Handle missing parameters
if [[ "$#" -eq 0 ]]; then
  echo "Usage: nanodisco preprocess -p <nb_threads> -f <path_fast5> -s <sample_name> -o <path_output> -r <path_reference_genome>" >&2
  echo -e "\tAdditional information can be found in our GitHub repository." >&2
fi

# Read options
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
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
    -f|--path_fast5)
      if [[ -d $2 ]]; then
        path_fast5="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-f/--path_fast5 directory doesn't exist." >&2
        echo -e "\t"$path_fast5_def
        exit 5
      fi
      ;;
    -s|--name_sample)
      name_sample="$2"
      shift # pass argument
      shift # pass value
      ;;
    -o|--path_output)
      if [[ -d $2 ]]; then
        # Directory already exist
        # echo "-o/--path_output directory already exist. Some data could be overwrite." >&2
        :
      fi
      path_output="${2%/}" # Strip potential slash
      shift # pass argument
      shift # pass value
      ;;
    -r|--path_reference_genome)
      if [[ -f $2 ]]; then
        path_reference_genome="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-r/--path_reference_genome file doesn't exist." >&2
        echo -e "\t"$path_reference_genome_def >&2
        exit 5
      fi
      ;;
    -c|--nb_chunks)
      if [[ $2 =~ ^[\-0-9]+$ ]] && (( $2 > 0)); then
        nb_chunks="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-c/--nb_chunks positive integer required." >&2
        echo -e "\t"$path_reference_genome_def >&2
        exit 5
      fi
      ;;
    --basecall_version)
      basecall_version="$2"
      shift # pass argument
      shift # pass value
      ;;
    -h|--help) # Print help
      echo "Usage: nanodisco preprocess -p <nb_threads> -f <path_fast5> -s <sample_name> -o <path_output> -r <path_reference_genome>" >&2
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

if [ -z "$nb_threads" ]; then
  flag_mp=1 && echo -e "-p/--nb_threads is missing. \n\t$nb_threads_def" >&2
fi
if [ -z "$path_fast5" ]; then
  flag_mp=1 && echo -e "-f/--path_fast5 is missing. \n\t$path_fast5_def" >&2
fi
if [ -z "$name_sample" ]; then
  flag_mp=1 && echo -e "-s/--name_sample is missing. \n\t$name_sample_def" >&2
fi
if [ -z "$path_output" ]; then
  flag_mp=1 && echo -e "-o/--path_output is missing. \n\t$path_output_def" >&2
fi
if [ -z "$path_reference_genome" ]; then
  flag_mp=1 && echo -e "-r/--path_reference_genome is missing. \n\t$path_reference_genome_def" >&2
fi
if [ -z "$basecall_version" ]; then
  basecall_version="default"
fi
if [ -z "$nb_chunks" ]; then
  nb_chunks=40
fi


# Exit if missing parameters
if [ $flag_mp -eq 1 ]; then
  echo -e "\nPlease provide the parameters mentioned above and run again." >&2
  exit 4
fi

# Setup final output path for the dataset
mkdir -p $path_output
path_output_sample="$path_output/$name_sample"

## Extract fasta files from fast5, necessary as nanopolish create internal read<->fast5 index
/home/nanodisco/code/extract.R -i $path_fast5 -o $path_output -b $name_sample -p $nb_threads -s fa -c $nb_chunks --basecall_version $basecall_version

# Check if some reads have been extracted
if [[ -f $path_output_sample".fasta" ]]; then
  if [[ $(wc -l < $path_output_sample".fasta") -eq 0 ]]; then
    echo "No reads were extracted. Please check that -f/--path_fast5 is correct."
    exit 6
  fi
else
  echo "Unexpected error during read extraction process."
  exit 7
fi

## Prepare reference
# if [[ ! -f $path_reference_genome".amb" || ! -f $path_reference_genome".ann" || ! -f $path_reference_genome".bwt" || ! -f $path_reference_genome".pac" || ! -f $path_reference_genome".sa" ]]; then
#   print_message "Prepare bwa index"
#   bwa index $path_reference_genome 2> $path_output_sample"_log" # Needed for mapping reads
#   check_error "$?" $path_output_sample # Hide bwa index stderr by default and show only if an error is identified
# fi
if [[ ! -f $path_reference_genome".mmi" ]]; then
  print_message "Prepare minimap2 index"
  minimap2 -x map-ont -d ${path_reference_genome/.fasta/.mmi} $path_reference_genome > $path_output_sample"_log" 2>&1 # Needed for mapping reads
  check_error "$?" $path_output_sample # Hide bwa index stderr by default and show only if an error is identified
fi

if [[ ! -f $path_reference_genome".fai" ]]; then
  print_message "Prepare samtools index"
  samtools faidx $path_reference_genome 2> $path_output_sample"_log" # Needed for indexing reads
  check_error "$?" $path_output_sample # Hide bwa index stderr by default and show only if an error is identified
fi

## Map reads on reference
print_message "Map reads"
# bwa mem -t $nb_threads -x ont2d $path_reference_genome $path_output_sample".fasta" 2> $path_output_sample"_log" | samtools view -b - | samtools sort -T $path_output_sample"_tmp" > $path_output_sample".sorted.bam"
minimap2 -ax map-ont -t $nb_threads ${path_reference_genome/.fasta/.mmi} $path_output_sample".fasta" 2> $path_output_sample"_log" | samtools sort -o $path_output_sample".sorted.bam" -T $path_output_sample"_tmp" -
check_error "$?" $path_output_sample # Hide bwa mem stderr by default and show only if an error is identified

## Index alignment
print_message "Index alignment"
samtools index $path_output_sample".sorted.bam"

print_message "Done"
