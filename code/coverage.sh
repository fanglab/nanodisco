#!/bin/bash

# Example command using provide dataset with get_data_bacteria.
# nanodisco coverage -b /home/nanodisco/analysis/EC_NAT.sorted.bam -r /home/nanodisco/reference/Ecoli_K12_MG1655_ATCC47076.fasta -o /home/nanodisco/analysis/

###
# Compute average coverage in each contig from the metagenome
###

## Describe parameters
path_mapping_def="Specify the path of mapping data (.sorted.bam)."
metagenome_def="Specify the path of reference metagenome (.fasta)."
path_output_def="Specify the path of output directory (.sorted.bam suffix replaced by .cov)."

# Read options
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -b|--path_mapping)
      if [[ -f $2 ]]; then
        if [[ $2 =~ ^.*\.sorted.bam$ ]]; then
          path_mapping="$2"
          shift # pass argument
          shift # pass value
        else
          echo "-b/--path_mapping not matching expected file name (<name_sample>.sorted.bam). Please check -b/--path_mapping and/or run nanodisco preprocessing." >&2
          echo -e "\t"$path_mapping_def
          exit 5
        fi
      else
        echo "-b/--path_mapping file doesn't exist. Please check -b/--path_mapping and/or run nanodisco preprocessing." >&2
        echo -e "\t"$path_mapping_def
        exit 5
      fi
      ;;
    -r|--metagenome)
      if [[ -f $2 ]]; then
        metagenome="$2"
        shift # pass argument
        shift # pass value
      else
        echo "-r/--metagenome file doesn't exist." >&2
        echo -e "\t"$metagenome_def >&2
        exit 5
      fi
      ;;
    -o|--path_output)
      path_output="${2%/}" # Strip potential slash
      shift # pass argument
      shift # pass value
      ;;
    -h|--help) # Print help
      echo "Usage: nanodisco coverage -b <path_mapping> -r <path_metagenome> -o <path_output>" >&2
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
 
# Datasets processing parameters
if [ -z "$path_mapping" ]; then
  flag_mp=1 && echo -e "-b/--path_mapping is missing. \n\t$path_mapping_def" >&2
fi
if [ -z "$metagenome" ]; then
  flag_mp=1 && echo -e "-r/--metagenome is missing. \n\t$metagenome_def" >&2
fi
if [ -z "$path_output" ]; then
  flag_mp=1 && echo -e "-o/--path_output is missing. \n\t$path_output_def" >&2
fi

# Exit if missing parameters
if [ $flag_mp -eq 1 ]; then
  echo -e "\nPlease provide the parameters mentioned above and run again." >&2
  exit 4
fi

# Define output path and file name
input_file_name=$(basename $path_mapping)
path_coverage="$path_output/${input_file_name/.sorted.bam/.cov}"
if [[ -d $path_output ]]; then
  if [[ -f $path_coverage ]]; then
    echo "$path_coverage already exist and will be overwritten." >&2
  fi
else
  mkdir -p $path_output
fi

bedtools genomecov -ibam $path_mapping -g $metagenome > $path_coverage
