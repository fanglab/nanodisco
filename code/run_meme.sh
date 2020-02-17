#!/bin/bash

fasta=$1
output_dir=$2
mod=$3
nbMotif=$4
minMotifSize=$5
maxMotifSize=$6
nb_threads=$7

meme $fasta -dna -oc $output_dir -mod $mod -nmotifs $nbMotif -minw $minMotifSize -maxw $maxMotifSize -p $nb_threads -maxsize 1000000
