#!/bin/bash

export HDF5_PLUGIN_PATH=/nanopolish-0.13.3/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin

function check_ref_fasta (){
  ref_fasta=$1

  if [ ! -f $ref_fasta ]; then
    echo "Problem with reference: "$ref_fasta
    exit 2
  fi
}

function check_bwa_index (){
  ref=$1

  if [ ! -e "$ref.amb" ] || [ ! -e "$ref.ann" ] || [ ! -e "$ref.pac" ] || [ ! -e "$ref.bwt" ] || [ ! -e "$ref.sa" ]; then
    bwa index $ref
  fi
  if [[ ! -f $path_reference_genome".mmi" ]]; then
    minimap2 -x map-ont -d ${ref/.fasta/.mmi} $ref
  fi
}

function correct_events (){
  subset_nb_threads=$1
  in_fasta=$2
  in_bam=$3
  tmp_genome=$4
  out_eventalign=$5
  out_error_file=$6

  if [ -s $in_fasta ]; then
    nanopolish index $in_fasta 2> $out_error_file
    nanopolish eventalign -t $subset_nb_threads --scale-events -n -r $in_fasta -b $in_bam -g $tmp_genome > $out_eventalign 2>> $out_error_file
  else
    echo "contig" > $out_eventalign # Create empty file to handle by R
  fi
}

genome=$1
fasta=$2
nb_threads=$3
chunk=$4
strand_bias=$5
map_type=$6

if [ "$strand_bias" = "revc" ]; then
  revcgenome=${genome/fasta/rev_comp.fasta}
  st_rev_flag="-F"
elif [ "$strand_bias" = "ori" ]; then
  revcgenome=$genome
  st_rev_flag="-f"
else
  echo "Problem with strand bias correction type." # TODO Handle by R?
  exit 1
fi

check_ref_fasta $genome
check_ref_fasta $revcgenome

check_bwa_index $genome
check_bwa_index $revcgenome

# 16 mapped on rev
# 2048 not representative (not primary) supplementary alignment
# XA tag for alternative alignment, might not be 
# SA tag for supplementary alignment (on primary & other SA), chimeric alignment (linear aln withut large overlaps) 
# 2048 & 16 may remove real data, de novo assembly could resolve somes
# grep not elegant but seems only solution to remove multimapped reads.
# -v 1 to remove unwanted output on std
cmd_bwa="bwa mem -x ont2d -t $nb_threads"
cmd_mm="minimap2 -ax map-ont -t $nb_threads"
cmd_mapping=$cmd_mm
if [ "$map_type" == "all" ]; then
  $cmd_mapping $genome ${fasta/fasta/fwd.fasta} | samtools view -bF 16 - | samtools sort -T "fwd.$chunk" > ${fasta/fasta/fwd.sorted.bam}
  $cmd_mapping $revcgenome ${fasta/fasta/rev.fasta} | samtools view -b"$st_rev_flag" 16 - | samtools sort -T "rev.$chunk" > ${fasta/fasta/rev.sorted.bam}
elif [ "$map_type" == "noMulti" ]; then
  $cmd_mapping $genome ${fasta/fasta/fwd.fasta} | grep -v 'XA:Z:' | samtools view -bF 16 - | samtools sort -T "fwd.$chunk" > ${fasta/fasta/fwd.sorted.bam}
  $cmd_mapping $revcgenome ${fasta/fasta/rev.fasta} | grep -v 'XA:Z:' | samtools view -b"$st_rev_flag" 16 - | samtools sort -T "rev.$chunk" > ${fasta/fasta/rev.sorted.bam}
elif [ "$map_type" == "noAddSupp" ]; then
  $cmd_mapping $genome ${fasta/fasta/fwd.fasta} | samtools view -bF 2048 -F 16 - | samtools sort -T "fwd.$chunk" > ${fasta/fasta/fwd.sorted.bam}
  $cmd_mapping $revcgenome ${fasta/fasta/rev.fasta} | samtools view -bF 2048 "$st_rev_flag" 16 - | samtools sort -T "rev.$chunk" > ${fasta/fasta/rev.sorted.bam}
elif [ "$map_type" == "noSupp" ]; then
  $cmd_mapping $genome ${fasta/fasta/fwd.fasta} | grep -v 'SA:Z:' | samtools view -bF 16 - | samtools sort -T "fwd.$chunk" > ${fasta/fasta/fwd.sorted.bam}
  $cmd_mapping $revcgenome ${fasta/fasta/rev.fasta} | grep -v 'SA:Z:' | samtools view -b"$st_rev_flag" 16 - | samtools sort -T "rev.$chunk" > ${fasta/fasta/rev.sorted.bam}
elif [ "$map_type" == "noMultiAddSupp" ]; then
  $cmd_mapping $genome ${fasta/fasta/fwd.fasta} | grep -v 'XA:Z:' | samtools view -bF 2048 -F 16 - | samtools sort -T "fwd.$chunk" > ${fasta/fasta/fwd.sorted.bam}
  $cmd_mapping $revcgenome ${fasta/fasta/rev.fasta} | grep -v 'XA:Z:' | samtools view -bF 2048 "$st_rev_flag" 16 - | samtools sort -T "rev.$chunk" > ${fasta/fasta/rev.sorted.bam}
elif [ "$map_type" == "none" ]; then
  $cmd_mapping $genome ${fasta/fasta/fwd.fasta} | grep -v 'XA:Z:\|SA:Z:' | samtools view -bF 16 - | samtools sort -T "fwd.$chunk" > ${fasta/fasta/fwd.sorted.bam}
  $cmd_mapping $revcgenome ${fasta/fasta/rev.fasta} | grep -v 'XA:Z:\|SA:Z:' | samtools view -b"$st_rev_flag" 16 - | samtools sort -T "rev.$chunk" > ${fasta/fasta/rev.sorted.bam}
else
  echo "Problem with mapping type." # TODO Handle by R?
  exit 2
fi

wait
samtools index ${fasta/fasta/fwd.sorted.bam}
samtools index ${fasta/fasta/rev.sorted.bam}

correct_events $((nb_threads/2)) ${fasta/fasta/fwd.fasta} ${fasta/fasta/fwd.sorted.bam} $genome ${fasta/fasta/fwd.eventalign} ${fasta/fasta/fwd.eventalign.err}
correct_events $((nb_threads/2)) ${fasta/fasta/rev.fasta} ${fasta/fasta/rev.sorted.bam} $revcgenome ${fasta/fasta/rev.eventalign} ${fasta/fasta/rev.eventalign.err}
wait
