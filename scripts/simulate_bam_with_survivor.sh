#!/bin/bash

name="$1"
survivor_path="../$2" # path to SURVIVOR/Debug/SURVIVOR
ref="../$3" # path to the reference fasta to have SVs added to it
param_file="../$4"
error_profile="../$5"

mkdir $name
cd $name
echo "Creating fasta."
$survivor_path simSV $ref $param_file 0.0 0 $name
$survivor_path simreads $name.fasta $error_profile 30 $name.reads
minimap2 -a $ref $name.reads > $name.sam
samtools sort -o $name.bam $name.sam
samtools index $name.bam
cd ..