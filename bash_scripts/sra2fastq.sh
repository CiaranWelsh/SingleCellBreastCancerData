#!/bin/bash

# takes a directory containing sra files as input and puts fastq files 
#  into the same directory. Uses fastq-dump

dir=$(realpath "$1")
echo "converting all sra files in $dir to fastq"
sra_files=$(find $dir -name "*.sra")

for f in $sra_files
do
  output_directory="$(dirname $f)"
  output_filename="$(basename $f)"
  output_filename="${output_filename%.*}"
  out="$output_directory/$output_filename.fastq"
  echo "> Converting $output_filename to fastq format. Output is $out"
  fastq-dump $f --outdir $out --split-files --readids -origfmt --clip --read-filter pass 
done

