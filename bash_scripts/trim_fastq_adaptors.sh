#!bin/bash

# Usage
# . trim_fastq_adaptors.sh directory output threads
# Where
#  - directory: a directory containing fastq files. Program recursively finds all fastq in subfolders
#  - output: 	folder to store the output
#  - threads: 	number of threads to use in trim_galore

dir=$1
fastq_folders=$(find $dir -name "*.fastq")
declare -a fastq_files
str=""


for i in $fastq_folders
do
  full_path=$(realpath $i)
  if test -f "$full_path"; then
    fastq_files+=("$full_path")
    str="$str $full_path"
  fi
done

mkdir -p $2
echo $str

# -o output directory
# -j number of cores
trim_galore -o $2 -j $3 $str

