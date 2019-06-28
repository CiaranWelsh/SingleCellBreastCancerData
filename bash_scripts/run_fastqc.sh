#!bin/bash


# Find all folders with the fastq extension and passes the list to fastqc
# Arguments
# ---------
# - $1: directory to search for fastq files. Also looks in subdirectories
# - $2: directory for output files
# - $3: argument to the -t option in fastqc. Number of threads to use.  

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

fastqc -o $2 --extract -t $3 $str

