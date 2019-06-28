#!bin/bash

# Usage
# . trim_fastq_adaptors.sh directory output threads
# Where
#  - directory: a directory containing fastq files. Program recursively finds all fastq in subfolders
#  - output: 	folder to store the output
#  - threads: 	number of threads to use in trim_galore
source activate py36
dir="$(realpath $1)"
cd "$dir"
fastq_folders=$(find -name "SRR*" -type d)
output_folder="$dir/trimmed"

trim_galore="module add Trim_Galore/0.6.1-foss-2018b-Python-3.6.6"
cutadapt="module add cutadapt/1.18-foss-2018b-Python-3.6.6"
fastqc="module add FastQC/0.11.8-Java-1.8.0_144"

for i in $fastq_folders; 
do
  cd $i
  #printf "> i is: $i\n"
  one=$(find -name "SRR*_*1*.fastq")
  two=$(find -name "SRR*_*2*.fastq")
  #printf "> one is: $one\n"
  submit_file=$(echo $one | grep -Po '(SRR\d*)')".submit"
  submit_file=$(realpath "$submit_file")
  #printf "> submit file is: $submit_file\n"
  printf "#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --time=120:00\n#SBATCH -c 1\n$trim_galore\n$cutadapt\n$fastqc\ntrim_galore -j 1 --paired --fastqc -o $output_folder $one $two" > $submit_file
  sbatch $submit_file
  rm $submit_file
  cd $dir
done





