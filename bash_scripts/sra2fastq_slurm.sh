#!/bin/bash

# Convert sra files to fastq using slurm
# 
# 
# Usage
# -----
# . sra2fastq_slurm.sh <sra_directory> 
# 
# where: 
#    - <sra_directory> is any directory that contains sra files. 

dir=$(realpath "$1")
printf "converting all sra files in $dir to fastq"
sra_files=$(find $dir -name "*.sra")

for f in $sra_files
do
  output_directory="$(dirname $f)"
  output_filename="$(basename $f)"
  output_filename="${output_filename%.*}"
  out="$output_directory/$output_filename"
  sbatch_script="$output_directory/$output_filename.sbatch"
  printf "\n>Converting '$f' to fastq format.\n>Batch script is '$sbatch_script'.\n>Output is '$out'"
  printf "#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --time=120:00\n#SBATH -c 1 \n#SBATCH \n\nmodule load SAMtools/1.9-foss-2018b\nfastq-dump $f --outdir $out --split-files --readids -origfmt --clip --read-filter pass" > $sbatch_script
  sbatch $sbatch_script
  rm $sbatch_script
done

