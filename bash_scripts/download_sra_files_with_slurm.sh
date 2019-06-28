#!/bin/bash

# Usage
# -----
# . download_sra_files_with_slurm.sh <accession_list_file> <output_folder>
# 
# where: 
#    - <accession_list_file> contains SRR accession numbers, one per line
#    - <output_folder> place to put the sra files. Create if not exists. 

accession_file=$1
accession_file=$(realpath $accession_file)
accession_dir=$(dirname $accession_file)
output_dir="$accession_dir/$2"

printf "> The accession file we are using now us '$accession_file'\n"
printf "> The directory of the accession file is '$accession_dir'\n"
printf "> The output directory will be '$output_dir'\n"

mkdir -p $output_dir

cwd=$(pwd)
cd $accession_dir
while read line; do
  echo "fetching $line"
  output_file="$output_dir/$line.sra"
  sbatch_filename="$output_dir/$line.sbatch"
  if [ ! -f "$output_file" ]; then 
    printf $"#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --time=120:00\n#SBATCH -c 1\n#SBATCH -J $accession_file\nprefetch $line --output-file $output_file" > "$sbatch_filename"
    sbatch "$sbatch_filename"

  else
    echo "Already have a file called $output_file. Skipping $line"
  fi
done < "$accession_file"

cd $cwd






