#!/bin/bash
source activate py36

dir=$(realpath $1)
temp_dir="\$TMPDIR"

ones=$(find SRR*pass*1*.fq)
twos=$(find SRR*pass*2*.fq)

# convert to array
ones=($ones)
twos=($twos)

ref_dir="/mnt/storage/nobackup/b3053674/scRNA-seq/reference_genome/transcriptome/human_transcriptome_index"
i=0

#echo $ones
for f in "${ones[@]}";
do
  ((i++))
  one=${ones[$i]}
  two=${twos[$i]}
  accession_number=$(echo "$one" | grep -Po '(SRR\d*)')
  out="$dir/quants/$accession_number"

  temp_one="$temp_dir/$one"
  temp_two="$temp_dir/$two"
  temp_ref="$temp_dir/human_transcriptome_index"
  temp_out="$temp_dir/accession_number"

  slurm="$dir/$accession_number.sbatch"

  cp_string="cp $one $temp_one\ncp $two $temp_two\ncp -R $ref_dir $temp_ref"
  sbatch_commands="#SBATCH --nodes=1\n#SBATCH --time=120:00\n#SBATCH -c 4"
  salmon_command="salmon quant -p 4 -l A -i $temp_ref --output $temp_out ba--validateMappings --mates1 '$temp_one' --mates2 '$temp_two' "
  cp_results_command="cp -R $temp_out $out"
  printf "#!/bin/bash\n$sbatch_commands\n$cp_string\n$salmon_command\n$cp_results_command" > $slurm
#  sbatch $slurm
##  #rm $slurm
done




