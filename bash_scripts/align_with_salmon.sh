#!bin/bash

# get all 3 prime files in a space separated array
dir=$1

ones=$(find SRR*pass*1*.fq)
twos=$(find SRR*pass*2*.fq)

# convert to array
ones=($ones)
twos=($twos)

ref_dir="/media/ncw135/DATA/PublicRepositories/geodatasets/reference_genome/transcriptome/human_transcripome_index"
i=0

#echo $ones
for f in "${ones[@]}";
do
  ((i++))
  one=${ones[$i]}
  two=${twos[$i]}
  accession_number=$(echo "$one" | grep -Po '(SRR\d*)')
  out="$dir/quants/$accession_number"

  salmon quant -l A -i "$ref_dir" --validateMappings -o $out -1 "$one" -2 "$two"

done



#ones=$(find SRR*pass*1*.fq)
#twos=$(find SRR*pass*2*.fq)
#
#for i in $ones;
#do
#  fle=$(basename $i)
#  out=$(echo $fle | grep -Po '(SRR\d*)')"_quant"
#  echo "$out"
#done
#
#ref_dir="/media/ncw135/DATA/PublicRepositories/geodatasets/reference_genome/transcriptome/human_transcripome_index"
#
##echo $ones
##echo $out
##parallel --link -j 3 salmon quant -l A -i $ref_dir --validateMappings -o {3} -1 {1} -2 {2} ::: $ones ::: $twos ::: $out
#salmon quant -l A -i $ref_dir --validateMappings -o {3} -1 {1} -2 {2} ::: $ones ::: $twos ::: $out

