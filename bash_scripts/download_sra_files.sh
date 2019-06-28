#!bin/bash


accession_file=$1
echo $accession_file
#realpath $accession_file
accession_file=$(realpath $accession_file)
accession_dir=$(dirname $accession_file)

echo $accession_dir
echo "doing $accession_file"
cwd=$(pwd)
cd $accession_dir
while read line; do
  echo "fetching $line"
  output_file="$accession_dir/$line.sra"
  if [ ! -f "$output_file" ]; then 
    prefetch $line --output-file "$output_file"
    echo "output to $accession_dir"
  else
    echo "Already have a file called $output_file. Skipping $line"
  fi
done < "$accession_file"

cd $cwd






