#!/bin/bash

while read line; do
  echo "trying to download $line"
  out="$(dirname $(realpath $1))/published_sam_files/$line.sam"
  if [ ! -f "$out" ]; then
    echo "$line not found. Downloading"
    sam-dump $line --output-file $out -v 
  else
    echo "$line was found already. Not Downloading"
  fi
done < $1


