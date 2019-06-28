#!/bin/bash

while read line; do
  echo "trying to download $line"
  out="$(dirname $(realpath $1))/$line.fastq"
  if [ ! -f "$out" ]; then
    echo "$line not found. Downloading"
    fastq-dump $line --outdir $out -v --split-files --readids -origfmt --clip --read-filter pass &
  else
    echo "$line was found already. Not Downloading"
  fi
done < "$1"


