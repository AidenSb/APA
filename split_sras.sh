#!/bin/bash
for f in *; 
do 
fastq-dump --gzip --split-files $f; 
done
