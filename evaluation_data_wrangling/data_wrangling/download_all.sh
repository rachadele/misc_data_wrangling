#!/bin/bash


file=$1

cat $file | while read line; do

study=$(echo $line | cut -d ' ' -f1)

echo $study

gemma-cli-sc getSingleCellDataMatrix -e $study --format mex --scale-type count --use-ensembl-ids -o "/space/scratch/gemma-single-cell-data-ensembl-id/$study"
curl -u "$GEMMA_USERNAME:$GEMMA_PASSWORD" "https://gemma.msl.ubc.ca/rest/v2/datasets/$study/samples" | jq > "/space/scratch/gemma-single-cell-data-ensembl-id/$study/${study}_meta.json" 

done



