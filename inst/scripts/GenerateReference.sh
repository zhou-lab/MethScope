#!/bin/bash

### change these parameters accordingly
path="/home/fuh1/tmp/samplesheet/20220923_Loyfer_WGBS.tsv" # path for the cell type label tsv
cg_path="/mnt/isilon/zhou_lab/projects/20230727_all_public_WGBS/hg38/2023_Loyfer.cg" # path for the .cg data
COLUMN_INDEX=7 # column index for the cell type label
SAMPLE_INDEX=1 # column index for the sample id that match with .cg file
min_CG=70 # the minimum CpG # a pattern must have (Default: 50)

if [[ ! -f "$path" ]]; then
    echo "Error: File $path does not exist."
    exit 1
fi

if [[ ! -f "$cg_path" ]]; then
    echo "Error: File $cg_path does not exist."
    exit 1
fi

mergeCG2 () 
{ 
    fdr=$1;
    fn=$2;
    : > $fn;
    find $fdr -name '*.cac' -o -name '*.cg' | LC_ALL=C sort | while read f; do
        cat $f >> $fn;
        yame index -1 $(basename $(basename $f .cg) .ca) $fn;
    done
}

mkdir -p ./tmp_folder/
awk -F'\t' -v col="$COLUMN_INDEX" '{print $col}' $path | sort -u | \
while read value; do
awk -F'\t' -v val="$value" -v col="$COLUMN_INDEX" '$col == val' $path | cut -f $SAMPLE_INDEX | yame subset -l - $cg_path | yame rowop -o musum - ./tmp_folder/${value}.cg
done
mergeCG2 ./tmp_folder/ ./tmp.cg
rm -rf ./tmp_folder/
yame rowop -o binstring ./tmp.cg > ./tmp
rm -f ./tmp.cg
rm -f ./tmp.cg.idx
Rscript ./GenerateReference.R ./tmp $min_CG
yame pack -f s ./patterns.txt ./patterns.cm