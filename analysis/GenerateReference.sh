#!/bin/bash
#SBATCH --job-name=generate_ref
#SBATCH --output=logs/generate_ref_%j.out
#SBATCH --error=logs/generate_ref_%j.err
#SBATCH --mem=420G
#SBATCH -c 24
#SBATCH -t 30:00:00

### change these parameters accordingly
#path="/home/fuh1/cell_type/20250425_eckerhuman.tsv" # path for the cell type label tsv
#cg_path="/home/fuh1/zhou_lab/projects/20220324_SingleCellMeth/hg38/2025_Zhou.cg" # path for the .cg data

path="$1"     
cg_path="$2" 

### Note: row names might be the first index, so you should add 1 to the index if that's the case 
### use cut -f to check that's the right index column 
COLUMN_INDEX=30 # column index for the cell type label
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

export path cg_path COLUMN_INDEX SAMPLE_INDEX

mkdir -p ./tmp_folder/
# awk -F'\t' -v col="$COLUMN_INDEX" '{print $col}' $path | sort -u | \
# while read value; do
# awk -F'\t' -v val="$value" -v col="$COLUMN_INDEX" '$col == val' $path | cut -f $SAMPLE_INDEX | yame subset -l - $cg_path | yame rowop -o musum - ./tmp_folder/${value}.cg
# done

### use parallel
#cut -f"$COLUMN_INDEX" "$path" skip header line
awk -F'\t' -v col="$COLUMN_INDEX" 'NR>1 {print $col}' "$path" | sort -u | \
parallel --jobs 12 --env path --env cg_path --env COLUMN_INDEX --env SAMPLE_INDEX '
val={}
safe_val=$(echo "$val" | sed "s/[^a-zA-Z0-9_.-]/_/g")
awk -F"\t" -v val="$val" -v col="$COLUMN_INDEX" '\''$col == val'\'' "$path" | \
cut -f "$SAMPLE_INDEX" | \
yame subset -l - "$cg_path" | \
yame rowop -o musum - ./tmp_folder/"$safe_val".cg
'

mergeCG2 ./tmp_folder/ ./tmp.cg
rm -rf ./tmp_folder/
rm -f ./patterns.cm
yame rowop -o binstring ./tmp.cg > ./tmp
rm -f ./tmp.cg
rm -f ./tmp.cg.idx
Rscript ./GenerateReference.R ./tmp $min_CG
yame pack -f s ./patterns.txt ./patterns.cm
#Rscript ./FilterArchetype.R ./tmp $min_CG
#yame pack -f s ./patterns_updated.txt ./patterns_updated.cm
