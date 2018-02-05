#!/bin/bash

if [ $# != 3 ]; then
    echo "Usage: ./plot_cnv.sh sample.cnv.txt sample.normalized.coverage.bed sample.models.bed"
    exit 1
fi

sample=`echo "$1" | awk '{ n=split($1, arr, "/"); print arr[n] }' | cut -d '.' -f 1`
cnv=$1
cov=$2
models=$3

SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cut -f -4 $cov >tmp.$sample.cov.txt

awk '{
    if ($4 >= 0) {
        n++;
        window[n] = $1 "\t" $2 "\t" $3;
        start[n] = $2;
        end[n] = $3;
    }
} END {
    for (i = 4; i <= n-3; i++) {
        print window[i] "\t" start[i-3] "\t" end[i+3]
    }
}' $models \
| sort -k1,1 -k2,2n \
>tmp.$sample.prev-and-next.windows.bed

cut -f -7 $cnv \
| bedtools map -a - -b tmp.$sample.prev-and-next.windows.bed -c 4,5 -o first,last \
| awk 'BEGIN { OFS="\t"; } {
    if ($8 != "." && $9 != ".") { $2 = $8; $3 = $9 }
    $1 = $1;
    printf "%s", $1;
    for (i = 2; i <= 7; i++) printf "\t%s", $i;
    printf "\n";
}' \
| bedtools intersect -a - -b tmp.$sample.cov.txt -wb \
| sort -k1,1 -k2,2n \
| bedtools map -a - -b $models -c 4,9,10 -o first,first,first \
| awk '{
    split($4, arr1, ":"); split(arr1[2], arr2, "-");
    x1 = 0; if ($2 >= arr2[1]-1 && $3 <  arr2[2]) x1 = 1;
    x2 = 0; if ($2 >= arr2[1]-1 && $3 <= arr2[2]) x2 = 1;
    print $5 "\t" $4 "\t" $6 "\t" $7 "\t" $1 ":" $2 "-" $3 "\t" x1 "\t" x2 "\t" $11 "\t" $12 "\t" $13 "\t" $14 }' >tmp.plot_cnv.$sample.txt

mkdir -p clamms_cnv_plots/$sample
Rscript $SCRIPT_DIR/plot_cnv.R tmp.plot_cnv.$sample.txt

rm tmp.$sample.cov.txt 
rm tmp.$sample.prev-and-next.windows.bed
rm tmp.plot_cnv.$sample.txt

