#!/bin/bash

if [ -z "$CLAMMS_DIR" ]; then
    echo "Environment variable CLAMMS_DIR must be set."
    exit 1
fi  

if [[ $# != 5 ]]; then
    echo "Usage: ./annotate_windows targets.bed genome.fa mappability.bed insert_size special_regions.bed >windows.bed"
    exit 1
fi

awk -f $CLAMMS_DIR/split_targets_into_windows.awk $1 \
| awk -v INSERT_SIZE=$4 '{
    window_len = $3 - $2;
    if (INSERT_SIZE > window_len) window_len = INSERT_SIZE;
    left_bases = int(window_len / 2);
    right_bases = left_bases;
    if (window_len % 2 == 1) right_bases++;
    mid = int(($3+$2)/2);
    printf "%s\t%d\t%d\t%s\n",
            $1, mid - left_bases, mid + right_bases, $4 }' \
| bedtools nuc -fi $2 -bed - \
| tail -n +2 | cut -f 1-4,6 \
| awk '{
    split($4, arr1, ":");
    split(arr1[2], arr2, "-");
    printf "%s\t%d\t%d\t%s\t%.0f\t%.3f\n",
            arr1[1], arr2[1], arr2[2], $4, $5*200, $5  }' \
| $CLAMMS_DIR/calc_window_mappability.py $3 \
| bedtools map -a - -b $5 -c 4 -o min | awk '
    BEGIN { OFS="\t"; } $8 == "." { $8 = 3; } { print; }
'

