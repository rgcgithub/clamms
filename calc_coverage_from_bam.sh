#!/bin/bash

if [[ $# != 3 ]]; then
    echo "Usage: ./calc_coverage windows.bed sample.bam min_map_qscore"
    exit 1
fi

samtools bedcov -Q $3 $1 $2 \
| awk '{ printf "%s\t%d\t%d\t%.6g\n", $1, $2, $3, $NF/($3-$2); }'

