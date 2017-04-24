#!/bin/bash


for line in $(ls ../data/*.mtx); do
    echo "===Matrix: $line==="
    ./pgmres $line
    printf "\n\n"
done

# ./pgmres ../data/bcspwr01.mtx
# ./pgmres ../data/bcspwr03.mtx
# ./pgmres ../data/bcspwr06.mtx
# ./pgmres ../data/cage4.mtx
# ./pgmres ../data/cryg2500.mtx
# ./pgmres ../data/mhd4800a.mtx
# ./pgmres ../data/cryg10000.mtx
