#!/bin/bash

rm logs/*

for line in cage4.mtx cage5.mtx cage7.mtx cage8.mtx cage9.mtx cage10.mtx; do
    echo "===Matrix: $line==="

    for i in {1..10}; do
        ./dense_seq/dense_seq ../data/$line 2>&1 | tee -a logs/dense_seq_${line%.mtx}.log
    done

    printf "\n"

    for i in {1..10}; do
        ./dense_omp/dense_omp ../data/$line 2>&1 | tee -a logs/dense_omp_${line%.mtx}.log
    done

    printf "\n"

    for i in {1..10}; do
        ./sparse_omp/sparse_omp ../data/$line 2>&1 | tee -a logs/sparse_omp_${line%.mtx}.log
    done

    printf "\n"
done

for line in bcspwr01.mtx  bcspwr03.mtx  bcspwr06.mtx  bcspwr10.mtx; do
    echo "===Matrix: $line==="

    for i in {1..10}; do
        ./dense_seq/dense_seq ../data/bcspwr/$line 2>&1 | tee -a logs/dense_seq_${line%.mtx}.log
    done

    printf "\n"

    for i in {1..10}; do
        ./dense_omp/dense_omp ../data/bcspwr/$line 2>&1 | tee -a logs/dense_omp_${line%.mtx}.log
    done

    printf "\n"

    for i in {1..10}; do
        ./sparse_omp/sparse_omp ../data/bcspwr/$line 2>&1 | tee -a logs/sparse_omp_${line%.mtx}.log
    done

    printf "\n"
done
