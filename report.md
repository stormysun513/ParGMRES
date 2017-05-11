---
title: Project Report
permalink: report.html
---


Project Report
======================

## Title

ParGMRES: A parallel linear solver

## Team member

- Yu-Lun Tsai (yulunt)
- Chih-Wei Chang (cchang3)

## Summary

short (no more than a paragraph) project summary. If applicable, the summary 
should list your project deliverables (including what you plan to show at the 
parallelism competition) and what machines they ran on.

## Background

Describe the algorithm, application, or system you parallelized in computer 
science terms. (Recall our discussion from the last day of class.) Figure(s) 
would be really useful here.

## Approach

Tell us how your implementation works. Your description should be sufficiently 
detailed to provide the course staff a basic understanding of your approach. 
Again, it might be very useful to include a figure here illustrating components 
of the system and/or their mapping to parallel hardware.

## Result

We run our code on two sets of matrix: BCSPWR and Cage. The size of these matrix
varies cross a very large range, in which we can better test the scalability of
our algorithm. Also, the difference between the size scale can also show us the
effectiveness of parallization.

Matrix details for BCSPWR:
- BCSPWR01: 39-by-39, 85 non-zero entries.
- BCSPWR03: 118-by-118, 297 non-zero entries.
- BCSPWR06: 1454-by-1454, 3377 non-zero entries.

The results of BCSPWR is shown as following, the y-axis is linear-scale.
![CPU_CAGE](imgs/cpu_bcspwr.png) 

Matrix details for Cage:
- Cage4: 9-by-9, 49 non-zero entries.
- Cage5: 37-by-37, 233 non-zero entries.
- Cage7: 340-by-340, 3084 non-zero entries.
- Cage8: 1015-by-1015, 11,003 non-zero entries.
- Cage9: 3534-by-3534, 41,594 non-zero entries.
- Cage10: 11,397-by-11,397, 150,645 non-zero entries.

The results of Cage is shown as following, the y-axis is log-scale.
![CPU_CAGE](imgs/cpu_cage.png) 

## References

Please provide a list of references used in the project.
