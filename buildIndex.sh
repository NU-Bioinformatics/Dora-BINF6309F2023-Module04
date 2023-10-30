#!/usr/bin/env bash
# buildIndex.sh
salmon index --transcripts /work/courses/BINF6309/data_BINF6309/Module4/Trinity.fasta \
             --index AipIndex \
             --kmerLen 25
