#!/usr/bin/env bash
# alignAll.sh
outDir='results/quant/'

# Iterate over all files that match the pattern Aip*.R1.fastq in the specified directory
for file in /work/courses/BINF6309/AiptasiaMiSeq/fastq/Aip*.R1.fastq; do
    # Extract the sample name from the file name
    sample=$(basename "${file}" .R1.fastq)

    # Run Salmon for each sample
    function align {
        salmon quant -l IU \
            -1 /work/courses/BINF6309/AiptasiaMiSeq/fastq/${sample}.R1.fastq \
            -2 /work/courses/BINF6309/AiptasiaMiSeq/fastq/${sample}.R2.fastq \
            -i AipIndex \
            --validateMappings \
            -o ${outDir}${sample}
    }

    align
done


