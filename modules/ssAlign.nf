#!/usr/bin/env nextflow

process ssAlign {
    tag "$sample_id"
    cpus 16
    memory '64 GB'
    time '12h'

    input:
    tuple val(sample_id), path(fastqs), path(fq_list)
    val genome_dir1
    val genome_dir2

    output:
    tuple val(sample_id), path("${sample_id}/"), emit: outdir
    tuple val(sample_id), path("${sample_id}_g2/"), emit: outdir2
    path "seq_batch_${sample_id}.tsv"

    script:
    """
    set -euo pipefail
    input_file=$fq_list
    output_file="fqsUse.txt"
    while IFS= read -r line; do
        if [[ "\$line" =~ R1\\.fastq\\.gz ]]; then
            prefix=\$(echo "\$line" | sed 's/R1\\.fastq\\.gz//')
            r1_file="\${prefix}R1.fastq.gz"
            r2_file="\${prefix}R2.fastq.gz"
            if grep -q -F "\$r1_file" "\$input_file" && grep -q -F "\$r2_file" "\$input_file"; then
                echo "\$r1_file" >> "\$output_file"
                echo "\$r2_file" >> "\$output_file"
            fi
        fi
    done < "\$input_file"
    sed 's/...........\$//' < fqsUse.txt | sort | uniq | awk '{cell=\$1; print cell"R1.fastq.gz", cell"R2.fastq.gz", cell}' OFS="\\t" > seq_batch_${sample_id}.tsv
    echo "Aligning with genome 1..."
    STAR --genomeDir ${genome_dir1} \
        --runThreadN 16 \
        --runDirPerm All_RWX \
        --readFilesCommand zcat \
        --outSAMtype None \
        --soloType SmartSeq \
        --readFilesManifest seq_batch_${sample_id}.tsv \
        --soloUMIdedup Exact \
        --soloStrand Unstranded \
        --soloFeatures GeneFull_ExonOverIntron \
        --outFileNamePrefix ./${sample_id}/ \
        --outTmpDir tmp_${sample_id} \
        --soloMultiMappers Unique \
        --soloOutFileNames ${sample_id} features.tsv barcodes.tsv matrix.mtx
    if [ -n "${genome_dir2}" ]; then
        STAR --genomeDir ${genome_dir2} \
          --runThreadN 16 \
          --runDirPerm All_RWX \
          --readFilesCommand zcat \
          --outSAMtype None \
          --outFileNamePrefix ./${sample_id}_g2/ \
          --soloType SmartSeq \
          --readFilesManifest seq_batch_${sample_id}.tsv \
          --soloUMIdedup Exact \
          --soloStrand Unstranded \
          --soloFeatures Gene \
          --outTmpDir tmp2_${sample_id} \
          --soloMultiMappers Unique \
          --soloOutFileNames ${sample_id} features.tsv barcodes.tsv matrix.mtx
    fi
    """
}
