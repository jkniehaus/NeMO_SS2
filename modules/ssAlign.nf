#!/usr/bin/env nextflow

process ssAlign {
    tag "$batch_prefix"
    cpus 16
    memory '64 GB'
    time '12h'

    input:
    path fastq_dir
    path batch_file
    val genome_dir1
    val genome_dir2
    val batch_prefix

    output:
    path "${batch_prefix}*", optional: true

    script:
    """
    set -euo pipefail
    echo "Aligning with genome 1..."
    STAR --genomeDir ${genome_dir1} \
        --runThreadN 16 \
        --runDirPerm All_RWX \
        --readFilesCommand zcat \
        --outSAMtype None \
        --soloType SmartSeq \
        --readFilesManifest ${batch_file}
        --soloUMIdedup Exact \
        --soloStrand Unstranded \
        --soloFeatures GeneFull_ExonOverIntron \
        --outFileNamePrefix ./${batch_prefix}_
        --outTmpDir tmp_${batch_prefix} \
        --soloMultiMappers Unique \
        --soloOutFileNames ${batch_prefix} features.tsv barcodes.tsv matrix.mtx
    if [ -n "${genome_dir2}" ]; then
      STAR --genomeDir ${genome_dir2} \
          --runThreadN 16 \
          --runDirPerm All_RWX \
          --readFilesCommand zcat \
          --outSAMtype None \
          --outFileNamePrefix ./${batch_prefix}2
          --soloType SmartSeq \
          --readFilesManifest ${batch_file}
          --soloUMIdedup Exact \
          --soloStrand Unstranded \
          --soloFeatures GeneFull_ExonOverIntron \
          --outTmpDir tmp2_${batch_prefix} \
          --soloMultiMappers Unique \
          --soloOutFileNames ${batch_prefix}2 features.tsv barcodes.tsv matrix.mtx
    fi
    """
}
