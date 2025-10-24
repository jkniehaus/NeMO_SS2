#!/usr/bin/env nextflow

params.genome1 = '/proj/gs25/users/Jesse/references/mm39_exonSpecificSMARTseq'
params.genome2 = '/proj/gs25/users/Jesse/references/mm39_starsolo_ssv4_MOR1UTR'

include { DownloadAndPrepareFastq } from './modules/dlNemo.nf'
include { ssAlign } from './modules/ssAlign.nf'
include { MapCells } from './modules/aligned2mmc.nf'
include { download_sra } from './modules/dlsra.nf'
include { trim_reads } from './modules/trimreads.nf'
include { alignQuant } from './modules/normAlign.nf'
include { altalignQuant } from './modules/altAlign.nf'

workflow {
    Channel.fromPath("manifests/*.txt")
        .map { file -> tuple(file, file.baseName) }  // e.g. ACA.txt -> ("manifests/ACA.txt", "ACA")
        |>
        DownloadAndPrepareFastq
    ssAlign(
        DownloadAndPrepareFastq.fastq_dir,
        DownloadAndPrepareFastq.batch_file,
        params.genome1,
        params.genome2,
        DownloadAndPrepareFastq.batch_file.map { it.simpleName }
    )
    MapCells(
        DownloadAndPrepareFastq.batch_file.map { it.simpleName }
    )
}
