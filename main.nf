#!/usr/bin/env nextflow

params.genome1 = '/proj/gs25/users/Jesse/references/mm39_exonSpecificSMARTseq'
params.genome2 = '/proj/gs25/users/Jesse/references/mm39_starsolo_ssv4_MOR1UTR'

include { DownloadAndPrepareFastq } from './modules/dlNemo.nf'
include { ssAlign } from './modules/ssAlign.nf'

workflow {
    Channel.fromPath("manifests/*.txt")
        .map { file -> tuple(file, file.baseName) }  // e.g. ACA.txt -> ("manifests/ACA.txt", "ACA")
        |>
        DownloadAndPrepareFastq
        |>
        ssAlign
}
