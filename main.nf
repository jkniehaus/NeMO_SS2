#!/usr/bin/env nextflow

params.genome1 = '/proj/gs25/users/Jesse/references/mm39_exonSpecificSMARTseq'
params.genome2 = '/proj/gs25/users/Jesse/references/mm39_starsolo_ssv4_MOR1UTR'

include { dlNemo } from './modules/dlNemo.nf'
include { ssAlign } from './modules/ssAlign.nf'
include { mapCells } from './modules/aligned2mmc.nf'
workflow {
    manifest_ch = Channel.fromPath("manifests/*.txt")
                         .map { file -> tuple(file, file.baseName) }

    dl_out_ch = dlNemo(manifest_ch)

    (align_out1_ch, align_out2_ch) = ssAlign(dl_out_ch, params.genome1, params.genome2)

    map_in_ch = align_out1_ch.join(align_out2_ch, by: 0)
                             .map { sample_id, dir1, dir2 ->
                                 tuple(sample_id, dir1, dir2, file("${projectDir}/resources/ssRaw2h5ad.R"), file("${projectDir}/resources/mmc2Seurat.R"))
                             }

    mapCells(map_in_ch)
}
