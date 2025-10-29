#!/usr/bin/env nextflow

process mapCells {
    scratch true
    tag "$batch_prefix"
    cpus 16
    memory '32 GB'
    time '12h'

    publishDir 'output', mode: 'copy'

    input:
    tuple val(batch_prefix), path(batch_prefix_path), path(batch_prefix2), path(ssRaw2h5ad)

    output:
    tuple val(batch_prefix), path("rds/${batch_prefix}.rds")

    script:
    """
    set -euo pipefail

    echo "Running R script to generate Seurat and h5ad files..."
    Rscript ${ssRaw2h5ad} -n 16 -m 32 -o ${batch_prefix_path} -d ${batch_prefix2}

    echo "Running mapmycells annotation..."
    python -m cell_type_mapper.cli.from_specified_markers \\
        --query_path h5ads/${batch_prefix}.h5ad \\
        --extended_result_path rds/${batch_prefix}.json \\
        --csv_result_path rds/${batch_prefix}.csv \\
        --drop_level CCN20230722_SUPT \\
        --cloud_safe False \\
        --query_markers.serialized_lookup /proj/gs25/users/Jesse/references/mapmycells/mouse_markers_230821.json \\
        --precomputed_stats.path /proj/gs25/users/Jesse/references/mapmycells/precomputed_stats_ABC_revision_230821.h5 \\
        --type_assignment.normalization raw \\
        --type_assignment.n_processors 16
    """
}
