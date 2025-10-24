#!/usr/bin/env nextflow

process MapCells {
    tag "$batch_prefix"
    cpus 16
    memory '32 GB'
    time '12h'

    input:
    val batch_prefix

    output:
    path "rds/${batch_prefix}.rds"
    path "rds/${batch_prefix}.json"
    path "rds/${batch_prefix}.csv"
    path "h5ads/${batch_prefix}.h5ad"

    script:
    """
    set -euo pipefail

    echo "Running R script to generate Seurat and h5ad files..."
    Rscript ssRaw2h5ad.R -n 16 -m 32 -p ${batch_prefix}

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
    rm -r tarfiles_${batch_prefix}
    rm -r fq_${batch_prefix}
    """
}
