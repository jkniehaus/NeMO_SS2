#!/usr/bin/env nextflow

process dlNemo {
    scratch true
    tag "$sample_id"
    cpus 2
    time '10d'

    input:
    tuple path(manifest_file), val(sample_id)

    output:
    tuple val(sample_id), path ("*.fastq.gz"), path ("fqs.txt"), emit: dl_out
    script:
    """
    set -euo pipefail

    cp ${manifest_file} downloadlist.txt

    if [ \$(wc -l < downloadlist.txt) -gt 1 ]; then
        xargs -P 2 -n 1 curl -O < downloadlist.txt
        echo "Untarring files"
        ls *.tar | xargs -n1 tar -xvf
        find . -name '*.fastq.gz' -exec mv -t . {} +
        find . -iname "*fastq.gz" > fqs.txt
    fi
    """
}
