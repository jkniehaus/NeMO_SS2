#!/usr/bin/env nextflow

process dlNemo {
    tag "$sample_id"
    cpus 2
    time '10d'
    queue 'general'

    // Limit to 5 concurrent downloads
    maxForks 5

    // Each sample always uses its own fixed output folder
    storeDir "downloads/${sample_id}"
    publishDir "downloads/${sample_id}", mode: 'copy'
    cache 'deep'

    input:
    tuple path(manifest_file), val(sample_id)

    output:
    tuple val(sample_id), path ("*.fastq.gz"), path ("fqs.txt"), emit: dl_out

    script:
    """
    set -euo pipefail
    cp ${manifest_file} downloadlist.txt
    if [ \$(wc -l < downloadlist.txt) -gt 1 ]; then
        xargs -P 2 -n 1 curl -C - -O < downloadlist.txt
        echo "Untarring files"
        for f in *.tar; do
            tar -xvf "\$f"
            rm -f "\$f"
        done

        echo "Collecting FASTQ files"
        find . -type f -name '*.fastq.gz' | while read -r fq; do
            base=\$(basename "\$fq")
            if [ "\$fq" != "./\$base" ]; then
                if [ -e "./\$base" ]; then
                    uuid=\$(uuidgen | cut -c1-8)
                    mv "\$fq" "./\${uuid}_\$base"
                else
                    mv "\$fq" .
                fi
            fi
        done

        find . -maxdepth 1 -name '*.fastq.gz' > fqs.txt
    fi
    """
}
