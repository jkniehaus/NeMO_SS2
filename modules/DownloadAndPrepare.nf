#!/usr/bin/env nextflow

process DownloadAndPrepareFastq {
    tag "$sample_id"
    cpus 2
    time '12h'

    input:
    path manifest_file        // .txt file with URLs
    val sample_id             // e.g. ACA

    output:
    path "fq_${sample_id}", emit: fastq_dir
    path "seq_batch_${sample_id}.tsv", emit: batch_file

    script:
    """
    set -euo pipefail
    mkdir fq_${sample_id} tarfiles_${sample_id}
    cp ${manifest_file} tarfiles_${sample_id}/downloadlist.txt
    cp ${manifest_file} fq_${sample_id}/downloadlist.txt
    parent=\$PWD
    cd tarfiles_${sample_id}
    # Up to 5 download attempts
    for i in {1..5}; do
        if [[ \$(wc -l < downloadlist.txt) -gt 1 ]]; then
            xargs -P 2 -n 1 curl -O < downloadlist.txt
            echo "Untarring files"
            ls *.tar | xargs -n1 tar -xvf
            find . -name '*.fastq.gz' -exec mv -t \$parent/fq_${sample_id} {} +

            cd \$parent/fq_${sample_id}

            find * -iname "*fastq.gz" > fqs.txt

            # Keep only R1/R2 pairs
            input_file="fqs.txt"
            output_file="fqsUse.txt"
            rm -f "\${output_file}"
            while IFS= read -r line; do
                if [[ \${line} =~ R1\\.fastq\\.gz ]]; then
                    prefix=\$(echo "\${line}" | sed 's/R1\\.fastq\\.gz//')
                    r1_file="\${prefix}R1.fastq.gz"
                    r2_file="\${prefix}R2.fastq.gz"
                    if grep -q -F "\${r1_file}" "\${input_file}" && grep -q -F "\${r2_file}" "\${input_file}"; then
                        echo "\${r1_file}" >> "\${output_file}"
                        echo "\${r2_file}" >> "\${output_file}"
                    fi
                fi
            done < "\${input_file}"

            # Construct seq_batch.tsv
            sed 's/...........$//' < fqsUse.txt | sort | uniq | awk '{cell=\$1; print cell"R1.fastq.gz", cell"R2.fastq.gz", cell}' OFS="\\t" > seq_batch_${sample_id}.tsv

            if [[ \$(wc -l < downloadlist.txt) == \$(wc -l < seq_batch_${sample_id}.tsv) ]]; then
                echo "Downloads complete"
                break
            else
                echo "Download error: regenerating downloadlist.txt"
                cut -f1 seq_batch_${sample_id}.tsv > temp.txt
                sed -i 's/.\$//' temp.txt
                grep -v -f temp.txt -i downloadlist.txt > newlist.txt
                mv newlist.txt downloadlist.txt
                cp downloadlist.txt ${parent}/tarfiles_${sample_id}/downloadlist.txt
            fi
        fi
    done
    """
}
