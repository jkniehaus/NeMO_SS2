#!/usr/bin/env nextflow

process DownloadAndPrepareFastq {
    tag "$sample_id"

    input:
    path manifest_file        // .txt file with URLs
    val sample_id             // e.g. ACA

    output:
    path "fq", emit: fastq_dir
    path "seq_batch.tsv", emit: batch_file

    script:
    """
    mkdir fq tarfiles
    cp ${manifest_file} downloadlist.txt
    cp downloadlist.txt tarfiles/

    parent=\$PWD

    # Up to 5 download attempts
    for i in {1..5}; do
        if [[ \$(wc -l < downloadlist.txt) -gt 1 ]]; then
            cd tarfiles
            xargs -P 2 -n 1 curl -O < downloadlist.txt
            echo "Untarring files"
            ls *.tar | xargs -n1 tar -xvf
            find . -name '*.fastq.gz' -exec mv -t \$parent/fq {} +

            cd \$parent/fq

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
            sed 's/...........$//' < fqsUse.txt | sort | uniq | awk '{cell=\$1; print cell"R1.fastq.gz", cell"R2.fastq.gz", cell}' OFS="\\t" > seq_batch.tsv

            if [[ \$(wc -l < downloadlist.txt) == \$(wc -l < seq_batch.tsv) ]]; then
                echo "Downloads complete"
                break
            else
                echo "Download error: regenerating downloadlist.txt"
                cut -f1 seq_batch.tsv > temp.txt
                sed -i 's/.\$//' temp.txt
                grep -v -f temp.txt -i downloadlist.txt > newlist.txt
                mv newlist.txt downloadlist.txt
                cp downloadlist.txt tarfiles/
            fi
        fi
    done
    """
}
