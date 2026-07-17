process get_UsageReport{
    shell '/bin/bash', '-u'

    input:
        val usage_input

    output:
        path "usage_report.done", emit: usage_report_done

    script:
    """
    ${params.python3} ${params.generateUsageReport_py} ${params.outdir}/../
    cp ${params.outdir}/../runit_report.html ${params.outdir}/upload
    touch usage_report.done
    """
}

process cleanup_intermediate{
    shell '/bin/bash', '-ue'

    input:
        path cleanup_input

    script:
    """
    sample_ids=\$(awk 'BEGIN{FS="\\t"} /^[[:space:]]*(#|\$)/ {next} {print \$1}' ${params.sampleList} | sort -u)
    work_root=\$(realpath ${params.outdir}/../work 2>/dev/null || true)

    while IFS= read -r sample_id; do
        [ -z "\$sample_id" ] && continue

        rm -f \
            "${params.outdir}/00.FQ/\${sample_id}/\${sample_id}.merge.fq.1.gz" \
            "${params.outdir}/00.FQ/\${sample_id}/\${sample_id}.merge.fq.2.gz" \
            "${params.outdir}/01.BAM/\${sample_id}/\${sample_id}.${params.ref}.PanVariants.sorted.bam" \
            "${params.outdir}/01.BAM/\${sample_id}/\${sample_id}.${params.ref}.PanVariants.sorted.bam.bai"

        if [ -n "\$work_root" ] && [ -d "\$work_root" ]; then
            find "\$work_root" -type f '(' \
                -name "\${sample_id}.merge.fq.1.gz" -o \
                -name "\${sample_id}.merge.fq.2.gz" -o \
                -name "\${sample_id}.${params.ref}.PanVariants.sorted.bam" -o \
                -name "\${sample_id}.${params.ref}.PanVariants.sorted.bam.bai" \
            ')' -delete
        fi
    done <<< "\$sample_ids"
    """
}
