process mergePeReads {
    input:
        tuple val(sampleId), path(fastq1), path(fastq2)
    
    output:
        tuple val(sampleId), path("mergepe_${sampleId}.fq"), emit: mergePe
    
    script:
    """
    echo "Merging paired-end reads for sample ${sampleId}..."
    ${params.seqtk} mergepe "${fastq1}" "${fastq2}" > "mergepe_${sampleId}.fq"
    """
}

process runPangenie {
    input:
        tuple val(sampleId), 
        path(mergedFastq) 
    
    output:
        tuple val(sampleId), path("${sampleId}_genotyping.vcf"), emit: pangenie_out
    
    script:
    def pangenie_index_path = params["${params.ref}_pangenie_index_path"]
    """
    echo "Running PanGenie for sample ${sampleId}..."
    ${params.pangenie} \
        -f "${pangenie_index_path}" \
        -i "${mergedFastq}" \
        -t "${params.threads}" \
        -j "${params.threads}" \
        -o "${sampleId}"
    """
}

process Pangenie_postprocess{
    publishDir "${params.outdir}/03.SVs/${sampleId}", mode: 'link', overwrite: true
    input:
        tuple val(sampleId),path(pangenie_vcf)
    output:
        tuple val(sampleId),path("${sampleId}_${params.ref}_pangenie_PanVariants_SVs.vcf.gz"),path("${sampleId}_${params.ref}_pangenie_PanVariants_SVs.vcf.gz.tbi"), emit: pangenie_PanVariants_SVs
    script:
        def refpath = params["${params.ref}_refpath"]
        """
        mkdir -p \$PWD/tmp
        export TMPDIR=\$PWD/tmp
        # 压缩和索引输入VCF
        ${params.bgzip} -c ${pangenie_vcf} > tmp.vcf.gz
        ${params.tabix} -p vcf tmp.vcf.gz

        # 过滤和排序
        ${params.bcftools} view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 tmp.vcf.gz \
            | ${params.bcftools} view -i 'STRLEN(REF)>50 || MAX(STRLEN(ALT))>50' \
            | ${params.bcftools} view -e 'GT[*]="0/0"' \
            | ${params.bcftools} sort -Oz -o "tmp_gt50.vcf.gz"

        # 标准化VCF文件
        ${params.bcftools} norm -m -any "tmp_gt50.vcf.gz" \
            | ${params.bcftools} view -e 'GT[*]="0/0"' > "tmp_gt50_bialle.vcf"
        
        # 生成GRCh38头部
        echo '##fileformat=VCFv4.2' > new_header.txt 
        # 重构VCF头部并处理后续步骤
        awk '{print "##contig=<ID="\$1",length="\$2">"}' "${refpath}.fai" \
            | grep -v _decoy | grep -v _random | grep -v chrUn_ \
            | grep -v chrEBV | grep -v chrM | sort -t '#' -k3,3V >> new_header.txt
        ${params.bcftools} view -h tmp_gt50_bialle.vcf | grep -v "##contig" >> new_header.txt
        output_PREFIX=${sampleId}_${params.ref}_pangenie
        ${params.bcftools} reheader -h new_header.txt tmp_gt50_bialle.vcf > \${output_PREFIX}_rawSV.vcf
        ${params.bgzip} \${output_PREFIX}_rawSV.vcf
        ${params.tabix} -p vcf \${output_PREFIX}_rawSV.vcf.gz

        # Truvari处理流程
        ${params.truvari} collapse \
            -i "\${output_PREFIX}_rawSV.vcf.gz" \
            -o "\${output_PREFIX}_truvari.filter.vcf" \
            -c "\${output_PREFIX}_collapsed.vcf" \
            -r 500 -p 0.95 -P 0.95 -s 50 -S 100000

        ${params.bcftools} view -i 'STRLEN(REF)>50 || MAX(STRLEN(ALT))>50' "\${output_PREFIX}_truvari.filter.vcf" \
            | ${params.bcftools} view -i 'FORMAT/GQ[*] != "."' \
            | ${params.bcftools} sort -Oz -o "\${output_PREFIX}_truvari.filter_GQ.filter.vcf.gz"
        ${params.tabix} -p vcf "\${output_PREFIX}_truvari.filter_GQ.filter.vcf.gz"

        ${params.bcftools} view "\${output_PREFIX}_truvari.filter_GQ.filter.vcf.gz" | \
        awk -F'\t' -v threshold=10 '
          /^#/ { print; next }
          {
            d = length(\$4) - length(\$5)
            if (d < 0) d = -d
            if (d > threshold) print
          }' > "\${output_PREFIX}_truvari.filter_GQ.filter_length.filter.vcf"

        ${params.bgzip} "\${output_PREFIX}_truvari.filter_GQ.filter_length.filter.vcf"

        mv "\${output_PREFIX}_truvari.filter_GQ.filter_length.filter.vcf.gz" "\${output_PREFIX}_PanVariants_SVs.vcf.gz"
        ${params.tabix} -p vcf "\${output_PREFIX}_PanVariants_SVs.vcf.gz"
        """
}

process run_manta{
    publishDir "${params.outdir}/03.SVs/${sampleId}",
               mode: 'link', overwrite: true
    input:
        tuple val(sampleId),path(bam),path(bam_bai)
    output:
        tuple val(sampleId),path("${sampleId}_sv_calls.vcf.gz"),path("${sampleId}_sv_calls.vcf.gz.tbi"), emit: manta_sv_vcf
        val "run_manta_done", emit: run_manta_done
    script:
        def only_chromosome_fa = params["${params.ref}_only_chromosome_fa"]
        """
        mkdir -p manta_out
        ${params.configManta_py} \
            --bam ${bam} \
            --referenceFasta ${only_chromosome_fa} \
            --runDir manta_out
        ./manta_out/runWorkflow.py -m local -j ${task.cpus}
        mv manta_out/results/variants/diploidSV.vcf.gz ${sampleId}_sv_calls.vcf.gz
        mv manta_out/results/variants/diploidSV.vcf.gz.tbi ${sampleId}_sv_calls.vcf.gz.tbi
        """
}

process truvariBench {
    publishDir "${params.outdir}/031.SVs_eval",
               mode: 'link', overwrite: true

    input:
        tuple val(sampleId), path(SVs_query), path(SVs_query_tbi)
    
    output:
        tuple val(sampleId), path("${sampleId}_vs_${params.svs_basename}/*"), emit:SVs_eval
        val "truvariBench_done", emit: truvariBench_done
    script:
    def refpath = params["${params.ref}_refpath"]
    def svs_basevcf = params["${params.ref}_svs_basevcf"]
    def svs_basebed = params["${params.ref}_svs_basebed"]
    """
    echo "Running truvari benchmark for sample ${sampleId}..."
    mkdir -p \$PWD/tmp
    export TMPDIR=\$PWD/tmp
    ${params.truvari} bench \
        -f "${refpath}" \
        -b "${svs_basevcf}" \
        -c "${SVs_query}" \
        -o "${sampleId}_vs_${params.svs_basename}" \
        --includebed "${svs_basebed}" \
        -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 --passonly
    """
}
