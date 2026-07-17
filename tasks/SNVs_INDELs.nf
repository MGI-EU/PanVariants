// KMC process
process kmcProcess {
    publishDir "${params.outdir}/00.FQ/${sampleId}",
               mode: 'link', overwrite: true,
               saveAs: { fn -> fn.endsWith('.fq.paths') ? fn : null }

    input:
        tuple val(sampleId), path(fastq1), path(fastq2)

    output:
        path("${sampleId}.fq.paths"), emit: fq_paths
        tuple val(sampleId), path(fastq1), path(fastq2), path("${sampleId}.fq.kmc.kff"),   emit: kmc_out

    script:
    """
    cat > ${sampleId}.fq.paths <<EOF
    ${fastq1}
    ${fastq2}
    EOF

    ${params.kmc} -k29 -m32 -okff -t${params.threads} \
        @${sampleId}.fq.paths \
        ${sampleId}.fq.kmc \
        \$PWD
    """
}

// vg paths creation
// need large memory, cannot run in local mode
process vgPaths {
    publishDir "${params.outdir}/00.FQ/${sampleId}",
               mode: 'link', overwrite: true

    input:
    tuple val(sampleId), path(fastq1), path(fastq2)

    output:
    tuple val(sampleId), path("${sampleId}.${params.ref}.path_list.txt"), emit: path_list  // add val(sampleId)

    script:
    def gbz = params["${params.ref}_gbz"]
    """
    ${params.vg} paths -x ${gbz} -L -Q ${params.ref} \
        | grep -v _decoy | grep -v _random | grep -v chrUn_ | grep -v chrEBV | sort -t '#' -k3,3V \
        > ${sampleId}.${params.ref}.path_list.txt
    """
}
// grep -x 'GRCh38#0#chr1'
// delete chr1 if you want to run on all chromosomes
// keep chrM

// vg giraffe mapping
process vgGiraffe {
    input:
    tuple val(sampleId), path(fastq1), path(fastq2), path(kmc_out), path(path_list)

    output:
    tuple val(sampleId), path("${sampleId}.unsorted.bam"), emit: unsort_bam

    script:
    def gbz = params["${params.ref}_gbz"]
    def hapl = params["${params.ref}_hapl"]
    """
    ${params.vg} giraffe --progress \
    --read-group "ID:${params.id} LB:${params.lb} SM:${sampleId} PL:${params.pl} PU:${params.pu}" \
    --sample "${sampleId}" \
    -o BAM \
    --ref-paths "${path_list}" \
    -P -L 3000 \
    -f "${fastq1}" -f "${fastq2}" \
    -Z "${gbz}" \
    --kff-name "${kmc_out}" \
    --haplotype-name "${hapl}" \
    --index-basename "${sampleId}" \
    -t "${task.cpus}" > ${sampleId}.unsorted.bam
    """
}

// samtools sort
process samtoolsSort {
    publishDir "${params.outdir}/01.BAM/${sampleId}", mode:'link'
    input:
    tuple val(sampleId), path(unsort_bam)

    output:
    tuple val(sampleId), 
    path("${sampleId}.${params.ref}.PanVariants.sorted.bam"),
    path("${sampleId}.${params.ref}.PanVariants.sorted.bam.bai"),  emit: sorted_bam

    script:
    def sorted_BAM="${sampleId}.${params.ref}.PanVariants.sorted.bam"
    """
    echo "Sample prefix: ${sampleId}"
    echo "Sorted BAM file name: ${sorted_BAM}"
    ${params.samtools} view -h ${unsort_bam} \
      | sed -e 's/${params.ref}#0#//g' \
      | ${params.samtools} sort -@ ${params.threads} -m 1G -O BAM -o ${sorted_BAM}
    ${params.samtools} index ${sorted_BAM}
    """
}
// Extract paired reads where both ends are unmapped
process extract_unmapped_reads {
    publishDir "${params.outdir}/01.BAM/${sampleId}", mode:'link'

    input:
    tuple val(sampleId),path(sorted_bam),path(sorted_bam_bai)

    output:
    tuple val(sampleId),path("${sampleId}_both_unmapped.bam"), emit: unmapped_reads
    script:
    """
    ${params.samtools} view -@ ${params.threads} -f 0xC -hb ${sorted_bam} > ${sampleId}_both_unmapped.bam
    """

}
// Sort by read name and convert to FASTQ
process sort_and_convert_to_fastq{
    publishDir "${params.outdir}/01.BAM/${sampleId}", mode:'link'

    input:
    tuple val(sampleId),path(unmapped_bam)

    output:
    tuple val(sampleId),path("${sampleId}_R1.fastq"),path("${sampleId}_R2.fastq"),stdout, emit: unmapped_fq

    script:
    """
        BOTH_COUNT=\$(${params.samtools} view -c ${unmapped_bam})
        if [ \$BOTH_COUNT -gt 0 ]; then
            ${params.samtools} index ${unmapped_bam} &&  \
            ${params.samtools} sort -n -@ ${params.threads} -o ${sampleId}_unmapped_sorted.bam ${unmapped_bam} && \
            ${params.samtools} fastq -@ ${params.threads} \
            -1 ${sampleId}_R1.fastq \
            -2 ${sampleId}_R2.fastq \
            -0 /dev/null \
            -s /dev/null \
            -n ${sampleId}_unmapped_sorted.bam
        else
            touch ${sampleId}_R1.fastq ${sampleId}_R2.fastq
        fi
        echo \$BOTH_COUNT

    """
}
// Realign with bwa
process bwa_align{
    publishDir "${params.outdir}/01.BAM/${sampleId}", mode:'link'
    input:
    tuple val(sampleId),path(unmapped_fq1),path(unmapped_fq2),val(unmapped_reads_num)
    output:
    tuple val(sampleId),path("${sampleId}_bwa_realigned.bam"), emit: bwa_realigned_bam
    script:
    def refpath = params["${params.ref}_bwa_refpath"]
    """
    ${params.bwa_mem2} mem -t ${params.threads} \
    -Y \
    ${refpath} \
    ${unmapped_fq1} \
    ${unmapped_fq2} 2> bwa.log | \
    ${params.samtools} view -@ ${params.threads} -b -o ${sampleId}_bwa_realigned.bam

    """
}
// Extract originally mapped reads
process extract_mapped_reads {
    publishDir "${params.outdir}/01.BAM/${sampleId}", mode:'link'
    input:
    tuple val(sampleId),path(sorted_bam),path(sorted_bam_bai),val(unmapped_reads_num)

    output:
    tuple val(sampleId),path("${sampleId}_vg_mapped.bam"), emit: vg_mapped_bam

    script:
    """
    ${params.samtools} view -@ ${params.threads} -F 0xC -b ${sorted_bam} > ${sampleId}_vg_mapped.bam
    """
}
// Merge BAM files
process merge_bam {
    publishDir "${params.outdir}/01.BAM/${sampleId}", mode:'link'

    input:
    tuple val(sampleId),path(vg_mapped_bam),path(bwa_realigned_bam)

    output:
    tuple val(sampleId),path("${sampleId}_merged.bam"), emit: merged_bam

    script:
    """
    ${params.samtools} merge -@ ${params.threads} -f ${sampleId}_merged.bam \
    ${vg_mapped_bam} \
    ${bwa_realigned_bam}
    """
}
// Sort merged BAM by coordinate
process sort_merge_bam {
    publishDir "${params.outdir}/01.BAM/${sampleId}", mode:'link'

    input:
    tuple val(sampleId),path(merged_bam)

    output:
    tuple val(sampleId),path("${sampleId}_merged_sorted.bam"),path("${sampleId}_merged_sorted.bam.bai"), emit: merged_sorted_bam

    script:
    """
    ${params.samtools} sort -@ ${params.threads} -o ${sampleId}_merged_sorted.bam ${merged_bam}
    ${params.samtools} index -@ ${params.threads} ${sampleId}_merged_sorted.bam
    """
}
// Count alignment statistics
process bam_stat{
    publishDir "${params.outdir}/01.BAM/${sampleId}", mode:'link'

    input:
    tuple val(sampleId),path(sorted_bam),path(sorted_bai),path(merged_bam),path(merged_bam_bai)

    output:
    tuple val(sampleId),path("${sampleId}_vg_bam_stat.txt"),path("${sampleId}_merge_bam_stat.txt"), emit: bamstat
    
    script:
    """
    ${params.samtools} flagstat ${sorted_bam} | head -5 > ${sampleId}_vg_bam_stat.txt
    ${params.samtools} flagstat ${merged_bam} | head -5 > ${sampleId}_merge_bam_stat.txt
    """

}

process samtools_markdup_dv{
    input:
        tuple val(sample_name),path(sort_bam),path(sort_bam_bai)

    output:
        tuple val(sample_name),path("${sample_name}_markdup.bam"),path("${sample_name}_markdup.bam.bai"), emit: markdup_bam

    script:
        """
        ${params.samtools} sort -@ ${params.threads} -n ${sort_bam} | \
        ${params.samtools} fixmate -@ ${params.threads} -m - - | \
        ${params.samtools} sort -@ ${params.threads} - | \
        ${params.samtools} markdup -@ ${params.threads} - ${sample_name}_markdup.bam && \
        ${params.samtools} index ${sample_name}_markdup.bam
        # rm -f `realpath ${sort_bam}`
        # rm -f `realpath ${sort_bam_bai}`
        """
}

process split_bam{
    input:
        tuple val(sampleId),path(bam),path(bam_bai)
        each chr_num

    output:
        tuple val(sampleId),val(chr_num),path("${sampleId}_${chr_num}.bam"),path("${sampleId}_${chr_num}.bam.bai"), emit: chr_bam_out

    script:
    """
    ${params.samtools} view -b \
        -@ ${task.cpus} \
        ${bam} \
        -o ${sampleId}_${chr_num}.bam \
        ${chr_num} && \
    ${params.samtools} index ${sampleId}_${chr_num}.bam
    """
}


// DeepVariant
process deepVariant_chr{
    publishDir "${params.outdir}/02.SNVs_INDELs/${sampleId}",mode: 'link', overwrite: true

    input:
    tuple val(sampleId), val(chr_num), path(chr_bam), path(chr_bam_bai)

    output:
    tuple val(sampleId),
    path("${sampleId}.${params.ref}_${chr_num}.PanVariants.SNVs_INDELs.vcf.gz"),
    path("${sampleId}.${params.ref}_${chr_num}.PanVariants.SNVs_INDELs.vcf.gz.tbi"), emit: vcf_chr_query
    tuple val(sampleId),
    path("${sampleId}.${params.ref}_${chr_num}.PanVariants.SNVs_INDELs.gvcf.gz"),
    path("${sampleId}.${params.ref}_${chr_num}.PanVariants.SNVs_INDELs.gvcf.gz.tbi"), emit: gvcf_chr_query
    val "deepVariant_chr_done", emit: deepVariant_chr_done
    script:
    def gbz = params["${params.ref}_gbz"]
    def refpath = params["${params.ref}_refpath"]
    """
    make_examples_extra_args="min_mapping_quality=0,keep_legacy_allele_counter_behavior=true,normalize_reads=true"
    if [[ ${params.ref} == "CHM13" ]];then
        make_examples_extra_args="ref_name_pangenome=CHM13,\$make_examples_extra_args"
    fi
    mkdir -p \$PWD/tmp
    export TMPDIR=\$PWD/tmp
    /opt/deepvariant/bin/run_pangenome_aware_deepvariant \
      --model_type "${params.dvmodel}" \
      --ref "${refpath}" \
      --reads "${chr_bam}" \
      --regions "${chr_num}" \
      --customized_model "${params.customized_model}" \
      --pangenome "${gbz}" \
      --gbz_shared_memory_name "${sampleId}_share" \
      --output_vcf "\$PWD/${sampleId}.${params.ref}_${chr_num}_vg1.66.0_pandv1.9.0.vcf.gz" \
      --output_gvcf "\$PWD/${sampleId}.${params.ref}_${chr_num}_vg1.66.0_pandv1.9.0.gvcf.gz" \
      --runtime_report True \
      --intermediate_results_dir \$PWD/tmp \
      --vcf_stats_report True \
      --num_shards ${task.cpus} \
      --postprocess_variants_extra_args="only_keep_pass=true" \
      --make_examples_extra_args=\"\$make_examples_extra_args\"
    
    gunzip -c "\$PWD/${sampleId}.${params.ref}_${chr_num}_vg1.66.0_pandv1.9.0.vcf.gz" \
    | awk '/^##DeepVariant_version=1.9.0\$/ && !seen {print "##PanVariants_version=1.0.0"; seen=1} 1' \
    | bgzip -c > "\$PWD/${sampleId}.${params.ref}_${chr_num}.PanVariants.SNVs_INDELs.vcf.gz" && \
    tabix -p vcf "\$PWD/${sampleId}.${params.ref}_${chr_num}.PanVariants.SNVs_INDELs.vcf.gz"

    gunzip -c "\$PWD/${sampleId}.${params.ref}_${chr_num}_vg1.66.0_pandv1.9.0.gvcf.gz" \
    | awk '/^##DeepVariant_version=1.9.0\$/ && !seen {print "##PanVariants_version=1.0.0"; seen=1} 1' \
    | bgzip -c > "\$PWD/${sampleId}.${params.ref}_${chr_num}.PanVariants.SNVs_INDELs.gvcf.gz" && \
    tabix -p vcf "\$PWD/${sampleId}.${params.ref}_${chr_num}.PanVariants.SNVs_INDELs.gvcf.gz"
    """
}

process collect_vcf{
    executor 'local'
    input:
        tuple val(sampleId),path(vcfs),path(vcf_tbis)
    
    output:
        tuple val(sampleId),val(merge_vcf_input), emit: merge_vcf_input_out
    
    exec:
        merge_vcf_input = ""
        for (vcf in vcfs){
            vcf_input = " -V ${vcf}"
            merge_vcf_input = merge_vcf_input + vcf_input
        }

}

process merge_vcf{
    publishDir "${params.outdir}/02.SNVs_INDELs/${sampleId}",mode: 'link', overwrite: true
    input:
        tuple val(sampleId),val(merge_vcf_input),path(vcfs),path(vcf_tbis)
    
    output:
        tuple val(sampleId),path("${sampleId}.${params.ref}.PanVariants.SNVs_INDELs.vcf.gz"),path("${sampleId}.${params.ref}.PanVariants.SNVs_INDELs.vcf.gz.tbi"), emit: merge_vcf
        tuple val(sampleId),val("${sampleId}_merge_vcf_done"), emit: merge_vcf_done

    script:
        def refpath = params["${params.ref}_refpath"]
        """
        ${params.java} -Xmx${task.memory.giga}g \
         -Djava.io.tmpdir=./tmp \
         -cp ${params.GATK} \
         org.broadinstitute.gatk.tools.CatVariants \
         -R ${refpath} \
         ${merge_vcf_input} \
         -out ${sampleId}.${params.ref}.PanVariants.SNVs_INDELs.vcf.gz  && \
        ${params.tabix} -p vcf ${sampleId}.${params.ref}.PanVariants.SNVs_INDELs.vcf.gz --force
        """
}

process deepVariant {
    publishDir "${params.outdir}/02.SNVs_INDELs/${sampleId}",mode: 'link', overwrite: true

    input:
    tuple val(sampleId), path(sorted_bam), path(sorted_bai)

    output:
    tuple val(sampleId),
    path("${sampleId}.${params.ref}.PanVariants.SNVs_INDELs.vcf.gz"),
    path("${sampleId}.${params.ref}.PanVariants.SNVs_INDELs.vcf.gz.tbi"), emit: vcf_query
    tuple val(sampleId),
    path("${sampleId}.${params.ref}.PanVariants.SNVs_INDELs.gvcf.gz"),
    path("${sampleId}.${params.ref}.PanVariants.SNVs_INDELs.gvcf.gz.tbi"), emit: gvcf_query
    val "deepVariant_done", emit: deepVariant_done
    script:
    def gbz = params["${params.ref}_gbz"]
    def refpath = params["${params.ref}_refpath"]
    """
    make_examples_extra_args="min_mapping_quality=0,keep_legacy_allele_counter_behavior=true,normalize_reads=true"
    if [[ ${params.ref} == "CHM13" ]];then
        make_examples_extra_args="ref_name_pangenome=CHM13,\$make_examples_extra_args"
    fi
    mkdir -p \$PWD/tmp
    export TMPDIR=\$PWD/tmp
    /opt/deepvariant/bin/run_pangenome_aware_deepvariant \
      --model_type "${params.dvmodel}" \
      --ref "${refpath}" \
      --reads "${sorted_bam}" \
      --customized_model "${params.customized_model}" \
      --pangenome "${gbz}" \
      --gbz_shared_memory_name "${sampleId}_share" \
      --output_vcf "\$PWD/${sampleId}.${params.ref}_vg1.66.0_pandv1.9.0.vcf.gz" \
      --output_gvcf "\$PWD/${sampleId}.${params.ref}_vg1.66.0_pandv1.9.0.gvcf.gz" \
      --runtime_report True \
      --intermediate_results_dir \$PWD/tmp \
      --vcf_stats_report True \
      --num_shards ${task.cpus} \
      --postprocess_variants_extra_args="only_keep_pass=true" \
      --make_examples_extra_args=\"\$make_examples_extra_args\"
    
    gunzip -c "\$PWD/${sampleId}.${params.ref}_vg1.66.0_pandv1.9.0.vcf.gz" \
    | awk '/^##DeepVariant_version=1.9.0\$/ && !seen {print "##PanVariants_version=1.0.0"; seen=1} 1' \
    | bgzip -c > "\$PWD/${sampleId}.${params.ref}.PanVariants.SNVs_INDELs.vcf.gz" && \
    tabix -p vcf "\$PWD/${sampleId}.${params.ref}.PanVariants.SNVs_INDELs.vcf.gz"

    gunzip -c "\$PWD/${sampleId}.${params.ref}_vg1.66.0_pandv1.9.0.gvcf.gz" \
    | awk '/^##DeepVariant_version=1.9.0\$/ && !seen {print "##PanVariants_version=1.0.0"; seen=1} 1' \
    | bgzip -c > "\$PWD/${sampleId}.${params.ref}.PanVariants.SNVs_INDELs.gvcf.gz" && \
    tabix -p vcf "\$PWD/${sampleId}.${params.ref}.PanVariants.SNVs_INDELs.gvcf.gz"
    """
}
// --regions chr1
// delete chr1 if you want to run on all chromosomes
// --gbz_shared_memory_size_gb 16 when using hprcv2.0

// vcfeval
process vcfEval {
    publishDir "${params.outdir}/021.SNVs_INDELs_eval/${sampleId}",mode: 'link', overwrite: true

    input:
    tuple val(sampleId), path(vcf_query), path(vcf_query_index)

    output:
    tuple val(sampleId), path("${sampleId}.${params.ref}.PanVariants.SNVs_INDELs.summary.csv"), emit: vcf_eval

    script:
    def refpath = params["${params.ref}_refpath"]
    def snvs_indels_basevcf = params["${params.ref}_snvs_indels_basevcf"]
    def snvs_indels_basebed = params["${params.ref}_snvs_indels_basebed"]
    def snvs_indels_basetsv = params["${params.ref}_snvs_indels_basetsv"]
    def snvs_indels_basesdf = params["${params.ref}_snvs_indels_basesdf"]
    """
    ${params.python2} ${params.happy} \
        --threads "$params.threads" \
        "${snvs_indels_basevcf}" "${vcf_query}" \
        -f "${snvs_indels_basebed}" \
        --stratification ${snvs_indels_basetsv} \
        --engine=vcfeval \
        --engine-vcfeval-template "${snvs_indels_basesdf}" \
        -r "${refpath}" \
        -o "${sampleId}.${params.ref}.PanVariants.SNVs_INDELs"
    """
}

process run_rtgtools{
    publishDir "${params.outdir}/021.SNVs_INDELs_eval/${sampleId}",mode: 'link', overwrite: true
    input:
        tuple val(sampleId), path(vcf_query), path(vcf_query_index)
    output:
        tuple val(sampleId),path("*"), emit: rtg_out
        val "run_rtgtools_done", emit: run_rtgtools_done
    script:
        def snvs_indels_basevcf = params["${params.ref}_snvs_indels_basevcf"]
        def snvs_indels_basebed = params["${params.ref}_snvs_indels_basebed"]
        def snvs_indels_basetsv = params["${params.ref}_snvs_indels_basetsv"]
        def snvs_indels_basesdf = params["${params.ref}_snvs_indels_basesdf"]
        """
        ${params.rtg} vcfeval \
          -b ${snvs_indels_basevcf} \
          -c ${vcf_query} \
          -e ${snvs_indels_basebed} \
          -o ${sampleId}.PanVariants.SNVs_INDELs \
          -t ${snvs_indels_basesdf}
        """
}
