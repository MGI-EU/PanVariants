/*
 * CNV calling workflow using CNVnator
 */


process runCNVnator{
    publishDir "${params.outdir}/04.CNVs/${sample_id}", mode: 'link'
    input:
        tuple val(sample_id),path(BAM),path(BAI)
    output:
        tuple val(sample_id),path("${sample_id}.cnv.txt"), emit: cnv_txt
    script:
        def refpath = params["${params.ref}_refpath"]
        """
        source ${params.thisroot_sh}
        # 步骤1: 从BAM构建读分布树
        echo "[1/5] Building read distribution tree..."
        ${params.cnvnator} -root ${sample_id}.root -tree ${BAM} -chrom \$(seq -f 'chr%g' 1 22) chrX chrY

         步骤2: 生成读深度直方图
        echo "[2/5] Generating read depth histogram (bin_size=${params.cnv_bin_size})..."
        ${params.cnvnator} -root ${sample_id}.root -his ${params.cnv_bin_size} -fasta ${refpath}

         步骤3: 计算统计量
        echo "[3/5] Calculating statistics..."
        ${params.cnvnator} -root ${sample_id}.root -stat ${params.cnv_bin_size}

         步骤4: 分割染色体区域
        echo "[4/5] Partitioning chromosomes..."
        ${params.cnvnator} -root ${sample_id}.root -partition ${params.cnv_bin_size}

         步骤5: 调用CNV并转换为VCF
        echo "[5/5] Calling CNVs and generating VCF..."
        ${params.cnvnator} -root ${sample_id}.root -call ${params.cnv_bin_size} > ${sample_id}.cnv.txt
        """
}

process cnv_filter{
    publishDir "${params.outdir}/04.CNVs/${sample_id}", mode: 'link'
    input:
        tuple val(sample_id),path(cnv_txt)
    output:
        tuple val(sample_id),path("${sample_id}.PanVariants_CNVs.vcf.gz"),path("${sample_id}.PanVariants_CNVs.vcf.gz.tbi"), emit: cnv_filtered_vcf
        val "cnv_filter", emit: cnv_filter_done
    script:
        """
        # 过滤CNV结果
        ${params.perl} ${params.cnvnator_filter_pl} ${cnv_txt} > ${sample_id}.PanVariants_CNVs.txt

        # 转换为VCF格式
        ${params.perl} ${params.cnvnator2VCF_pl} ${sample_id}.PanVariants_CNVs.txt -prefix ${sample_id} > ${sample_id}.PanVariants_CNVs.vcf

        # 压缩和索引VCF
        echo "Compressing and indexing VCF..."
        ${params.bgzip} -f ${sample_id}.PanVariants_CNVs.vcf
        ${params.tabix} -p vcf ${sample_id}.PanVariants_CNVs.vcf.gz

        # 过滤gap区域的CNV
        if [ -f "${sample_id}_gap_region.bed" ]; then
            echo "Filtering CNVs in gap regions..."
            # 提取CNV位置信息
            ${params.bcftools} query -f '%CHROM\t%POS\t%END\n' ${sample_id}.PanVariants_CNVs.vcf.gz > ${sample_id}_cnvs.bed

            # 找出与gap区域重叠的CNV
            if command -v ${params.bedtools} &> /dev/null; then
                ${params.bedtools} intersect -a ${sample_id}_cnvs.bed -b ${sample_id}_gap_region.bed -wa > ${sample_id}_overlapping_cnvs.bed

                # 过滤掉重叠的CNV
                if command -v ${params.vcftools} &> /dev/null; then
                    ${params.vcftools} --gzvcf ${sample_id}.PanVariants_CNVs.vcf.gz --exclude-bed ${sample_id}_overlapping_cnvs.bed --recode --recode-INFO-all --out ${sample_id}_filtered_cnvs

                    # 替换原始VCF文件
                    mv ${sample_id}_filtered_cnvs.recode.vcf ${sample_id}.PanVariants_CNVs.vcf
                    ${params.bgzip} -f ${sample_id}.PanVariants_CNVs.vcf
                    ${params.tabix} -p vcf "${sample_id}.PanVariants_CNVs.vcf.gz"
                fi
            fi
        fi
        """
}
process cnvBench {
    tag "CNV Benchmark - ${sample_id}"
    publishDir "${params.outdir}/04.CNVs/${sample_id}", mode: 'link'
    
    input:
    tuple val(sample_id), path(cnv_query_vcf), path(cnv_query_vcf_tbi)
    
    output:
    tuple val(sample_id), path("${sample_id}_cnv_benchmark.txt"), emit: cnv_benchmark
    
    script: 
    """
    # 使用truvari进行CNV benchmark分析
    ${params.wittyer} \
    -i ${cnv_query_vcf} \
    -t ${params.cnvs_basevcf} \
    --configFile ${params.truvari_config} \
    --includeBed ${params.cnvs_basebed} \
    -em CrossTypeAndSimpleCounting \
    -o "\$PWD"
     
    """
}
