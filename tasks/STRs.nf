process runExpansionHunter{
    publishDir "${params.outdir}/05.STRs/${sample_id}", mode: 'link'
    input:
        tuple val(sample_id),path(bam),path(bam_bai)
    output:
        tuple val(sample_id),path("${sample_id}.json"), emit: str_json
        tuple val(sample_id),path("${sample_id}.vcf"), emit: str_vcf
    script:
        def refpath = params["${params.ref}_refpath"]
        def variant_catalog = params["${params.ref}_variant_catalog"]
        """
        ${params.ExpansionHunter} --reads ${bam} \
                --reference ${refpath} \
                --variant-catalog ${variant_catalog} \
                --output-prefix ${sample_id} \
                --sex female \
                --threads 16 \
                --analysis-mode seeking
        """
}