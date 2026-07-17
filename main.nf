include {timeStamp;readInput;merge_fq_PE;fqfilter} from './tasks/prepare.nf'
include {kmcProcess;vgPaths;vgGiraffe;samtoolsSort;split_bam;deepVariant_chr;
collect_vcf;merge_vcf;deepVariant;run_rtgtools;
extract_unmapped_reads;sort_and_convert_to_fastq;bwa_align;extract_mapped_reads
merge_bam;sort_merge_bam;bam_stat;samtools_markdup_dv} from './tasks/SNVs_INDELs.nf'
include {mergePeReads;runPangenie;Pangenie_postprocess;run_manta;truvariBench
} from './tasks/SVs.nf'
include {runCNVnator;cnvBench;cnv_filter} from './tasks/CNVs.nf'
include {runExpansionHunter} from './tasks/STRs.nf'
include {STR_filter;adjformat;prepare_merge_vcfs;filter_cnv_str_by_pangenie;
    filter_pangenie_by_new_variants;finalize_merged_sv
} from './tasks/merge_SV_vcf.nf'
include {get_UsageReport;cleanup_intermediate} from './tasks/generateUsageReport.nf'

input_list = readInput(params.sampleList)

// Main workflow
workflow {
    input_ctuple = Channel.fromList(input_list)
    process_done =  Channel.from("process done!")
    merge_fq_PE(input_ctuple)
    if ("${params.fq_filter}" == "yes"){
        fqfilter(merge_fq_PE.out.merge_fq)
        fq_input = fqfilter.out.clean_fq
    }else{
        fq_input = merge_fq_PE.out.merge_fq
    }
    
    kmcProcess(fq_input)
    vgPaths(fq_input)
    vgGiraffe(kmcProcess.out.kmc_out.join(vgPaths.out.path_list, by: 0))
    samtoolsSort(vgGiraffe.out.unsort_bam)
    if ("${params.mark_dup}" == "yes"){
        samtools_markdup_dv(samtoolsSort.out.sorted_bam)
        bam_input = samtools_markdup_dv.out.markdup_bam
    }else{
        bam_input = samtoolsSort.out.sorted_bam
    }
    if ("${params.re_alignment}" == "yes"){
        extract_unmapped_reads(bam_input)
        sort_and_convert_to_fastq(extract_unmapped_reads.out.unmapped_reads)
        sort_and_convert_to_fastq.out.unmapped_fq
            .filter(it -> it[3] != "0")
            .set{bwa_align_input}
        bwa_align(bwa_align_input)
        bwa_align_input.map{v -> tuple(v[0], v[3])}.set{unmapped_reads_num_out}
        bam_input
            .combine(unmapped_reads_num_out,by:0)
            .set{extract_mapped_reads_input}
        extract_mapped_reads(extract_mapped_reads_input)
        extract_mapped_reads.out.vg_mapped_bam
            .combine(bwa_align.out.bwa_realigned_bam,by:0)
            .set{merge_bam_input}
        merge_bam(merge_bam_input)
        sort_merge_bam(merge_bam.out.merged_bam)
        bam_input
            .combine(sort_merge_bam.out.merged_sorted_bam,by:0)
            .set{bam_stat_input}
        bam_stat(bam_stat_input)
        sort_merge_bam.out.merged_sorted_bam
            .ifEmpty(bam_input)
            .set{final_bam}
    }else{
        final_bam = bam_input
    }
    if (params.split_bam == "yes"){
        chr_list = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
        "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
        "chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"]
        chr_num_out = Channel.fromList(chr_list)
        split_bam(final_bam,chr_num_out)
        deepVariant_chr(split_bam.out.chr_bam_out)
        deepVariant_chr.out.vcf_chr_query
            .groupTuple(by:0, size:25)
            .set{collect_vcf_input}
        collect_vcf(collect_vcf_input)
        collect_vcf.out.merge_vcf_input_out
            .combine(collect_vcf_input,by:0)
            .set{merge_vcf_input}
        merge_vcf(merge_vcf_input)
        deepVariant_done = merge_vcf.out.merge_vcf_done
        deepVariant_vcf = merge_vcf.out.merge_vcf
    }else{
        deepVariant(final_bam)
        deepVariant_done = deepVariant.out.deepVariant_done
        deepVariant_vcf = deepVariant.out.vcf_query
    }
    

    // SNVs_INDELs benchmark - only executed when standard_sample = 1
    if (params.standard_sample == 1) {
        //vcfEval(deepVariant.out.vcf_query)
        run_rtgtools(deepVariant_vcf)
        run_rtgtools_done = run_rtgtools.out.run_rtgtools_done
    }else{
        run_rtgtools_done = process_done
    }

    // SV analysis workflow - only executed when run_sv_analysis is true
    if (params.run_sv == "yes") {
        mergePeReads(fq_input)
        runPangenie(mergePeReads.out.mergePe)
        Pangenie_postprocess(runPangenie.out.pangenie_out)
        Pangenie_postprocess.out.pangenie_PanVariants_SVs
            .map { "Pangenie_postprocess_done" }
            .set { Pangenie_postprocess_done }
        run_manta(bam_input)
        run_manta_done = run_manta.out.run_manta_done
        // SV benchmark - only executed when standard_sample = 1
        // if (params.standard_sample == 1) {
        //     truvariBench(runPangenie.out.pangenie_out)
        // }
    }else{
        Pangenie_postprocess_done = process_done
        run_manta_done = process_done
    }
        
    // CNV analysis workflow - only executed when cnv_calling is true
    if (params.run_cnv == "yes") {
        runCNVnator(final_bam)
        cnv_filter(runCNVnator.out.cnv_txt)
        cnv_filter_done = cnv_filter.out.cnv_filter_done
        // CNV benchmark - only executed when standard_sample = 1
        // CNV benchmark is currently disabled.
        // if (params.standard_sample == 1) {
        //     cnvBench(cnv_filter.out.cnv_filtered_vcf)
        // }
    }else{
        cnv_filter_done = process_done
    }

    // STR analysis workflow - only executed when run_str_analysis is true
    if (params.run_str == "yes"){
        runExpansionHunter(final_bam)
        STR_filter(runExpansionHunter.out.str_vcf)
        adjformat(STR_filter.out.str_filtered_vcf)
        adjformat_done = adjformat.out.adjformat_done
    }else{
        adjformat_done = process_done
    }

    // Merge SV, STR, CNV
    if (params.run_sv == "yes" && params.run_str == "yes" && params.run_cnv == "yes"){
        run_manta.out.manta_sv_vcf
            .combine(cnv_filter.out.cnv_filtered_vcf, by:0)
            .combine(adjformat.out.str_adjformat_vcf, by:0)
            .combine(Pangenie_postprocess.out.pangenie_PanVariants_SVs, by:0)
            .set{prepare_merge_vcfs_input}
        prepare_merge_vcfs(prepare_merge_vcfs_input)

        prepare_merge_vcfs.out.cnvnator_annotated_vcf
            .combine(prepare_merge_vcfs.out.str_annotated_vcf, by:0)
            .combine(prepare_merge_vcfs.out.pangenie_annotated_vcf, by:0)
            .set{filter_cnv_str_by_pangenie_input}
        filter_cnv_str_by_pangenie(filter_cnv_str_by_pangenie_input)

        prepare_merge_vcfs.out.manta_annotated_vcf
            .combine(filter_cnv_str_by_pangenie.out.cnv_str_supported_vcf, by:0)
            .combine(filter_cnv_str_by_pangenie.out.pangenie_after_cnv_str_vcf, by:0)
            .set{filter_pangenie_by_new_variants_input}
        filter_pangenie_by_new_variants(filter_pangenie_by_new_variants_input)

        filter_pangenie_by_new_variants.out.new_variants_vcf
            .combine(filter_pangenie_by_new_variants.out.pangenie_filtered_vcf, by:0)
            .combine(prepare_merge_vcfs.out.merged_header, by:0)
            .set{finalize_merged_sv_input}
        finalize_merged_sv(finalize_merged_sv_input)
        merge_SV_vcf_done = finalize_merged_sv.out.merge_SV_vcf_done
        if (params.standard_sample == 1) {
            truvariBench(finalize_merged_sv.out.merged_vcf)
            truvariBench_done = truvariBench.out.truvariBench_done
        }else{
            truvariBench_done = process_done
        }
    }else{
        merge_SV_vcf_done = process_done
        truvariBench_done = process_done
    }
    deepVariant_done
        .combine(run_rtgtools_done)
        .combine(Pangenie_postprocess_done)
        .combine(run_manta_done)
        .combine(cnv_filter_done)
        .combine(adjformat_done)
        .combine(merge_SV_vcf_done)
        .combine(truvariBench_done)
        .collect()
        .set{get_UsageReport_input}
    get_UsageReport(get_UsageReport_input)
    cleanup_intermediate(get_UsageReport.out.usage_report_done)
}
