process STR_filter{
    //STR filter: PASS, non-0/0, |len change| >= 50 bp
    publishDir "${params.outdir}/05.STRs/${sample_id}", mode:'link'
    input:
        tuple val(sample_id),path(str_vcf)
    output:
        //tuple val(sample_id),path("str.adjformat.${sample_id}.vcf.gz"),path("str.adjformat.${sample_id}.vcf.gz.tbi"), emit: str_adjformat_vcf
        tuple val(sample_id),path("str.filtered.${sample_id}.vcf"), emit: str_filtered_vcf
    script:
        """
        ${params.bcftools} view -f PASS -i 'GT!="0/0"' ${str_vcf} \
        | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {
        refcnt = ""
        if (match(\$8, /REF=[0-9]+/)) {
            tmp = substr(\$8, RSTART, RLENGTH)
            gsub(/REF=/, "", tmp)
            refcnt = tmp
        }

        ru = ""
        if (match(\$8, /RU=[^;]+/)) {
            tmp = substr(\$8, RSTART, RLENGTH)
            gsub(/RU=/, "", tmp)
            ru = tmp
        }

        if (refcnt == "" || ru == "") next

        rul = length(ru)
        maxabs = 0
        nalt = split(\$5, alt, ",")

        for (i = 1; i <= nalt; i++) {
            if (alt[i] ~ /^<STR[0-9]+>\$/) {
                tmp_str = alt[i]
                gsub(/[<>STR]/, "", tmp_str)
                d = (tmp_str - refcnt) * rul
                if (d < 0) d = -d
                if (d > maxabs) maxabs = d
            }
        }

        if (maxabs >= 50) print
        }' > str.filtered.${sample_id}.vcf
        """
}

process adjformat{
    //Converts and normalizes an ExpansionHunter VCF into 
    // an adjusted bgzipped VCF with index.
    publishDir "${params.outdir}/05.STRs/${sample_id}", mode:'link'
    input:
        tuple val(sample_id),path(str_filtered_vcf)
    output:
        tuple val(sample_id),path("str.adjformat.${sample_id}.vcf.gz"),path("str.adjformat.${sample_id}.vcf.gz.tbi"), emit: str_adjformat_vcf
        val "adjformat_done", emit: adjformat_done
    script:
        def refpath = params["${params.ref}_refpath"]
        """
        mkdir -p \$PWD/tmp
        export TMPDIR=\$PWD/tmp
        # Step 1: Convert VCF
        # ===========================
        ${params.python3} "${params.converter_script}" -i "${str_filtered_vcf}" -o "${sample_id}_convert.vcf"

        # ===========================
        # Step 2: Normalize and filter
        # ===========================
        ${params.bgzip} -f "${sample_id}_convert.vcf"
        ${params.tabix} -f -p vcf "${sample_id}_convert.vcf.gz"

        ${params.bcftools} filter -i 'FILTER="PASS"' "${sample_id}_convert.vcf.gz" \
        | ${params.bcftools} norm -N -d any \
        | ${params.bcftools} norm -N -m-any -o "${sample_id}_split.vcf"

        # ===========================
        # Step 3: Truvarizer annotate
        # ===========================
        ${params.python3} "${params.truvarizer_script}" -i "${sample_id}_split.vcf" -o "${sample_id}_fillin.vcf"

        ${params.bcftools} sort "${sample_id}_fillin.vcf" -o "${sample_id}_sorted.vcf"
        ${params.bgzip} -f "${sample_id}_sorted.vcf"
        ${params.tabix} -f -p vcf "${sample_id}_sorted.vcf.gz"

        # ===========================
        # Step 4: Fill REF from FASTA
        # ===========================
        ${params.bcftools} +fill-from-fasta "${sample_id}_sorted.vcf.gz" -- -c REF -f "${refpath}" > "${sample_id}_fix.vcf"

        ${params.bgzip} -f "${sample_id}_fix.vcf"
        ${params.tabix} -f -p vcf "${sample_id}_fix.vcf.gz"

        # ===========================
        # Step 5: Move final outputs, cleanup handled by trap
        # ===========================
        mv -f "${sample_id}_fix.vcf.gz"     "str.adjformat.${sample_id}.vcf.gz"
        mv -f "${sample_id}_fix.vcf.gz.tbi" "str.adjformat.${sample_id}.vcf.gz.tbi"
        """
}

process SV_add_type{
    // Ensure SVTYPE in panvariants
    publishDir "${params.outdir}/03.SVs/${sample_id}", mode:'link'
    input:
        tuple val(sample_id),path(SV_vcf),path(SV_vcf_tbi)
    output:
        tuple val(sample_id),path("${sample_id}_pv.addSVTYPE.vcf.gz"),path("${sample_id}_pv.addSVTYPE.vcf.gz.tbi"), emit: pv_with_svtype_vcf
        val "SV_add_type_done", emit: SV_add_type_done
    script:
        """
        echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' > svtype.hdr
        ${params.bcftools} annotate -h svtype.hdr -O z -o ${sample_id}_pv.addSVTYPE.vcf.gz ${SV_vcf}
        ${params.tabix} -p vcf ${sample_id}_pv.addSVTYPE.vcf.gz
        """
}

process sample_concat{
    // Unify sample names and Concat manta + cnvnator + STR
    input:
        tuple val(sample_id),path(manta_vcf),path(manta_vcf_tbi),path(cnvnator_vcf),path(cnvnator_vcf_tbi),path(str_adjformat_vcf),path(str_adjformat_vcf_tbi),path(pv_with_svtype),path(pv_with_svtype_tbi)
    output:
        tuple val(sample_id),path("${sample_id}_manta_CNV_STR.concat.vcf.gz"),path("${sample_id}_manta_CNV_STR.concat.vcf.gz.tbi"), emit: concat_vcf
        tuple val(sample_id),path("${sample_id}_pv.s.vcf.gz"),path("${sample_id}_pv.s.vcf.gz.tbi"), emit: pv_s_vcf
    script:
        """
        echo "${sample_id}" > "sample.name"
        reheader_if_needed () {
          in="\$1"; out="\$2"
          ${params.bcftools} reheader -s "sample.name" "\$in" | ${params.bcftools} view -Oz -o "\$out"
          ${params.tabix} -f -p vcf "\$out"
        }
        reheader_if_needed ${manta_vcf} ${sample_id}_manta.s.vcf.gz
        reheader_if_needed ${cnvnator_vcf} ${sample_id}_cnv.s.vcf.gz
        reheader_if_needed ${str_adjformat_vcf} ${sample_id}_str.s.vcf.gz
        reheader_if_needed ${pv_with_svtype} ${sample_id}_pv.s.vcf.gz

        ${params.bcftools} concat -a \
            -O z \
            -o ${sample_id}_manta_CNV_STR.concat.vcf.gz \
            ${sample_id}_manta.s.vcf.gz \
            ${sample_id}_cnv.s.vcf.gz \
            ${sample_id}_str.s.vcf.gz
        ${params.tabix} -p vcf ${sample_id}_manta_CNV_STR.concat.vcf.gz
        """
}

process vcf_normalize_collapse{
    // Normalize + collapse + sort
    input:
        tuple val(sample_id),path(concat_vcf),path(concat_vcf_tbi)
    output:
        tuple val(sample_id),path("${sample_id}_manta_CNV_STR.sorted.vcf.gz"),path("${sample_id}_manta_CNV_STR.sorted.vcf.gz.tbi"), emit: manta_CNV_STR_sorted_vcf
    script:
        def refpath = params["${params.ref}_refpath"]
        """
        ${params.bcftools} norm \
            --check-ref s \
            --fasta-ref ${refpath} \
            -N -m-any \
            -O z \
            -o ${sample_id}_manta_CNV_STR.norm.vcf.gz ${concat_vcf}
        ${params.tabix} -p vcf ${sample_id}_manta_CNV_STR.norm.vcf.gz
        ${params.truvari} collapse \
            -i ${sample_id}_manta_CNV_STR.norm.vcf.gz \
            -o ${sample_id}_manta_CNV_STR.collapse.vcf \
            -c ${sample_id}_manta_CNV_STR.collapse.log
        ${params.bcftools} sort \
            -T ./ \
            -m 4G -Oz \
            -o ${sample_id}_manta_CNV_STR.sorted.vcf.gz \
           ${sample_id}_manta_CNV_STR.collapse.vcf
        ${params.tabix} -p vcf ${sample_id}_manta_CNV_STR.sorted.vcf.gz
        """
}

process exclude_overlap{
    //Exclude pv TPs overlapping manta/cnv/str TPs, then concat pv + (mcs)
    input:
        tuple val(sample_id),path(pv_s_vcf),path(pv_s_vcf_tbi),path(manta_CNV_STR_sorted_vcf),path(manta_CNV_STR_sorted_vcf_tbi)
    output:
        tuple val(sample_id),path("${sample_id}_pv_plus_mantaCNVSTR.concat.vcf.gz"),path("${sample_id}_pv_plus_mantaCNVSTR.concat.vcf.gz.tbi"), emit: pv_mantaCNVSTR_vcf
    script:
        """
        mkdir -p \$PWD/tmp
        export TMPDIR=\$PWD/tmp
        ${params.truvari} bench \
            -b ${pv_s_vcf} \
            -c ${manta_CNV_STR_sorted_vcf} \
            -o ${sample_id}_bench --passonly
        ${params.bcftools} isec -C -p ${sample_id}_isec ${pv_s_vcf} ${sample_id}_bench/tp-base.vcf.gz
        ${params.bgzip} -c ${sample_id}_isec/0000.vcf > ${sample_id}_pv.filtered.vcf.gz
        ${params.tabix} -p vcf ${sample_id}_pv.filtered.vcf.gz

        echo "${sample_id}" > "sample.name"
        reheader_if_needed () {
          in="\$1"; out="\$2"
          ${params.bcftools} reheader -s "sample.name" "\$in" | ${params.bcftools} view -Oz -o "\$out"
          ${params.tabix} -f -p vcf "\$out"
        }
        reheader_if_needed ${sample_id}_pv.filtered.vcf.gz ${sample_id}_pv.filtered.s.vcf.gz
        ${params.bcftools} concat -a \
            -O z \
            -o ${sample_id}_pv_plus_mantaCNVSTR.concat.vcf.gz \
            ${sample_id}_pv.filtered.s.vcf.gz \
            ${sample_id}_manta_CNV_STR.sorted.vcf.gz
        ${params.tabix} -p vcf ${sample_id}_pv_plus_mantaCNVSTR.concat.vcf.gz
        """
}

process reheader_merge_vcf{
    publishDir "${params.outdir}/03.SVs/${sample_id}", mode:'link'
    input:
        tuple val(sample_id),path(pv_mantaCNVSTR_vcf),path(pv_mantaCNVSTR_vcf_tbi)
    output:
        tuple val(sample_id),path("${sample_id}.merged.vcf.gz"),path("${sample_id}.merged.vcf.gz.tbi"), emit: merged_vcf
        val "reheader_merge_vcf_done", emit: reheader_merge_vcf_done
    script:
        def refpath = params["${params.ref}_refpath"]
        """
        ${params.bcftools} norm \
            --check-ref s \
            --fasta-ref ${refpath} \
            -N -m-any \
            -O z \
            -o ${sample_id}_pv_plus_mantaCNVSTR.norm.vcf.gz ${pv_mantaCNVSTR_vcf}
        ${params.tabix} -p vcf ${sample_id}_pv_plus_mantaCNVSTR.norm.vcf.gz
        
        ${params.truvari} collapse \
            -i ${sample_id}_pv_plus_mantaCNVSTR.norm.vcf.gz \
            -o ${sample_id}_pmcs.collapse.vcf \
            -c ${sample_id}_pmcs.collapse.log
        ${params.bcftools} sort \
            -T ./ \
            -m 4G \
            -Oz \
            -o ${sample_id}_pv_plus_mantaCNVSTR.norm.sorted.vcf.gz \
            ${sample_id}_pmcs.collapse.vcf
        gzip -d ${sample_id}_pv_plus_mantaCNVSTR.norm.sorted.vcf.gz

        ${params.bcftools} view ${sample_id}_pv_plus_mantaCNVSTR.norm.sorted.vcf \
        | awk 'BEGIN{OFS="\t"}
          /^##/ {print; next}
          /^#CHROM/{
            print "##INFO=<ID=ShortSV,Number=1,Type=String,Description=\\"TRUE if |REF-ALT| length <=50bp, otherwise FALSE\\">";
            print;
            next
          }
          {
            ref=\$4; alt=\$5;
            d = length(ref) - length(alt);
            if (d < 0) d = -d;
            tag = (d < 50 ? "TRUE" : "FALSE");

            if (\$8 == ".")
              \$8 = "ShortSV=" tag;
            else
              \$8 = \$8 ";ShortSV=" tag;

            print
          }' | ${params.bgzip} -c > ${sample_id}_pv_plus_mantaCNVSTR.norm.sorted.ShortSV.vcf.gz
        ${params.tabix} -p vcf ${sample_id}_pv_plus_mantaCNVSTR.norm.sorted.ShortSV.vcf.gz

        ${params.bcftools} view -h ${sample_id}_pv_plus_mantaCNVSTR.norm.sorted.ShortSV.vcf.gz \
        | grep -v -E '^##bcftools|^##cmdline=|^##reference=|^##source=|^##fileDate=' \
        | awk '
          /^##fileformat/ {f[1]=f[1]\$0"\\n"; next}
          /^##contig/     {f[2]=f[2]\$0"\\n"; next}
          /^##FILTER/     {f[3]=f[3]\$0"\\n"; next}
          /^##ALT/        {f[4]=f[4]\$0"\\n"; next}
          /^##INFO/       {f[5]=f[5]\$0"\\n"; next}
          /^##FORMAT/     {f[6]=f[6]\$0"\\n"; next}
          /^##/           {f[7]=f[7]\$0"\\n"; next}
          /^#CHROM/       {c=\$0}
          END{printf "%s%s%s%s%s%s%s%s", f[1],f[2],f[3],f[4],f[5],f[6],f[7],c"\\n"}
        ' > ${sample_id}_clean_header.txt
        
        ${params.bcftools} reheader \
            -h ${sample_id}_clean_header.txt \
            ${sample_id}_pv_plus_mantaCNVSTR.norm.sorted.vcf \
            | ${params.bgzip} -c > ${sample_id}.merged.vcf.gz
        ${params.tabix} -p vcf ${sample_id}.merged.vcf.gz
        """
}
