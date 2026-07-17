process STR_filter{
    // STR filter: PASS and non-0/0, matching merge_260306.sh
    publishDir "${params.outdir}/05.STRs/${sample_id}", mode:'link'
    input:
        tuple val(sample_id),path(str_vcf)
    output:
        tuple val(sample_id),path("str.filtered.${sample_id}.vcf"), emit: str_filtered_vcf
    script:
        """
        ${params.bcftools} view -f PASS -i 'GT!="0/0"' ${str_vcf} > str.filtered.${sample_id}.vcf
        """
}

process adjformat{
    // Converts and normalizes an ExpansionHunter VCF into an adjusted bgzipped VCF with index.
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

        ${params.python3} "${params.converter_script}" -i "${str_filtered_vcf}" -o "${sample_id}_convert.vcf"

        ${params.bgzip} -f "${sample_id}_convert.vcf"
        ${params.tabix} -f -p vcf "${sample_id}_convert.vcf.gz"

        ${params.bcftools} filter -i 'FILTER="PASS"' "${sample_id}_convert.vcf.gz" \
        | ${params.bcftools} norm -N -d any \
        | ${params.bcftools} norm -N -m-any -o "${sample_id}_split.vcf"

        ${params.python3} "${params.truvarizer_script}" -i "${sample_id}_split.vcf" -o "${sample_id}_fillin.vcf"

        ${params.bcftools} sort "${sample_id}_fillin.vcf" -o "${sample_id}_sorted.vcf"
        ${params.bgzip} -f "${sample_id}_sorted.vcf"
        ${params.tabix} -f -p vcf "${sample_id}_sorted.vcf.gz"

        ${params.bcftools} +fill-from-fasta "${sample_id}_sorted.vcf.gz" -- -c REF -f "${refpath}" > "${sample_id}_fix.vcf"

        ${params.bgzip} -f "${sample_id}_fix.vcf"
        ${params.tabix} -f -p vcf "${sample_id}_fix.vcf.gz"

        mv -f "${sample_id}_fix.vcf.gz"     "str.adjformat.${sample_id}.vcf.gz"
        mv -f "${sample_id}_fix.vcf.gz.tbi" "str.adjformat.${sample_id}.vcf.gz.tbi"
        """
}

process prepare_merge_vcfs{
    // Build one compatible header, unify sample names, and annotate the source caller.
    input:
        tuple val(sample_id),
            path(manta_vcf),path(manta_vcf_tbi),
            path(cnvnator_vcf),path(cnvnator_vcf_tbi),
            path(str_adjformat_vcf),path(str_adjformat_vcf_tbi),
            path(pangenie_vcf),path(pangenie_vcf_tbi)
    output:
        tuple val(sample_id),path("${sample_id}.manta.annotate.vcf.gz"),path("${sample_id}.manta.annotate.vcf.gz.tbi"), emit: manta_annotated_vcf
        tuple val(sample_id),path("${sample_id}.cnvnator.annotate.vcf.gz"),path("${sample_id}.cnvnator.annotate.vcf.gz.tbi"), emit: cnvnator_annotated_vcf
        tuple val(sample_id),path("${sample_id}.str.annotate.vcf.gz"),path("${sample_id}.str.annotate.vcf.gz.tbi"), emit: str_annotated_vcf
        tuple val(sample_id),path("${sample_id}.pangenie.annotate.vcf.gz"),path("${sample_id}.pangenie.annotate.vcf.gz.tbi"), emit: pangenie_annotated_vcf
        tuple val(sample_id),path("${sample_id}.merged.header.txt"), emit: merged_header
    script:
        """
        { ${params.bcftools} view -h "${manta_vcf}"; \
          ${params.bcftools} view -h "${cnvnator_vcf}"; \
          ${params.bcftools} view -h "${str_adjformat_vcf}"; \
          ${params.bcftools} view -h "${pangenie_vcf}"; } \
        | grep -v -E '^##(bcftools|cmdline|source|reference|fileDate)' \
        | awk '
        function ID(p,s){s=\$0; sub("^"p"ID=","",s); sub(",.*","",s); sub(">.*","",s); return s}
        function KEEP(A,k){ if(!(k in A) || length(\$0)>length(A[k])) A[k]=\$0 }
        BEGIN{ haveC=haveS=0 }
        \$0~/^##INFO=<ID=SVLEN,/ {sub(/Number=[^,]*/,"Number=1")}
        \$0~/^##FORMAT=<ID=CN,/  {sub(/Number=[^,]*,Type=Integer/,"Number=1,Type=Float")}
        \$0~/^##contig=</ {KEEP(C,ID("##contig=<")); next}
        \$0~/^##FILTER=</ {KEEP(F,ID("##FILTER=<")); next}
        \$0~/^##ALT=</    {KEEP(A,ID("##ALT=<")); next}
        \$0~/^##INFO=</   {k=ID("##INFO=<"); KEEP(I,k); if(k=="CALLER")haveC=1; if(k=="ShortSV")haveS=1; next}
        \$0~/^##FORMAT=</ {KEEP(M,ID("##FORMAT=<")); next}
        \$0~/^##/ && \$0!~/^##fileformat=/ {O[\$0]=1; next}
        END{
          if(!haveC) I["CALLER"]="##INFO=<ID=CALLER,Number=1,Type=String,Description=\\"Source caller\\">";
          if(!haveS) I["ShortSV"]="##INFO=<ID=ShortSV,Number=1,Type=String,Description=\\"TRUE if |REF-ALT| length <=50bp, otherwise FALSE\\">";
          for(k in C) print C[k] > "H1"; for(k in F) print F[k] > "H2"; for(k in A) print A[k] > "H3";
          for(k in I) print I[k] > "H4"; for(k in M) print M[k] > "H5"; for(k in O) print k > "H6";
        }'

        touch H1 H2 H3 H4 H5 H6
        { echo "##fileformat=VCFv4.2"; \
          sort -u -V H1; sort -u H2; sort -u -V H3; sort -u H4; sort -u H5; sort -u H6; \
          printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t%s\\n' "${sample_id}"; } \
          > "${sample_id}.merged.header.txt"

        echo "${sample_id}" > sample.name

        reformat () {
          in=\"\$1\"; out=\"\$2\"; tag=\"\$3\"; caller=\"\$4\"
          rehdr=\"\${tag}.rehdr.vcf.gz\"

          ${params.bcftools} reheader -s sample.name -h "${sample_id}.merged.header.txt" "\$in" > "\$rehdr"
          ${params.tabix} -f -p vcf "\$rehdr"

          ${params.bcftools} view "\$rehdr" \
          | awk -v caller="\$caller" 'BEGIN{OFS="\\t"} /^#/ {print; next}
          {
            if(\$8=="."||\$8=="") \$8="CALLER="caller; else \$8=\$8";CALLER="caller;
            print
          }' | ${params.bgzip} -c > "\$out"
          ${params.tabix} -f -p vcf "\$out"
        }

        reformat "${manta_vcf}"         "${sample_id}.manta.annotate.vcf.gz"     manta manta
        reformat "${cnvnator_vcf}"      "${sample_id}.cnvnator.annotate.vcf.gz" cnv cnvnator
        reformat "${str_adjformat_vcf}" "${sample_id}.str.annotate.vcf.gz"       str expansionhunter
        reformat "${pangenie_vcf}"      "${sample_id}.pangenie.annotate.vcf.gz"  pg pangenie
        """
}

process filter_cnv_str_by_pangenie{
    // Keep only CNV/STR comparison calls matching the PanGenie baseline.
    input:
        tuple val(sample_id),
            path(cnv_annotated_vcf),path(cnv_annotated_vcf_tbi),
            path(str_annotated_vcf),path(str_annotated_vcf_tbi),
            path(pangenie_annotated_vcf),path(pangenie_annotated_vcf_tbi)
    output:
        tuple val(sample_id),path("${sample_id}.cnv_str.supported.vcf.gz"),path("${sample_id}.cnv_str.supported.vcf.gz.tbi"), emit: cnv_str_supported_vcf
        tuple val(sample_id),path("${sample_id}.pangenie.after_cnv_str.vcf.gz"),path("${sample_id}.pangenie.after_cnv_str.vcf.gz.tbi"), emit: pangenie_after_cnv_str_vcf
    script:
        """
        mkdir -p \$PWD/tmp
        export TMPDIR=\$PWD/tmp

        ${params.bcftools} concat -a -Ou "${cnv_annotated_vcf}" "${str_annotated_vcf}" \
        | ${params.bcftools} sort -m 4G -Oz -o "${sample_id}.cnv_str.sorted.vcf.gz"
        ${params.tabix} -f -p vcf "${sample_id}.cnv_str.sorted.vcf.gz"

        ${params.truvari} bench \
          -b "${pangenie_annotated_vcf}" \
          -c "${sample_id}.cnv_str.sorted.vcf.gz" \
          -o bench_pg_vs_cnv_str \
          -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 \
          --sizemax 100000 --passonly --pick multi

        cp -f bench_pg_vs_cnv_str/tp-comp.vcf.gz     "${sample_id}.cnv_str.supported.vcf.gz"
        cp -f bench_pg_vs_cnv_str/tp-comp.vcf.gz.tbi "${sample_id}.cnv_str.supported.vcf.gz.tbi"
        cp -f bench_pg_vs_cnv_str/fn.vcf.gz          "${sample_id}.pangenie.after_cnv_str.vcf.gz"
        cp -f bench_pg_vs_cnv_str/fn.vcf.gz.tbi      "${sample_id}.pangenie.after_cnv_str.vcf.gz.tbi"
        """
}

process filter_pangenie_by_new_variants{
    // Keep all Manta calls, add supported CNV/STR calls, and remove matching PanGenie calls.
    input:
        tuple val(sample_id),
            path(manta_annotated_vcf),path(manta_annotated_vcf_tbi),
            path(cnv_str_supported_vcf),path(cnv_str_supported_vcf_tbi),
            path(pangenie_after_cnv_str_vcf),path(pangenie_after_cnv_str_vcf_tbi)
    output:
        tuple val(sample_id),path("${sample_id}.new_variants.vcf.gz"),path("${sample_id}.new_variants.vcf.gz.tbi"), emit: new_variants_vcf
        tuple val(sample_id),path("${sample_id}.pangenie.filtered.vcf.gz"),path("${sample_id}.pangenie.filtered.vcf.gz.tbi"), emit: pangenie_filtered_vcf
    script:
        """
        mkdir -p \$PWD/tmp
        export TMPDIR=\$PWD/tmp

        ${params.bcftools} concat -a -Ou "${manta_annotated_vcf}" "${cnv_str_supported_vcf}" \
        | ${params.bcftools} sort -m 4G -Oz -o "${sample_id}.new_variants.vcf.gz"
        ${params.tabix} -f -p vcf "${sample_id}.new_variants.vcf.gz"

        ${params.truvari} bench \
          -b "${pangenie_after_cnv_str_vcf}" \
          -c "${sample_id}.new_variants.vcf.gz" \
          -o bench_pg_vs_new_variants \
          -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 \
          --sizemax 100000 --passonly --pick multi

        cp -f bench_pg_vs_new_variants/fn.vcf.gz     "${sample_id}.pangenie.filtered.vcf.gz"
        cp -f bench_pg_vs_new_variants/fn.vcf.gz.tbi "${sample_id}.pangenie.filtered.vcf.gz.tbi"
        """
}

process finalize_merged_sv{
    publishDir "${params.outdir}/03.SVs/${sample_id}", mode:'link', overwrite:true
    input:
        tuple val(sample_id),
            path(new_variants_vcf),path(new_variants_vcf_tbi),
            path(pangenie_filtered_vcf),path(pangenie_filtered_vcf_tbi),
            path(merged_header)
    output:
        tuple val(sample_id),path("${sample_id}.merged.vcf.gz"),path("${sample_id}.merged.vcf.gz.tbi"), emit: merged_vcf
        val "merge_SV_vcf_done", emit: merge_SV_vcf_done
    script:
        """
        mkdir -p \$PWD/tmp
        export TMPDIR=\$PWD/tmp

        ${params.bcftools} concat -a -Ou "${new_variants_vcf}" "${pangenie_filtered_vcf}" \
        | ${params.bcftools} sort -m 4G -Oz -o "${sample_id}.all.merged.sorted.vcf.gz"
        ${params.tabix} -f -p vcf "${sample_id}.all.merged.sorted.vcf.gz"

        ${params.bcftools} view "${sample_id}.all.merged.sorted.vcf.gz" \
        | awk 'BEGIN{OFS="\\t"} /^#/ {print; next}
        {
          n=split(\$8, info, ";");
          newinfo="";
          for(i=1; i<=n; i++){
            if(info[i] ~ /^(PctSeqSimilarity|PctSizeSimilarity|PctRecOverlap|SizeDiff|StartDistance|EndDistance|GTMatch|TruScore|MatchId)=/) continue;
            if(info[i] == "" || info[i] == ".") continue;
            newinfo = (newinfo=="" ? info[i] : newinfo ";" info[i]);
          }
          \$8 = (newinfo=="" ? "." : newinfo);

          ref=\$4; alt=\$5; d=length(ref)-length(alt); if(d<0)d=-d;
          tag=(d<=50 ? "TRUE" : "FALSE");
          if(\$8=="."||\$8=="") \$8="ShortSV="tag; else \$8=\$8";ShortSV="tag;

          print
        }' | ${params.bgzip} -c > "${sample_id}.all.merged.sorted.ShortSV.vcf.gz"
        ${params.tabix} -f -p vcf "${sample_id}.all.merged.sorted.ShortSV.vcf.gz"

        ${params.bcftools} reheader \
          -h "${merged_header}" \
          "${sample_id}.all.merged.sorted.ShortSV.vcf.gz" \
          > "${sample_id}.merged.vcf.gz"
        ${params.tabix} -f -p vcf "${sample_id}.merged.vcf.gz"
        """
}
