# Create a database folder within the workflow folder
DB=`pwd`/database
mkdir $DB

# ===================================================================
# Pangenome refererence genome: hprc-v1.1-mc-grch38
# ===================================================================
mkdir -p $DB/pangenome/grch38 && cd $DB/pangenome/grch38 
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.hapl

# ===================================================================
# linear reference genome: hg38/GRCh38
# ===================================================================
mkdir -p $DB/fa_path && cd $DB/fa_path
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz && \
gzip -d GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz &&
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
mkdir -p bwa_index && cd bwa_index && ln -s ../GCA_000001405.15_GRCh38_no_alt_analysis_set.fa ./ && \
bwa-mem2 index GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
cd $DB/fa_path
CHROMS=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrM" "chrX" "chrY")
REF_FILE="GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
for chrom in "${CHROMS[@]}"; do
    samtools faidx "$REF_FILE" "$chrom" >> GCA_000001405.15_GRCh38_only_chromosome.fa
    echo "$chrom"
done
sleep 10
samtools faidx GCA_000001405.15_GRCh38_only_chromosome.fa

# ===================================================================
# PanVariants DNBSEQ model: T1+
# ===================================================================
mkdir -p $DB/deepvariant_model && cd $DB/deepvariant_model
wget https://storage.googleapis.com/deepvariant/complete-case-study-testdata/complete-t1plus/2026/checkpoint-141056-0.98654-1.data-00000-of-00001
wget https://storage.googleapis.com/deepvariant/complete-case-study-testdata/complete-t1plus/2026/checkpoint-141056-0.98654-1.index
wget https://storage.googleapis.com/deepvariant/complete-case-study-testdata/complete-t1plus/2026/example_info.json

# ===================================================================
# PanGenie SV input VCF: HPRC-GRCh38 (88 haplotypes)
# ===================================================================
mkdir -p $DB/pangenie/grch38/reference && cd $DB/pangenie/grch38/reference
wget https://zenodo.org/record/6797328/files/cactus_filtered_ids.vcf.gz && \
gzip -d cactus_filtered_ids.vcf.gz

# ===================================================================
# PanGenie index files for the SV input VCF: HPRC-GRCh38 (88 haplotypes)
# ===================================================================
reference=$DB/fa_path/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
SV_input_vcf=$DB/pangenie/grch38/reference/cactus_filtered_ids.vcf
$pangenie_index -r ${reference} -v ${SV_input_vcf} -t 32 -o HPRC_GRCh38 -e 100000

# ===================================================================
# SNV/INDEL/SV truth set: HG002 T2T Q100 v1.1
# ===================================================================
mkdir -p $DB/HG2_T2TQ100_Truthset/GRCh38 && cd $DB/HG2_T2TQ100_Truthset/GRCh38
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed

# ===================================================================
# RTG format files for the reference genome: hg38/GRCh38
# ===================================================================
mkdir -p $DB/RTG/GRCh38.sdf && cd $DB/RTG
$rtg format -o GRCh38.sdf ${reference}

# ===================================================================
# RTG stratification files for the reference genome: hg38/GRCh38
# ===================================================================
mkdir -p $DB/stratifications/v3.6 && cd $DB/stratifications/v3.6
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v3.6/genome-stratifications-GRCh38@all.tar.gz && \
tar xzf genome-stratifications-GRCh38@all.tar.gz


# ===================================================================
# ExpansionHunter files for the reference genome: hg38/GRCh38
# ===================================================================
# Copy files from the workflow folder
mkdir -p $DB/STR/RepeatCatalogs/hg38 && cd $DB/STR/RepeatCatalogs/hg38
cp $DB/../ STR/RepeatCatalogs/hg38/variant_catalog.json ./