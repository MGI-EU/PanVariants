import java.text.SimpleDateFormat

def timeStamp(tip) {
    def now = new Date()
    def time_stamp = now.format("yyyy-MM-dd HH:mm:ss")
    println("$time_stamp: $tip.")
}


// read input data list
def readInput(input) {
    timeStamp("read input: $input")
    def input_line_num = 0
    def input_list = []
    file(input).withReader {
        String input_line
        while( input_line = it.readLine() ){
            if (!input_line.trim() || input_line.trim().startsWith('#')) {
                continue
            }
            
            input_line_num += 1
            input_line_split = input_line.split(' +|\t')
            
            if (input_line_split.size() < 3) {
                timeStamp("Invalid sample line format: $input_line. Expected: sampleID fastq1 fastq2")
                exit 1
            }
            
            timeStamp("check ${input_line_num} sample ${input_line_split[0]}")

            // Check FASTQ files
            def sampleId = input_line_split[0]
            def fastq1 = input_line_split[1]
            def fastq2 = input_line_split[2]

            input_list.add([sampleId, fastq1, fastq2])
        }
    }
    return input_list
}

// Check if file exists
def checkFile(file_path, miss_file_num) {
    if (!new File(file_path).exists()) {
        timeStamp("ERROR: File $file_path does not exist")
        miss_file_num += 1
    }
    return miss_file_num
}

process merge_fq_PE {
    input:
        val(sample_info)
    output:
        tuple val(sample_name),path("${sample_name}.merge.fq.1.gz"),path("${sample_name}.merge.fq.2.gz"), emit: merge_fq
    script:
        sample_name = sample_info[0]
        fq1path = sample_info[1]
        fq2path = sample_info[2]
        if (fq2path == null ) {
            allfqpath = "${fq1path}"
        } else {
            allfqpath = "${fq1path} ${fq2path}"
        }
        allfqpath = allfqpath.replace(';', ' ')
        allfqpath_split = allfqpath.split(' ')
        allfqpath_split.each { filepath ->
        if (!file(filepath).exists()) {
            error "There is no file: ${filepath}, please check the path: ${filepath}"
            System.exit(1)
            }
        }
        if (fq1path.contains(';')){
            allfq1path = fq1path.replace(';', ' ')
            allfq2path = fq2path.replace(';', ' ')
        } else {
            allfq1path = ""
            allfq2path = ""
        }
        """
        if [[ "$fq1path" =~ ";" ]]; then
            cat ${allfq1path} > ${sample_name}.merge.fq.1.gz &&
            cat ${allfq2path} > ${sample_name}.merge.fq.2.gz
        else
            ln -sf ${fq1path} ${sample_name}.merge.fq.1.gz &&
            ln -sf ${fq2path} ${sample_name}.merge.fq.2.gz
        fi
        """   
}

process fqfilter{
    input:
        tuple val(sample_name),path(read1),path(read2)
    
    output:
        tuple val(sample_name),path("${sample_name}/${sample_name}_clean_1.fq.gz"), path("${sample_name}/${sample_name}_clean_2.fq.gz"), emit: clean_fq
        tuple val(sample_name),path("${sample_name}/Basic_Statistics_of_Sequencing_Quality.txt"), emit: base_quality
        tuple val(sample_name),path("${sample_name}/Statistics_of_Filtered_Reads.txt"), emit: statistics_of_filtered_reads
        tuple val(sample_name),path("${sample_name}/*"), emit: other_out
   
    script:
        """
        ${params.soapnuke} filter \
         -n ${params.n_base_ratio} \
         -l ${params.qcthreshold} \
         -q ${params.low_quality_ratio} \
         -T ${task.cpus} \
         -f ${params.adapter1} \
         -r ${params.adapter2} \
         -1 ${read1} \
         -2 ${read2} \
         -C ${sample_name}_clean_1.fq.gz \
         -D ${sample_name}_clean_2.fq.gz \
         -o ${sample_name} && \
        read1_path=`realpath ${read1}` && \
        work_dir="\$(realpath ${params.outdir}/../ 2>/dev/null || pwd)" && \
        if [ "\$read1_path" = "\$work_dir"/* ]; then
            rm -f `realpath ${read1}` `realpath ${read2}`
        fi
        """
}