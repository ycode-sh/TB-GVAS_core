

// 1. Trimming with trimmomatic
// Inputs: reads, adapters; output:trimmed paired and unpaired reads
process trim_fastq {    

    input:
    val(trim_sample_script)
    tuple val(trim_fastq_reads_pair_id), path(trim_fastq_reads)
    val(adapter) 
   

    output:
    path "${trim_fastq_reads_pair_id}*.fastq", emit: trim_out

    script:
    
    """
    bash ${trim_sample_script} ${trim_fastq_reads[0]} ${trim_fastq_reads[1]} ${trim_fastq_reads_pair_id} ${adapter}
    """
}


