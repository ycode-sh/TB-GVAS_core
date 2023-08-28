// 1. Trimming with trimmomatic
// Inputs: reads, adapters; output:trimmed paired and unpaired reads

//reads = "reads"
//reads_pair_id = "pair_id"

if(params.in_data_type == "se_illumina_reads" | params.in_data_type == "minion_ont_reads"){
    pe_reads = "reads"
    reads_pair_id = "pair_id"
} 


process trim_fastq {    
    publishDir "${params.out_dir_int}/Trimmed_reads"
    
    input:
    val(trim_sample_script)
    path (trim_fastq_se_reads, stageAs: params.stage == "se_reads" ? "se_reads": null)
    params.what_input == "path_only" ? 'path(trim_fastq_se_reads)': 'tuple val(reads_pair_id), path(pe_reads)' // Dynamic input declaration
    val(adapter)
    val(input_read_type) 
   

    output:
    path "sample*.fastq", emit: trimmed_reads
    

    script:
    
    """
    if [[ ${params.go} == 0 ]]; then
        bash ${trim_sample_script} ${trim_fastq_se_reads} ${adapter} ${input_read_type}
    elif [[ ${params.go} == 1 ]]; then
        bash ${trim_sample_script} ${pe_reads[0]} ${pe_reads[1]} ${reads_pair_id} ${adapter} ${input_read_type}
    fi

    """
}


// 2. Kraken analysis
// Input: reads (raw or trimmed); outputs: Classified and unclassified reads, kraken reports and outputs

process kraken_contamination {
    publishDir "${params.out_dir_int}/Kraken_outputs"
    cpus = 8
    memory = 30.GB

    errorStrategy 'ignore'


    input:
    val (kraken_script)
    path (kc_se_reads, stageAs: params.stage == "se_reads" ? "se_reads": null)
    params.what_input == "path_only" ? 'path(kc_se_reads)': 'tuple val(reads_pair_id), path(pe_reads)' // Dynamic input declaration
    val (kraken_database_path)
    val (input_read_type)
    
    output:
        path "classified_*.fastq", emit: classified_fastq
        path "unclassified_*.fastq", emit: unclassified_fastq 
        path "*_kraken_*", emit: kraken_reports
        
        
    script:
        """
        if [[ ${params.go} == 0 ]]; then
            bash ${kraken_script} ${kc_se_reads} ${kraken_database_path} ${input_read_type}
        elif [[ ${params.go} == 1 ]]; then
            bash ${kraken_script} ${pe_reads[0]} ${pe_reads[1]} ${reads_pair_id} ${kraken_database_path} ${input_read_type}
        fi
        """

}


// 3. 16S Ribosomal RNA Analysis
// Input: Fasta contig; Output: Blast output
process find_16s_hit {
    publishDir "${params.out_dir_int}/16S_outputs"
    errorStrategy 'ignore'


    input:
        val find_16S_hits_script
        path fasta_contigs
        val (txdb_path_str) 


    output: 
        path "*_16s_output.tsv", emit: rRNA_16S_result


    script:
        """

        bash ${find_16S_hits_script} ${fasta_contigs} ${txdb_path_str}

        """


}

// 4. De novo genome assembly
// With flye or spades
// Inputs: raw or trimmed reads

process genome_assembly {
    publishDir "${params.out_dir_int}/Contigs"
    errorStrategy 'ignore'
    
    input:
        val (genome_assembly_script)
        path (ga_se_reads, stageAs: params.stage == "se_reads" ? "se_reads": null)
        params.what_input == "path_only" ? 'path(ga_se_reads)': 'tuple val(reads_pair_id), path(pe_reads)' // Dynamic input declaration
        val (input_read_type)

    output:
        path "*_contig.fasta"   

    script:
        """
        if [[ ${params.go} == 0 ]]; then
            bash ${genome_assembly_script}  ${ga_se_reads} ${input_read_type}
        elif [[ ${params.go} == 1 ]]; then
            bash ${genome_assembly_script} ${pe_reads[0]} ${pe_reads[1]} ${reads_pair_id} ${input_read_type}
        fi
        """

}


// 5. Map reads to Reference
// Emit sam with either bwa or minimap
process emit_sam {
    publishDir "${params.out_dir_int}/Sam_files"
    input:
        val (emit_sam_script)
        path (emit_sam_se_reads, stageAs: params.stage == "se_reads" ? "se_reads": null)
        params.what_input == "path_only" ? 'path(emit_sam_se_reads)': 'tuple val(reads_pair_id), path(pe_reads)' // Dynamic input declaration
        val (fastafile)
        val (input_read_type)
        


    output:
        path "*_sam", emit: reads_sam

    script:
        """
        if [[ ${params.go} == 0 ]]; then
            bash ${emit_sam_script} ${emit_sam_se_reads} ${fastafile} ${input_read_type}
        elif [[ ${params.go} == 1 ]]; then
            bash ${emit_sam_script} ${pe_reads[0]} ${pe_reads[1]} ${fastafile} ${reads_pair_id} ${input_read_type}
        fi
        """
}



// Emit bam
process coordsort_sam {
        publishDir "${params.out_dir_int}/Bam_files"

        input:
            val (coordsort_sam_script)
            path reads_sam



        output:
            path "*_sam.bam", emit: reads_bam
            

        script:
            """

            bash ${coordsort_sam_script} ${reads_sam}

            """
}

// Reads decontamination

process bamtofastq {
    publishDir "${params.out_dir_int}/Decontaminated_reads"
    input:
        val (bamtofastq_script)
        path bamfile
        val (input_read_type)

    output:
        path "*fastq"

    script:
        """

        bash ${bamtofastq_script} ${bamfile} ${input_read_type}

        """

}

// File parsing processes

process parse_16S {
    publishDir "${params.out_dir_int}/Parsed_16S_outputs"
    input:
        val (parse_16S_script)
        path output_16S_files
        val (perc_cov)
        val (perc_id) 

    output:
        path "parsed_16S_samples_output.tsv", emit: parsed_16S_out

    script:
        """

        python ${parse_16S_script} ${output_16S_files} ${perc_cov} ${perc_id}

        """

}

process parse_kraken {
    publishDir "${params.out_dir_int}/Parsed_kraken_outputs"
    input:
        val (parse_kraken_script)
        path output_kraken_files


    output:
        path "parsed_kraken_report.tsv", emit: parsed_kraken_report

    script:
        """

        python ${parse_kraken_script} ${output_kraken_files}

        """


}

