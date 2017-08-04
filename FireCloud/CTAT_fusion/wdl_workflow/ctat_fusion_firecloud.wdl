
task CTAT_FUSION_TASK {

    File input_bam
    File genome_lib_tar_gz
    String sample_name

    command {

    # revert aligned bam
    java -Xmx1000m -jar /usr/local/bin/picard.jar \
        RevertSam \
        INPUT=${input_bam} \
        OUTPUT_BY_READGROUP=false \
        VALIDATION_STRINGENCY=LENIENT \
        ATTRIBUTE_TO_CLEAR=FT \
        SORT_ORDER=queryname \
        OUTPUT=${sample_name}.reverted.bam 


    # bam to fastq
    java -jar /usr/local/bin/picard.jar \
        SamToFastq I=${sample_name}.reverted.bam \
        F=${sample_name}_1.fastq F2=${sample_name}_2.fastq \
        INTERLEAVE=false NON_PF=true \
        CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2

    # fusion pipe
    /usr/local/bin/fusion_pipe_runner.py \
            --left_fq ${sample_name}_1.fastq \
            --right_fq ${sample_name}_2.fastq \
            --genome_lib_tar_gz ${genome_lib_tar_gz} \
            --output ${sample_name} 

    mv ${sample_name}.tar.gz  ${sample_name}.ctat_fusion.tar.gz


    }
    
    output {
      File output_tar_gz="${sample_name}.ctat_fusion.tar.gz"
    }

    runtime {
            docker: "trinityctat/firecloud_ctatfusion:0.0.1"
            disks: "local-disk 200 SSD"
            memory: "50G"
            cpu: "16"
    }


}

workflow ctat_fusion_wf {

    String sample_name
    File rnaseq_aligned_bam
    File genome_lib_tar_gz

    call CTAT_FUSION_TASK {
        input: input_bam=rnaseq_aligned_bam,
	       sample_name=sample_name,
               genome_lib_tar_gz=genome_lib_tar_gz
    }

    

}


