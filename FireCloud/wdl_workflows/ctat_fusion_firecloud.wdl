
task PicardRevertSam {
  File input_bam
  String sample_name  

  command {
    java -Xmx1000m -jar /usr/local/bin/picard.jar \
    RevertSam \
    INPUT=${input_bam} \
    OUTPUT_BY_READGROUP=false \
    VALIDATION_STRINGENCY=LENIENT \
    ATTRIBUTE_TO_CLEAR=FT \
    SORT_ORDER=queryname \
    OUTPUT=${sample_name}.reverted.bam 
  }
  output {
    File unmapped_bam = "${sample_name}.reverted.bam"
  }

  runtime {
        docker: "trinityctat/firecloud_ctatfusion:0.0.1"
        disks: "local-disk 100 HDD"
        memory: "1200 MB"
    }
}


task PicardSamToFastq {
  File input_unmapped_bam
  String sample_name

  command {
    java -jar /usr/local/bin/picard.jar \
    SamToFastq I=${input_unmapped_bam} \
    F=${sample_name}_1.fastq F2=${sample_name}_2.fastq \
    INTERLEAVE=false NON_PF=true \
    CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2
  
    gzip ${sample_name}_1.fastq ${sample_name}_2.fastq

  }
  output {
	File left_fq_gz = "${sample_name}_1.fastq.gz"
        File right_fq_gz = "${sample_name}_2.fastq.gz"
  }

  runtime {
        docker: "trinityctat/firecloud_ctatfusion:0.0.1"
        disks: "local-disk 100 HDD"
        memory: "2G"
    }
}


task CTAT_FUSION_TASK {

    File left_fq_gz
    File right_fq_gz
    File genome_lib_tar_gz
    String sample_name

    command {
        /usr/local/bin/fusion_pipe_runner.py \
            --left_fq ${left_fq_gz} \
            --right_fq  ${right_fq_gz} \
            --genome_lib_tar_gz ${genome_lib_tar_gz} \
            --output ${sample_name} 
    }
    
    output {
      File output_tar_gz="${sample_name}.tar.gz"
    }

    runtime {
            docker: "trinityctat/firecloud_ctatfusion:0.0.1"
            disks: "local-disk 100 SSD"
            memory: "50G"
            cpu: "16"
    }


}

workflow ctat_fusion_wf {

    String sample_name
    File rnaseq_aligned_bam
    File genome_lib_tar_gz


    call PicardRevertSam {
        input: input_bam=rnaseq_aligned_bam,
	       sample_name=sample_name
    }

    call PicardSamToFastq {
        input: input_unmapped_bam=PicardRevertSam.unmapped_bam,
	       sample_name=sample_name
    }
	

    call CTAT_FUSION_TASK {
        input: left_fq_gz=PicardSamToFastq.left_fq_gz,
               right_fq_gz=PicardSamToFastq.right_fq_gz,
               sample_name=sample_name,
               genome_lib_tar_gz=genome_lib_tar_gz
    }

    

}


