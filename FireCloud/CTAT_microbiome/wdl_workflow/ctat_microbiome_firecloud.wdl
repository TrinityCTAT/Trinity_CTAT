
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
        docker: "trinityctat/firecloud_ctatmicrobiome:0.0.1"
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
        docker: "trinityctat/firecloud_ctatmicrobiome:0.0.1"
        disks: "local-disk 100 HDD"
        memory: "2G"
    }
}






workflow ctat_microbiome_workflow {

     String sample_name
     File rnaseq_aligned_bam
     File refindex_tar_gz
     Int threads
     Int primary_assignment


   call PicardRevertSam {
        input: input_bam=rnaseq_aligned_bam,
	       sample_name=sample_name
    }

    call PicardSamToFastq {
        input: input_unmapped_bam=PicardRevertSam.unmapped_bam,
	       sample_name=sample_name
    }


    call centrifuge {
       input:
          sample_name=sample_name,
          refindex_tar_gz=refindex_tar_gz,
          threads=threads,
          primary_assignment=primary_assignment,
          left_fq=PicardSamToFastq.left_fq_gz,
          right_fq=PicardSamToFastq.right_fq_gz
     } 
}


task centrifuge {

    String sample_name
    File refindex_tar_gz
    String refindex_name
    Int threads
    Int primary_assignment
    File left_fq
    File right_fq

   command {

     tar xvf ${refindex_tar_gz}

     centrifuge \
          -p ${threads} \
          -k ${primary_assignment} \
          -q -x ${refindex_name} \
          -1 ${left_fq} \
          -2 ${right_fq} \
          -S ${sample_name}.centrifuge.results.txt \
          --report-file ${sample_name}.centrifuge.report.txt
    }
	
   output {
       File classificationReport="${sample_name}.centrifuge.report.txt"
       File classificationResults="${sample_name}.centrifuge.results.txt"
   }

   runtime {
            docker: "trinityctat/firecloud_ctatfusion:0.0.1"
            disks: "local-disk 100 SSD"
            memory: "10G"
            cpu: "4"
    }



}

