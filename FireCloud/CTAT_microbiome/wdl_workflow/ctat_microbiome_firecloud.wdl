
workflow ctat_microbiome_workflow {

     
     String sample_name
     File rnaseq_aligned_bam
     File refindex_tar
     String refindex_name
     Int threads
     Int primary_assignment

    call centrifuge {
       input:
          sample_name=sample_name,
	  rnaseq_aligned_bam=rnaseq_aligned_bam,
          refindex_tar=refindex_tar,
          refindex_name=refindex_name,
          threads=threads,
          primary_assignment=primary_assignment,

     } 
}


task centrifuge {

    String sample_name
    File rnaseq_aligned_bam
    File refindex_tar
    String refindex_name
    Int threads
    Int primary_assignment


   command {

    set -e

    # revert bam file
    java -Xmx1000m -jar /usr/local/bin/picard.jar \
        RevertSam \
        INPUT=${rnaseq_aligned_bam} \
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


     # unpack ref index
     tar xvf ${refindex_tar}

     # run centrifuge
     centrifuge \
          -p ${threads} \
          -k ${primary_assignment} \
          -q -x ${refindex_name} \
          -1 ${sample_name}_1.fastq \
          -2 ${sample_name}_2.fastq \
          -S ${sample_name}.centrifuge.results.txt \
          --report-file ${sample_name}.centrifuge.report.txt

     tar -zcvf ${sample_name}.centrifuge.tar.gz ${sample_name}.centrifuge.results.txt ${sample_name}.centrifuge.report.txt

    }
	
   output {
       File centrifuge_results_tar_gz="${sample_name}.centrifuge.tar.gz"
   }

   runtime {
            docker: "trinityctat/firecloud_ctatmicrobiome:0.0.1"
            disks: "local-disk 200 SSD"
            memory: "10G"
            cpu: "4"
    }



}

