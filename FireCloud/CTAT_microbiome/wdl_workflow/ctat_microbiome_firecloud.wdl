workflow ctat_microbiome_workflow {

     String sample_name
     String refindex
     Int threads
     Int primary_assignment
     File left_fq
     File right_fq

     call centrifuge {
       input:
          sample_name=sample_name,
          refindex=refindex,
          threads=threads,
          primary_assignment=primary_assignment,
          left_fq=left_fq,
          right_fq=right_fq
     } 
}


task centrifuge {

    String sample_name
    String refindex
    Int threads
    Int primary_assignment
    File left_fq
    File right_fq

   command {
       centrifuge \
          -p ${threads} \
          -k ${primary_assignment} \
          -q -x ${refindex} \
          -1 ${left_fq} \
          -2 ${right_fq} \
          -S ${sample_name}.centrifuge.results.txt \
          --report-file ${sample_name}.centrifuge.report.txt
    }
	
   output {
       File classificationReport="${sample_name}.centrifuge.report.txt"
       File classificationResults="${sample_name}.centrifuge.results.txt"
   }

}

