
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
            docker: "trinityctat/firecloud_ctatfusion:latest"
            disks: "local-disk 100 SSD"
            memory: "50G"
            cpu: "16"
    }


}

workflow ctat_fusion_wf {

    String sample_name
    File input_left_fq_gz
    File input_right_fq_gz
    File genome_lib_tar_gz
    

    call CTAT_FUSION_TASK {
        input: left_fq_gz=input_left_fq_gz,
               right_fq_gz=input_right_fq_gz,
               sample_name=sample_name,
               genome_lib_tar_gz=genome_lib_tar_gz
    }


}


