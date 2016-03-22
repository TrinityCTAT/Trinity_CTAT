
task STAR_FUSION_TASK {

	File left_fq_gz
	File right_fq_gz
	File genome_lib_dir

	command {
			
			/usr/local/bin/CTAT_fusion_wrapper.pl \
            ${left_fq_gz} \
		    ${right_fq_gz} \
	        ${genome_lib_dir} FusionInspector DISCASM

	}
	
	output {

	    File output_ret = "ctat_out.tar.gz"
	
		#String myOutput = read_string(stdout())
 	}

	runtime {
	   docker: "trinityctat/ctatfusion:latest"
 	   disks: "local-disk 100 SSD"
	   memory: "45G"
	   cpu: "4"
    }


}


workflow CTAT_FUSION_WORKFOW {

    File input_left_fq_gz
	File input_right_fq_gz
	File genome_lib_dir

	call STAR_FUSION_TASK {  input: left_fq_gz=input_left_fq_gz,
					  		  	    right_fq_gz=input_right_fq_gz,
									genome_lib_dir=genome_lib_dir }

	

}
