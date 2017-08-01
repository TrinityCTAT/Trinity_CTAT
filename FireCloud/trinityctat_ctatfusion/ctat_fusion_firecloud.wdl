
task UNTAR {

    File tar_gz_file

	command {
       tar xvf ${tar_gz_file}
    }
}


task CTAT_FUSION_TASK {

    File left_fq_gz
    File right_fq_gz
    String genome_lib_dir
    String output_dir_name

    command {
        /usr/local/src/STAR-Fusion-v1.0.0/STAR-Fusion \
            --left_fq ${left_fq_gz} \
            --right_fq  ${right_fq_gz} \
            --genome_lib_dir ${genome_lib_dir} \
            --extract_fusion_reads \
            --FusionInspector validate \
            --denovo_reconstruct \
            --annotate \
            --examine_coding_effect \
            --output_dir ${output_dir_name} 
    }
	
    runtime {
            docker: "trinityctat/ctatfusion:latest"
            disks: "local-disk 100 SSD"
            memory: "50G"
            cpu: "16"
    }
}


task CAPTURE_OUTPUTS {

    String dirname

    command {

        tar -cvfz ${dirname}.tar.gz \
           ${dirname}/star-fusion.fusion_predictions.abridged.annotated.coding_effect.tsv \
           ${dirname}/star-fusion.fusion_evidence_reads_1.fq \
           ${dirname}/star-fusion.fusion_evidence_reads_2.fq \
           ${dirname}/FusionInspector/finspector.spanning_reads.bam \
           ${dirname}/FusionInspector/finspector.junction_reads.bam \
           ${dirname}/FusionInspector/finspector.fusion_predictions.final.abridged.FFPM \
           ${dirname}/FusionInspector/finspector.gmap_trinity_GG.fusions.gff3 \
           ${dirname}/FusionInspector/finspector.gmap_trinity_GG.fusions.fasta
           
                   
    }

    output {
       File out_tar_gz = ${dirname}.tar.gz
    }
}


workflow ctat_fusion_wf {

    String sample_name
    File input_left_fq_gz
    File input_right_fq_gz
    File genome_lib_tar_gz
	String genome_lib_dirname


	call UNTAR {
        input: tar_gz_file=genome_lib_tar_gz
    }


    call CTAT_FUSION_TASK {
        input: left_fq_gz=input_left_fq_gz,
               right_fq_gz=input_right_fq_gz,
               output_dir_name=sample_name,
               genome_lib_dirname=genome_lib_dirname
	}

	call CAPTURE_OUTPUTS {

	    input: dirname=${sample_name}

    }

}


