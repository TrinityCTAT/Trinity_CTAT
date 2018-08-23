
task CTAT_FUSION_TASK_BAM {

    File genome_lib_tar_gz
    String sample_name
    File input_bam
    
    command {

    set -e

    # initial potential cleanup of read names in the bam file
    /usr/local/bin/sam_readname_cleaner.py ${input_bam} ${input_bam}.cleaned.bam


    # revert aligned bam
    java -Xmx1000m -jar /usr/local/bin/picard.jar \
        RevertSam \
        INPUT=${input_bam}.cleaned.bam \
        OUTPUT_BY_READGROUP=false \
        VALIDATION_STRINGENCY=SILENT \
        SORT_ORDER=queryname \
        OUTPUT=${sample_name}.reverted.bam 


    # bam to fastq
    java -jar /usr/local/bin/picard.jar \
        SamToFastq I=${sample_name}.reverted.bam \
        F=${sample_name}_1.fastq F2=${sample_name}_2.fastq \
        INTERLEAVE=false NON_PF=true \
        CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2

    # untar the genome lib
    tar xvf ${genome_lib_tar_gz}

    # starFusion

    /usr/local/bin/STAR-Fusion \
         --left_fq ${sample_name}_1.fastq \
         --right_fq ${sample_name}_2.fastq \
     --CPU 16 \
         --genome_lib_dir $(echo $(basename ${genome_lib_tar_gz}) | sed s/.plug-n-play.tar.gz//)/ctat_genome_lib_build_dir \
         --output_dir ${sample_name} --run_STAR_only
     
    cp ${sample_name}/std.Chimeric.out.junction ${sample_name}.Chimeric.out.junction

    gzip ${sample_name}.Chimeric.out.junction

    }
    
    output {
      File star_fusion_results="${sample_name}.Chimeric.out.junction.gz"
    }

    runtime {
            docker: "trinityctat/firecloud_ctatfusion:0.0.3"
            disks: "local-disk 200 SSD"
            memory: "50G"
            cpu: "16"
    }


}



task CTAT_FUSION_TASK_FASTQ {

    File genome_lib_tar_gz
    String sample_name
    File left_fq
    File right_fq
    
    command {

    set -e

    # untar the genome lib
    tar xvf ${genome_lib_tar_gz}

    # starFusion

    /usr/local/bin/STAR-Fusion \
         --left_fq ${left_fq} \
         --right_fq ${right_fq} \
     --CPU 16 \
         --genome_lib_dir $(echo $(basename ${genome_lib_tar_gz}) | sed s/.plug-n-play.tar.gz//)/ctat_genome_lib_build_dir \
         --output_dir ${sample_name} --run_STAR_only
     
    cp ${sample_name}/std.Chimeric.out.junction ${sample_name}.Chimeric.out.junction

    gzip ${sample_name}.Chimeric.out.junction

    }
    
    output {
      File star_fusion_results="${sample_name}.Chimeric.out.junction.gz"
    }


    runtime {
            docker: "trinityctat/firecloud_ctatfusion:0.0.3"
            disks: "local-disk 200 SSD"
            memory: "50G"
            cpu: "16"
    }


}


task CTAT_FUSION_TASK_FQPAIRTARGZ {

    File genome_lib_tar_gz
    String sample_name
    File fastq_pair_tar_gz
    
    command {

    set -e
    
    # untar the fq pair
    tar xvf ${fastq_pair_tar_gz}

    left_fq="*_1.fastq"
    right_fq="*_2.fastq"

    if [ -z $left_fq ] || [ -z $right_fq ]; then
        ls -l *
        echo "ERROR, could not identify left and right fastq files"
        exit 1
    fi

    
    # untar the genome lib
    tar xvf ${genome_lib_tar_gz}

    # starFusion

    /usr/local/bin/STAR-Fusion \
         --left_fq $left_fq \
         --right_fq $right_fq \
     --CPU 16 \
         --genome_lib_dir $(echo $(basename ${genome_lib_tar_gz}) | sed s/.plug-n-play.tar.gz//)/ctat_genome_lib_build_dir \
         --output_dir ${sample_name}  --run_STAR_only
     
    cp ${sample_name}/std.Chimeric.out.junction ${sample_name}.Chimeric.out.junction

    gzip ${sample_name}.Chimeric.out.junction

    }
    
    output {
      File star_fusion_results="${sample_name}.Chimeric.out.junction.gz"
    }
     
    runtime {
            docker: "trinityctat/firecloud_ctatfusion:0.0.3"
            disks: "local-disk 200 SSD"
            memory: "50G"
            cpu: "16"
    }


}




workflow ctat_fusion_wf {

    String sample_name
    File genome_lib_tar_gz
    File? rnaseq_aligned_bam
    File? left_fq
    File? right_fq
    File? fastq_pair_tar_gz

    if (defined(rnaseq_aligned_bam)) {
        call CTAT_FUSION_TASK_BAM {
            input:
              input_bam=rnaseq_aligned_bam,
                 sample_name=sample_name,
                 genome_lib_tar_gz=genome_lib_tar_gz
        }
    }

    if (defined(left_fq)) {
        call CTAT_FUSION_TASK_FASTQ {
               input:
              sample_name=sample_name,
              genome_lib_tar_gz=genome_lib_tar_gz,
                 left_fq=left_fq,
              right_fq=right_fq
        }
    }

    if (defined(fastq_pair_tar_gz)) {
        call CTAT_FUSION_TASK_FQPAIRTARGZ {
            input:
              sample_name=sample_name,
              genome_lib_tar_gz=genome_lib_tar_gz,
                 fastq_pair_tar_gz=fastq_pair_tar_gz
        }

    }
    
}


