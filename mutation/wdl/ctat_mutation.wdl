workflow ctat_rnaseq_mutation {
    call single_paired_mutation
}

task single_paired_mutation {
    File cosmic_vcf
    File cosmic_vcf_index
    File cravat_header
    File darned
    File dbsnp_vcf
    File dbsnp_gz
    File dbsnp_gz_index
    File radar
    File reference
    File reference_bed
    File reference_fai
    File reference_dict
    File sample_left
    File sample_right
    File index
    String email
    String tissue_type

    command {
        rnaseq_mutation_pipeline.py \
                --left ${sample_left} \
                --right ${sample_right} \
                --tissue_type ${tissue_type} \
                --email ${email} \
                --index ${index} \
                --plot \
                --cosmic_vcf ${cosmic_vcf} \
                --reference ${reference} \
                --vcf ${dbsnp_vcf} \
                --darned ${darned} \
                --radar ${radar} \
                --bed ${reference_bed} \
                --cravat_annotation_header ${cravat_header} \
                --is_hg19 \
                --wdl_compatible_run \
                --update "GenomeAnalysisTK.jar:/usr/local/bin,AddOrReplaceReadGroups.jar:/usr/local/bin,MarkDuplicates.jar:/usr/local/bin,SortSam.jar:/usr/local/bin,java -:/usr/lib/jvm/java-7-oracle/bin,snpEff.jar:/usr/local/bin/snpEff"
    }
    output {
        File calling_bam = "misc/recalibrated.bam"
        File raw_variants = "variants.vcf"
        File filtered_variants = "variants_initial_filtering_clean_snp_RNAedit.vcf.gz"
        File filtered_annotated_variants = "variants_initial_filtering_clean_snp_RNAedit.vcf_dbsnp.vcf_snpeff_updated.vcf.gz"
        File final_variants = "cancer.vcf"
        File final_variants_tab = "cancer.tab"
    }
    runtime {
        docker: "ttickle/ctat_mutation:1"
    }
    parameter_meta {
        cosmic_vcf: "Bgzipped Cosmic VCF http://cancer.sanger.ac.uk/cosmic/analyses"
        cosmic_vcf_index: ".csi index of the bgzipped cosmic vcf file"
        cravat_header: "Technical file used to normalize CRAVAT headers"
        darned: "VCF of Known RNAEditing events. http://rnaedit.com"
        dbsnp_vcf: "DBSNP variant VCF used for annotation. https://www.ncbi.nlm.nih.gov/projects/SNP"
        dbsnp_gz: "bgzipped dbsnp.vcf"
        dbsnp_gz_index: ".csi indexed bgzipped dbnp.vcf"
        radar: "VCF of known RNAEditing events. http://darned.ucc.ie"
        reference: "Reference genome fasta"
        reference_bed: "Reference genome BED"
        reference_fai: ".fai index of reference genome"
        reference_dict: ".dict index of reference genome"
        sample_left: "Left paired RNA-Seq sample"
        sample_right: "Right Paired RNA-Seq sample"
        index: "Star Aligner index of Reference genome"
        email: "Your email, if CRAVAT has an error your email can be used to troubleshoot."
        tissue_type: "Used in CRAVAT analysis, any tissue type name found in the left column in the Analysis section of https://www.cravat.us/CRAVAT/help.jsp"
    }
}
