workflow ctat_rnaseq_mutation {
    File in_reference
    File in_alignment_index
    File in_cosmic_vcf
    File in_cravate_header
    File in_dbsnp_vcf
    File in_left
    File in_right
    String in_email
    String in_misc = "misc"
    String in_name
    String in_sample_type

    String in_gatk_path
    String in_picard_path
    String in_snpeff_path

    call align {
        input: sample_name=in_name,
            align_index=in_alignment_index,
            align_left_sample=in_left,
            align_right_sample=in_right,
            align_misc_dir=in_misc
    }
    call index_aligned_bam {
        input: sample_name=in_name,
            in_bam=align.aligned_bam,
            index_misc_dir=in_misc
    }
    call add_rg {
        input: misc=in_misc,
            path=in_picard_path,
            sample_name=in_name,
            in_bam=align.aligned_bam
    }
    call mark_duplicates {
        input: misc=in_misc,
            path=in_picard_path,
            sample_name=in_name,
            in_bam=add_rg.out_bam
    }
    call gatk_split_cigar_reads {
        input: misc=in_misc,
            path=in_gatk_path,
            sample_name=in_name,
            reference=in_reference,
            in_bam=mark_duplicates.out_bam
    }
    call gatk_realigner_target_creator {
        input: misc=in_misc,
            path=in_gatk_path,
            sample_name=in_name,
            reference=in_reference,
            in_bam=gatk_split_cigar_reads.out_bam,
             known_sites_db=in_dbsnp_vcf
    }
    call gatk_indel_realigner {
        input: misc=in_misc,
            path=in_gatk_path,
            sample_name=in_name,
            reference=in_reference,
            in_bam=gatk_split_cigar_reads.out_bam,
            intervals=gatk_realigner_target_creator.intervals,
            known_sites_db=in_dbsnp_vcf
    }
    call gatk_base_recalibrator {
        input: misc=in_misc,
            path=in_gatk_path,
            sample_name=in_name,
            realigned_bam=gatk_indel_realigner.out_bam,
            reference=in_reference,
            known_sites_db=in_dbsnp_vcf
    }
    call gatk_recalibrated_bam {
        input: sample_name=in_name,
            path=in_gatk_path,
            reference=in_reference,
            realigned_bam=gatk_indel_realigner.out_bam,
            recalibration_table=gatk_base_recalibrator.out_table
    }
    call gatk_haplotype_caller {
        input: sample_name=in_name,
            path=in_gatk_path,
            reference=in_reference,
            bam=gatk_recalibrated_bam.out_bam
    }
    call gatk_initial_filtering {
        input: sample_name=in_name,
            path=in_gatk_path,
            reference=in_reference,
            to_update_vcf=gatk_haplotype_caller.called_vcf
    }
    call groom_initial_vcf {
        input: sample_name=in_name,
            to_update_vcf=gatk_initial_filtering.updated_vcf
    }
    call subset_to_snps {
        input: sample_name=in_name,
            misc_dir=in_misc,
            to_update_vcf=groom_initial_vcf.updated_vcf
    }
    call edit_rna_editing {
        input: darned_db=in_darned_db,
            radar_db=in_radar_db,
            to_update_vcf=subset_to_snps.updated_vcf
    }
    call compress_bgzip as compress_editing_vcf {
        input: to_compress=edit_rna_editing.updated_vcf
    }
    call index_vcf as index_editing_vcf {
        input: in_vcf=compress_editing_vcf.compressed
    }
    call dbsnp_annotate {
        input: misc=in_misc,
            sample_name=in_name,
            dbsnp_vcf=in_dbsnp_vcf,
            to_update_vcf=index_editing_vcf.out_vcf
    }
    call snpeff {
        input: misc=in_misc,
            path=in_snpeff_path,
            sample_name=in_name,
            to_update_vcf=dbsnp_annotate.updated_vcf
    }
    call update_snpeff_annotations {
        input: misc=in_misc,
            sample_name=in_name,
            to_update_vcf=snpeff.updated_vcf
    }
    call compress_bgzip as compress_snpeff {
        input: to_compress=update_snpeff_annotations.updated_vcf
    }
    call index_vcf as index_snpeff {
        input: in_vcf=compress_snpeff.compressed
    }
    call cosmic_annotate {
        input: misc=in_misc,
            sample_name=in_name,
            cosmic_vcf=in_cosmic_vcf,
            to_annotate_vcf=index_snpeff.out_vcf
    }
    call filter_vcf_cancer {
        input: misc=in_misc,
            sample_name=in_name,
            in_vcf=cosmic_annotate.cosmic_annotated
    }
    call cravat_annotate_get {
        input: email=in_email,
            misc=in_misc,
            sample_name=in_name,
            sample_type=in_sample_type,
            variants_vcf=filter_vcf_cancer.out_vcf
    }
    call cravat_unzip_response {
        input: misc=in_misc,
            sample_name=in_name,
            to_unzip=cravat_annotate_get.out_annotations
    }
    call cravate_move_response {
        input: cravate_data=cravat_unzip_response.unzipped
    }
    call groom_cravat_coding_results {
        input: misc=in_misc,
            sample_name=in_name,
            in_tab=cravate_move_response.location_coding
    }
    call groom_cravat_noncoding_results {
        input: misc=in_misc,
            sample_name=in_name,
            in_tab=cravate_move_response.location_noncoding
    }
    call compress_bgzip as compress_coding_tab {
        input: to_compress=groom_cravat_coding_results.out_tab
    }
    call index_tab as index_coding_tab {
        input: in_tab=compress_coding_tab.compressed
    }
    call compress_bgzip as compress_noncoding_tab {
        input: to_compress=groom_cravat_noncoding_results.out_tab
    }
    call index_tab as index_noncoding_tab {
        input: in_tab=compress_noncoding_tab.compressed
    }
    call cravat_annotate_add_coding {
        input: misc=in_misc,
            sample_name=in_name,
            cravate_header=in_cravate_header,
            cravate_annotations=index_coding_tab.out_tabix,
            in_vcf=filter_vcf_cancer.out_vcf
    }
    call cravat_annotate_add_noncoding {
        input: misc=in_misc,
            sample_name=in_name,
            cravate_header=in_cravate_header,
            cravate_annotations=index_noncoding_tab.out_tabix,
            in_vcf=cravat_annotate_add_noncoding.out_vcf
    }
    call cravate_filter {
        input: misc=in_misc,
            sample_name=in_name,
            in_vcf=cravat_annotate_add_noncoding.out_vcf
    }
    call groom_cravate_vcf {
        input: sample_name=in_name,
            to_update_vcf=cravate_filter.out_vcf
    }
    call gatk_variants_table {
        input: sample_name=in_name,
            path=in_snpeff_path,
            reference=in_reference,
            in_vcf=groom_cravate_vcf.updated_vcf
    }
}

# Note assumes GZIPPED files
task align {
    String sample_name
    File align_index
    File align_left_sample
    File align_right_sample
    String align_misc_dir
    command {
        STAR --genomeDir ${align_index} \
        --runThreadN 8 \
        --readFilesCommand "gunzip -c" \
        --readFilesIn ${align_left_sample} \
        ${align_right_sample} \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --limitBAMsortRAM 30000000000 \
        --outFileNamePrefix ${align_misc_dir}
    }
    output {
        File aligned_bam = "${sample_name}/${align_misc_dir}/${sample_name}_aligned.sortedByCoord.out.bam"
    }
}

task compress_bgzip {
    File to_compress
    command {
        bgzip -c ${to_compress} > ${to_compress}.gz
    }
    output {
        File compressed = "${to_compress}.gz"
    }
}

task cosmic_annotate {
    String misc
    String sample_name
    File to_annotate_vcf
    File cosmic_vcf
    command {
        bcftools annotate --output-type z \
        --annotations ${cosmic_vcf} \
        --columns INFO/COSMIC_ID,INFO/TISSUE,INFO/TUMOR,INFO/FATHMM,INFO/SOMATIC \
        --output$ ${sample_name}/${misc}/${sample_name}_cosmic_ann.vcf \
        ${to_annotate_vcf}
    }
    output {
        File cosmic_annotated = "${sample_name}/${misc}/${sample_name}_cosmic_ann.vcf"
    }
}

task cravat_annotate_add_coding {
    String misc
    String sample_name
    File cravate_annotations
    File cravate_header
    File in_vcf
    command {
        bcftools annotate \
        --annotations ${cravate_annotations} \
        -h ${cravate_header} \
        --columns "CHROM,POS,CHASM_PVALUE,CHASM_FDR,VEST_PVALUE,VEST_FDR" \
        --output-type z \
        --output "${sample_name}/${misc}/${sample_name}_cravat_annotated_coding.vcf.gz" \
        ${in_vcf}
    }
    output {
        File out_vcf = "${sample_name}/${misc}/${sample_name}_cravat_annotated_coding.vcf.gz"
    }
}

task cravat_annotate_add_noncoding {
    String misc
    String sample_name
    File cravate_annotations
    File cravate_header
    File in_vcf
    command {
        bcftools annotate \
        --annotations ${cravate_annotations} \
        -h ${cravate_header} \
        --columns "CHROM,POS,CHASM_PVALUE,CHASM_FDR,VEST_PVALUE,VEST_FDR" \
        --output-type z \
        --output ${sample_name}/${misc}/${sample_name}_cravat_annotated_bothcoding.vcf.gz \
        ${in_vcf}
    }
    output {
        File out_vcf = "${sample_name}/${misc}/${sample_name}_cravat_annotated_bothcoding.vcf.gz"
    }
}

task cravat_annotate_get {
    String sample_type
    String email
    String misc
    String sample_name
    File variants_vcf
    command {
        annotate_with_cravat.py \
        --classifier ${sample_type} \
        --email ${email} \
        --max_attempts 180 \
        --wait 60 \
        ${variants_vcf} ${sample_name}/${misc}/${sample_name}_cravate_data.gz
    }
    output {
        File out_annotations = "${sample_name}/${misc}/${sample_name}_cravate_data.gz"
    }
}

task cravate_filter {
    String misc
    String sample_name
    File in_vcf
    command {
        bcftools filter \
        --include "CHASM_PVALUE < 0.3 || VEST_PVALUE < 0.3" \
        --output-type v \
        --output ${sample_name}/${misc}/${sample_name}_annotated_min_filtered.vcf.gz ${in_vcf}
    }
    output {
        File out_vcf = "${sample_name}/${misc}/${sample_name}_annotated_min_filtered.vcf.gz"
    }
}

task cravate_move_response {
    String cravate_data
    command {
        cp ${cravate_data}/*/Variant.Result.tsv ${cravate_data}
        cp ${cravate_data}/*/Variant_Non-coding.Result.tsv ${cravate_data}
    }
    output {
        File location_coding = "${cravate_data}/Variant.Result.tsv"
        File location_noncoding = "${cravate_data}/Variant_Non-coding.Result.tsv"
    }
}

task cravat_unzip_response {
    String misc
    String sample_name
    File to_unzip
    command {
        unzip -d ${sample_name}/${misc}/${sample_name}_extracted_cravate_data ${to_unzip}
    }
    output {
        File unzipped = "${sample_name}/${misc}/${sample_name}_extracted_cravate_data"
    }
}

task dbsnp_annotate {
    String misc
    String sample_name
    File dbsnp_vcf
    File to_update_vcf
    command {
        bcftools annotate \
        --output-type z \
        --annotations ${dbsnp_vcf} \
        --columns INFO/COMMON,INFO/PM,INFO/NSF,INFO/NSM,INFO/NSN,INFO/SAO,INFO/KGPROD,INFO/KGValidated,INFO/MUT,INFO/WTD,INFO/VLD,INFO/RS,INFO/PMC \
        --output ${sample_name}/${misc}/${sample_name}_dbsnp_annotated.vcf.gz \
        ${to_update_vcf}
    }
    output {
        File updated_vcf = "${sample_name}/${misc}/${sample_name}_dbsnp_annotated.vcf.gz"
    }
}

task edit_rna_editing {
    String misc
    String sample_name
    File darned_db
    File radar_db
    File to_update_vcf
    command {
        filter_snps_rna_editing.py \
        --darned ${darned_db} \
        --radar ${radar_db} \
        ${to_update_vcf} \
        ${sample_name}/${misc}/${sample_name}_editing_filtered.vcf
    }
    output {
        File updated_vcf = "${sample_name}/${misc}/${sample_name}_editing_filtered.vcf"
    }
}

task filter_vcf_cancer {
    String misc
    String sample_name
    File in_vcf
    command {
        filter_vcf_for_cancer.py ${in_vcf} ${sample_name}/${misc}/${sample_name}_cancer_filter.vcf
    }
    output {
        File out_vcf = "${sample_name}/${misc}/${sample_name}_cancer_filter.vcf"
    }
}

task add_rg {
    String misc
    String path
    String sample_name
    File in_bam
    command {
        java -jar ${path}/AddOrReplaceReadGroups.jar \
        I=${in_bam} \
        O=${sample_name}/${misc}/${sample_name}_rg.bam \
        SO=coordinate \
        RGID=id \
        RGLB=library \
        RGPL=ILLUMINA \
        RGPU=machine RGSM=${sample_name}
    }
    output {
        File out_bam = "${sample_name}/${misc}/${sample_name}_rg.bam"
    }
}

task gatk_base_recalibrator {
    String misc
    String path
    String sample_name
    File realigned_bam
    File reference
    File known_sites_db
    command {
        java -Xmx4g -jar ${path}/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -I ${realigned_bam} \
        -R ${reference} \
        --out ${sample_name}/${misc}/${sample_name}_recalibrated.tsv \
        -knownSites ${known_sites_db}
    }
    output {
        File out_table = "${sample_name}/${misc}/${sample_name}_recalibrated.tsv"
    }
}

task gatk_haplotype_caller {
    String sample_name
    String path
    File reference
    File bam
    command {
        java -jar ${path}/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R ${reference} \
        -I ${bam} \
        -recoverDanglingHeads \
        -dontUseSoftClippedBases \
        -stand_call_conf 20.0 \
        -stand_emit_conf 20.0 \
        --out ${sample_name}_original.vcf
    }
    output {
        File called_vcf = "${sample_name}_original.vcf"
    }
}

task gatk_indel_realigner {
    String misc
    String path
    String sample_name
    File reference
    File in_bam
    File intervals
    File known_sites_db
    command {
        java -Xmx4g -jar ${path}/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R ${reference} \
        -I ${in_bam} \
        -targetIntervals ${intervals} \
        --out ${sample_name}/${misc}/${sample_name}_realigned.bam \
        -known ${known_sites_db}
    }
    output {
        File out_bam = "${sample_name}/${misc}/${sample_name}_realigned.bam"
    }
}

task gatk_initial_filtering {
    String path
    String sample_name
    File reference
    File to_update_vcf
    command {
        java -jar ${path}/GenomeAnalysisTK.jar \
        -T VariantFiltration \
        -R ${reference} \
        -V ${to_update_vcf} \
        -window 35 \
        -cluster 3 \
        -filterName FS \
        -filter "FS > 30.0" \
        -filterName QD \
        -filter "QD < 2.0" \
        --out ${sample_name}_initial_filter.vcf
    }
    output {
        File updated_vcf = "${sample_name}_initial_filter.vcf"
    }
}

task gatk_realigner_target_creator {
    String misc
    String path
    String sample_name
    File reference
    File in_bam
    File known_sites_db
    command {
        java -Xmx2g -jar ${path}/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R ${reference} \
        -I ${in_bam} \
        --out ${sample_name}/${misc}/${sample_name}_realigner.intervals \
        --known ${known_sites_db}
    }
    output {
        File intervals = "${sample_name}/${misc}/${sample_name}_realigner.intervals"
    }
}

task gatk_recalibrated_bam {
    String path
    String sample_name
    File reference
    File realigned_bam
    File recalibration_table
    command {
        java -Xmx2g -jar ${path}/GenomeAnalysisTK.jar \
        -R ${reference} \
        -T PrintReads \
        --out ${sample_name}_recalibrated.bam \
        -I ${realigned_bam} \
        --BQSR ${recalibration_table}
    }
    output {
        File out_bam = "${sample_name}_recalibrated.bam"
    }
}

task gatk_split_cigar_reads {
    String misc
    String path
    String sample_name
    File reference
    File in_bam
    command {
        java -jar ${path}/GenomeAnalysisTK.jar \
        -T SplitNCigarReads \
        -R ${reference} \
        -I ${in_bam} \
        -o ${sample_name}/${misc}/${sample_name}_split.bam \
        -rf ReassignOneMappingQuality \
        -RMQF 255 \
        -RMQT 60 \
        -U ALLOW_N_CIGAR_READS
    }
    output {
        File out_bam = "${sample_name}/${misc}/${sample_name}_split.bam"
    }
}

task gatk_variants_table {
    String path
    String sample_name
    File reference
    File in_vcf
    command {
        java -jar ${path}/GenomeAnalysisTK.jar \
        -R ${reference} \
        -T VariantsToTable \
        -V ${in_vcf} \
        -F CHROM -F POS -F REF -F ALT -F GENE -F DP -F QUAL -F MQ -F SAO \
        -F NSF -F NSM -F NSN -F TUMOR -F TISSUE -F COSMIC_ID -F KGPROD -F RS \
        -F PMC -F CRAVAT_PVALUE -F CRAVAT_FDR -F VEST_PVALUE -F VEST_FDR \
        --allowMissingData \
        --unsafe LENIENT_VCF_PROCESSING \
        -o ${sample_name}_cancer.tab
    }
    output {
        File out_tab = "${sample_name}_cancer.tab"
    }
}

task groom_cravat_coding_results {
    String misc
    String sample_name
    File in_tab
    command {
        groom_cravat_annotation.py ${in_tab} ${sample_name}/${misc}/${sample_name}_variant_result_updated.tsv
    }
    output {
        File out_tab = "${sample_name}/${misc}/${sample_name}_variant_result_updated.tsv"
    }
}

task groom_cravat_noncoding_results {
    String misc
    String sample_name
    File in_tab
    command {
        groom_cravat_annotation.py ${in_tab} ${sample_name}/${misc}/${sample_name}_variant_non_coding_result_updated.tsv
    }
    output {
        File out_tab = "${sample_name}/${misc}/${sample_name}_variant_non_coding_result_updated.tsv"
    }
}

task groom_initial_vcf {
    String sample_name
    File to_update_vcf
    command {
        groom_vcf.py ${to_update_vcf} ${sample_name}_initial_calls_groomed.vcf
    }
    output {
        File updated_vcf = "${sample_name}_initial_calls_groomed.vcf"
    }
}

task groom_cravate_vcf {
    String sample_name
    File to_update_vcf
    command {
        groom_vcf.py ${to_update_vcf} ${sample_name}_cancer.vcf
    }
    output {
        File updated_vcf = "${sample_name}_cancer.vcf"
    }
}

task index_aligned_bam {
    String sample_name
    File in_bam
    String index_misc_dir
    command {
        samtools index ${in_bam}
    }
    output {
        File in_bai = "${sample_name}/${index_misc_dir}/${sample_name}_aligned.sortedByCoord.out.bam.bai"
    }
}

task index_tab {
    File in_tab
    command {
        tabix -f -s 1 -b 2 -e 2 -S 12 ${in_tab}
    }
    output {
        File out_tabix = "${in_tab}.tbi"
    }
}

task index_vcf {
    File in_vcf
    command {
        bcftools index ${in_vcf}
    }
    output {
        File out_vcf = "${in_vcf}.csi"
    }
}

task mark_duplicates {
    String misc
    String path
    String sample_name
    File in_bam
    command {
        java -jar ${path}/MarkDuplicates.jar \
        I=${in_bam} \
        O=${sample_name}/${misc}/${sample_name}_dedubbed.bam \
        CREATE_INDEX=true \
        M=${sample_name}/${misc}/${sample_name}_qc_metrics.txt
    }
    output {
        File out_bam = "${sample_name}/${misc}/${sample_name}_dedubbed.bam"
        File qc_metrics = "${sample_name}/${misc}/${sample_name}_qc_metrics.txt"
    }
}

task snpeff {
    String misc
    String path
    String sample_name
    File to_update_vcf
    command {
        bgzip -cd ${to_update_vcf} \
            | java -jar ${path}/snpEff.jar \
                -nostats \
                -noLof \
                -no-downstream \
                -no-upstream hg19 > ${sample_name}/${misc}/${sample_name}_snpeff.vcf
    }
    output {
        File updated_vcf = "${sample_name}/${misc}/${sample_name}_snpeff.vcf"
    }
}

task subset_to_snps {
    String sample_name
    String misc_dir
    File to_update_vcf
    command {
        reduce_vcf_to_snps.py ${to_update_vcf} ${sample_name}/${misc_dir}/${sample_name}_initial_snp.vcf
    }
    output {
        File updated_vcf = "${sample_name}/${misc_dir}/${sample_name}_initial_snp.vcf"
    }
}

task update_snpeff_annotations {
    String misc
    String sample_name
    File to_update_vcf
    command {
        update_snpeff_annotations.py ${to_update_vcf} ${sample_name}/${misc}/${sample_name}_snp_eff_ann.vcf
    }
    output {
        File updated_vcf = "${sample_name}/${misc}/${sample_name}_snp_eff_ann.vcf"
    }
}
