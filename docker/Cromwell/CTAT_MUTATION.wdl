workflow ctat_rnaseq_mutation {
    File in_reference
    File in_alignment_index
    File in_bed
    File in_cosmic_vcf
    File in_cravate_header
    File in_dbsnp_vcf
    File in_left
    File in_right
    String in_misc = "misc"
    String name

    call align {
        input: sample_name=name, align_index=in_alignment_index, align_left_sample=in_left, align_right_sample=in_right, align_misc_dir=in_misc
    }
    call index_aligned_bam {
        input: sample_name=name, in_bam=align.aligned_bam, index_misc_dir=in_misc
    }
    call gatk_add_rg {
        input: sample_name=name, in_bam=align.aligned_bam
    }
    call mark_duplicates {
        input: in_bam=gatk_add_rg.out_bam
    }
    call gatk_split_cigar_reads {
        input: sample_name=name, reference=in_reference, in_bam=mark_duplicates.out_bam
    }
    call gatk_realigner_target_creator {
        input: sample_name=name, reference=in_reference, known_sites_db=in_dbsnp_vcf, in_bam=gatk_split_cigar_reads.out_bam
    }
    call gatk_indel_realigner {
        input: sample_name=name, reference=in_reference, known_sites_db=in_dbsnp_vcf, intervals=gatk_realigner_target_creator.intervals, in_bam=gatk_split_cigar_reads.out_bam
    }
    call gatk_base_recalibrator {
        input: sample_name=name, reference=in_reference, known_sites_db=in_dbsnp_vcf, realigned_bam=gatk_indel_realigner.out_bam
    }
    call gatk_recalibrated_bam {
        input: sample_name=name, reference=in_reference, realigned_bam=gatk_indel_realigner.out_bam, recalibration_table=gatk_base_recalibrator.out_table
    }
    call gatk_haplotype_caller {
        input: sample_name=name, reference=in_reference, bam=gatk_recalibrated_bam.out_bam
    }
    call gatk_initial_filtering {
        input: sample_name=name, reference=in_reference, to_update_vcf=gatk_haplotype_caller.called_vcf
    }
    call groom_initial_vcf {
        input: to_update_vcf=gatk_initial_filtering.updated_vcf
    }
    call subset_to_snps {
        input: to_update_vcf=groom_initial_vcf.updated_vcf
    }
    call edit_rna_editing {
        input: darned_db=in_darned_db, radar_db=in_radar_db, to_update_vcf=subset_to_snps.updated_vcf
    }
    call compress_vcf as compress_editing_vcf {
        input: to_compress_vcf=edit_rna_editing.updated_vcf
    }
    call index_vcf as index_editing_vcf {
        input: in_vcf=compress_editing_vcf.compressed_vcf
    }
    call dbsnp_annotate {
        input: dbsnp_vcf=in_dbsnp_vcf, to_update_vcf=index_editing_vcf.out_vcf
    }
    call snpeff {
        input: to_update_vcf=dbsnp_annotate.updated_vcf
    }
    call update_snpeff_annotations {
        input: to_update_vcf=snpeff.updated_vcf
    }
    call compress_vcf as compress_snpeff {
        input: to_compress_vcf=update_snpeff_annotations.updated_vcf
    }
    call index_vcf as index_snpeff {
        input: in_vcf=compress_snpeff.out_vcf
    }
    call cosmic_annotate {
        input: cosmic_vcf=in_cosmic_vcf, to_annotate_vcf=index_snpeff.out_vcf
    }
    call filter_vcf_cancer {
        input: in_vcf=cosmic_annotate.cosmic_annotated
    }
    call cravat_annotate_get {
        input: sample_name=name, variants_vcf=filter_vcf_cancer.out_vcf
    }
    call cravat_unzip_response {
        input: sample_name=name, to_zip=cravat_annotate_get.out_annotations
    }
    call cravate_move_response {
        input: cravate_data=cravat_unzip_response.unzipped
    }
    call groom_cravat_coding_results {
        input: in_tab=
    }
    call groom_cravat_noncoding_results {
        input: in_tab=
    }
    call compress_tab {
        input: to_compress_tab=
    }
    call index_tab {
        input: in_tab
    }
    call compress_tab {
        input: to_compress_tab=
    }
    call index_tab {
        input: in_tab
    }
    call cravat_annotate_add {
        input: cravate_header=in_cravate_header,cravate_annotations=,in_vcf=
    }
    call cravat_annotate_add {
        input: cravate_header=in_cravate_header,cravate_annotations=,in_vcf=
    }
    call cravate_filter {
        input: in_vcf=
    }
    call groom_cravate_vcf {
        input: to_update_vcf=
    }
    call gatk_variants_table {
        input: reference=in_reference,in_vcf=
    }
    call make_mutation_json {
        input: cancer_bed=in_bed,cancer_dir=,cancer_idx=,cancer_tab=,cancer_bam=,cancer_bai=
    }
    call copy_bed {
        input: move_bed=in_bed
    }
}

task align {
    String sample_name
    File align_index
    File align_left_sample
    File align_right_sample
    String align_misc_dir
    command {
        STAR --genomeDir ${align_index} \
             --runThreadN 8 \
             --readFilesIn ${align_left_sample} \
                           ${align_right_sample} \
             --outSAMtype BAM SortedByCoordinate \
             --twopassMode Basic \
             --limitBAMsortRAM 30000000000 \
             --outFileNamePrefix ${align_misc_dir}
    }
    output {
        File aligned_bam = ${sample_name}/${align_misc_dir}/Aligned.out.bam
    }
}

task compress_tab {
    File to_compress_tab
    command {
        bgzip -c ${to_compress_tab} ${compressed_tab}
    }
    output {
        File compressed_tab
    }
}

task compress_vcf {
    File to_compress_vcf
    command {
        bgzip -c ${to_compress_vcf} > ${compressed_vcf}
    }
    output {
        File compressed_vcf = ${to_compress_vcf}.gz
    }
}

task copy_bed {
    String sample_name
    File move_bed
    command {
        cp ${move_bed} ${location}
    }
    output {
        File location =
    }
}

task cosmic_annotate {
    String sample_name
    File to_annotate_vcf
    File cosmic_vcf
    command {
        bcftools annotate --output-type z \
            --annotations ${cosmic_vcf} \
            --columns INFO/COSMIC_ID,INFO/TISSUE,INFO/TUMOR,INFO/FATHMM,INFO/SOMATIC \
            --output${cosmic_annotated} ${to_annotate_vcf}
    }
    output {
        File cosmic_annotated = ${sample_name}_cosmic_ann.vcf
    }
}

task cravat_annotate_add {
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
            --output ${out_vcf} \
            ${in_vcf}
    }
    output {
        File out_vcf =
    }
}

task cravat_annotate_get {
    String sample_name
    File variants_vcf
    command {
        annotate_with_cravat.py --classifier Blood-Lymphocyte \
            --email ttickle@broadinstitute.org \
            --max_attempts 180 \
            --wait 60 \
            ${variants_vcf} ${out_annotations}
    }
    output {
        File out_annotations = ${sample_name}_cravate_data.gz
    }
}

task cravate_filter {
    String sample_name
    File in_vcf
    command {
        bcftools filter \
            --include "CHASM_PVALUE < 0.3 || VEST_PVALUE < 0.3" \
            --output-type v \
            --output ${out_vcf} ${in_vcf}
    }
    output {
        File out_vcf =
    }
}

task cravate_move_response {
    String cravate_data
    command {
        cp {${cravate_data}/*/Variant.Result.tsv,${cravate_data}/*/Variant_Non-coding.Result.tsv} ${cravate_data}
    }
    output {
        File location_noncoding = ${cravate_data}/Variant.Result.tsv
        File location_noncoding = ${cravate_data}/Variant_Non-coding.Result.tsv
    }
}

task cravat_unzip_response {
    String sample_name
    File to_zip
    command {
        unzip -d ${to_zip} ${zipped}
    }
    output {
        File unzipped = ${sample_name}_extracted_data
    }
}

task dbsnp_annotate {
    String sample_name
    File dbsnp_vcf
    File to_update_vcf
    command {
        bcftools annotate \
            --output-type z \
            --annotations ${dbsnp_vcf} \
            --columns INFO/COMMON,INFO/PM,INFO/NSF,INFO/NSM,INFO/NSN,INFO/SAO,INFO/KGPROD,INFO/KGValidated,INFO/MUT,INFO/WTD,INFO/VLD,INFO/RS,INFO/PMC \
            --output ${updated_vcf} \
            ${to_update_vcf}
    }
    output {
        File updated_vcf = 
    }
}

task edit_rna_editing {
    String sample_name
    File darned_db
    File radar_db
    File to_update_vcf
    command {
        filter_snps_rna_editing.py \
            --darned ${darned_db} \
            --radar ${radar_db} \
            ${to_update_vcf} \
            ${updated_vcf}
    }
    output {
        File updated_vcf = ${sample_name}_editing_filtered.vcf
    }
}

task filter_vcf_cancer {
    String sample_name
    File in_vcf
    command {
        filter_vcf_for_cancer.py ${in_vcf} ${out_vcf}
    }
    output {
        File out_vcf = ${sample_name}_cancer_filter.vcf
    }
}

task gatk_add_rg {
    String sample_name
    File in_bam
    command {
        java -jar AddOrReplaceReadGroups.jar \
            I=${in_bam} \
            O=${out_bam} \
            SO=coordinate \
            RGID=id \
            RGLB=library \
            RGPL=ILLUMINA \
            RGPU=machine RGSM=FLI1_left
    }
    output {
        File out_bam = ${sample_name}_rg.bam
    }
}

task gatk_base_recalibrator {
    String sample_name
    File realigned_bam
    File reference
    File known_sites_db
    command {
        java -Xmx4g -jar GenomeAnalysisTK.jar \
            -T BaseRecalibrator \
            -I ${realigned_bam} \
            -R ${reference} \
            --out ${out_table} \
            -knownSites ${known_sites_db}
    }
    output {
        File out_table = ${sample_name}_recalibrated.tsv
    }
}

task gatk_haplotype_caller {
    String sample_name
    File reference
    File bam
    command {
        java -jar GenomeAnalysisTK.jar \
            -T HaplotypeCaller \
            -R ${reference} \
            -I ${bam} \
            -recoverDanglingHeads \
            -dontUseSoftClippedBases \
            -stand_call_conf 20.0 \
            -stand_emit_conf 20.0 \
            --out ${called_vcf}
    }
    output {
        File called_vcf = ${sample_name}_original.vcf
    }
}

task gatk_indel_realigner {
    String sample_name
    File reference
    File in_bam
    File intervals
    File known_sites_db
    command {
        java -Xmx4g -jar GenomeAnalysisTK.jar \
            -T IndelRealigner \
            -R ${reference} \
            -I ${in_bam} \
            -targetIntervals ${intervals} \
            --out ${out_bam} \
            -known ${known_sites_db}
    }
    output {
        File out_bam = ${sample_name}_realigned.bam
    }
}

task gatk_initial_filtering {
    String sample_name
    File reference
    File to_update_vcf
    command {
        java -jar GenomeAnalysisTK.jar \
            -T VariantFiltration \
            -R ${reference} \
            -V ${to_update_vcf} \
            -window 35 \
            -cluster 3 \
            -filterName FS \
            -filter "FS > 30.0" \
            -filterName QD \
            -filter "QD < 2.0" \
            --out ${updated_vcf}
    }
    output {
        File updated_vcf = ${sample_name}_initial_filter.vcf
    }
}

task gatk_realigner_target_creator {
    String sample_name
    File reference
    File in_bam
    File known_sites_db
    command {
        java -Xmx2g -jar GenomeAnalysisTK.jar \
            -T RealignerTargetCreator \
            -R ${reference} \
            -I ${in_bam} \
            --out ${intervales} \
            --known ${known_sites_db}
    }
    output {
        File intervals = ${sample_name}_realigner.intervals
    }
}

task gatk_recalibrated_bam {
    String sample_name
    File reference
    File realigned_bam
    File recalibration_table
    command {
        java -Xmx2g -jar GenomeAnalysisTK.jar \
            -R ${reference} \
            -T PrintReads \
            --out ${out_bam} \
            -I ${realigned_bam} \
            --BQSR ${recalibration_table}
    }
    output {
        File out_bam = ${sample_name}_recalibrated.bam
    }
}

task gatk_split_cigar_reads {
    String sample_name
    File reference
    File in_bam
    command {
        java -jar GenomeAnalysisTK.jar \
            -T SplitNCigarReads \
            -R ${reference} \
            -I ${in_bam} \
            -o ${out_bam} \
            -rf ReassignOneMappingQuality \
            -RMQF 255 \
            -RMQT 60 \
            -U ALLOW_N_CIGAR_READS
    }
    output {
        File out_bam = ${sample_name}_split.bam
    }
}

task gatk_variants_table {
    String sample_name
    File reference
    File in_vcf
    command {
        java -jar GenomeAnalysisTK.jar \
            -R ${reference} \
            -T VariantsToTable \
            -V ${in_vcf} \
            -F CHROM -F POS -F REF -F ALT -F GENE -F DP -F QUAL -F MQ -F SAO \
            -F NSF -F NSM -F NSN -F TUMOR -F TISSUE -F COSMIC_ID -F KGPROD -F RS \
            -F PMC -F CRAVAT_PVALUE -F CRAVAT_FDR -F VEST_PVALUE -F VEST_FDR \
            --allowMissingData \
            --unsafe LENIENT_VCF_PROCESSING \
            -o ${out_table}
    }
    output {
        File out_tab =
    }
}

task groom_cravat_coding_results {
    String sample_name
    File in_tab
    command {
        groom_cravat_annotation.py ${in_tab} ${out_tab}
    }
    output {
        File out_tab =
    }
}

task groom_cravat_noncoding_results {
    String sample_name
    File in_tab
    comand {
        groom_cravat_annotation.py ${in_tab} ${out_tab}
    }
    output {
        File out_tab =
    }
}

task groom_initial_vcf {
    String sample_name
    File to_update_vcf
    command {
        groom_vcf.py ${to_update_vcf} ${updated_vcf}
    }
    output {
        File updated_vcf = ${sample_name}_groomed.vcf
    }
}

task groom_cravate_vcf {
    String sample_name
    File to_update_vcf
    command {
        groom_vcf.py ${to_update_vcf} ${updated_vcf}
    }
    output {
        File updated_vcf = ${sample_name}_groomed_cravate.vcf
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
        File in_bai = ${sample_name}/${index_mimsc_dir}/Aligned.out.bai
    }
}

task index_tab {
    File in_tab
    command {
        tabix -f -s 1 -b 2 -e 2 -S 12 ${in_tab}
    }
    output {
        File out_tabix
    }
}

task index_vcf {
    File in_vcf
    command {
        bcftools index ${in_vcf}
    }
    output {
        File out_vcf
    }
}

task make_mutation_json {
    String sample_name
    File cancer_bam
    File cancer_bai
    File cancer_bed
    File cancer_dir
    File cancer_idx
    File cancer_tab
    command {
        make_mutation_inspector_json.py \
            --sample ${cancer_dir} \
            --tab ${cancer_tab} \
            --bam ${cancer_bam} \
            --bam_index ${cancer_bai} \
            --bed ${cancer_bed} \
            --bed_index ${cancer_idx} \
            ${out_json}
    }
    output {
        File out_json =
    }
}

task mark_duplicates {
    String sample_name
    File in_bam
    command {
        java -jar MarkDuplicates.jar \
            I=${in_bam} \
            O=${out_bam} \
            CREATE_INDEX=true \
            M=${qc_metrics}
    }
    output {
        File out_bam = ${sample_name}_dedubbed.bam 
        File qc_metrics = ${sample_name}_qc_metrics.txt
    }
}

task snpeff {
    String sample_name
    File to_update_vcf
    command {
        bgzip -cd ${to_update_vcf} \
            | java -jar snpEff.jar \
                -nostats \
                -noLof \
                -no-downstream \
                -no-upstream hg19 > ${updated_vcf}
    }
    output {
        File updated_vcf =
    }
}

task subset_to_snps {
    String sample_name
    File to_update_vcf
    command {
        reduce_vcf_to_snps.py ${to_update_vcf} ${updated_vcf}
    }
    output {
        File updated_vcf =
    }
}

task update_snpeff_annotations {
    String sample_name
    File to_update_vcf
    command {
        update_snpeff_annotations.py ${to_update_vcf} ${updated_vcf}
    }
    output {
        File updated_vcf = ${sample_name}_snp_eff_ann.vcf
    }
}
