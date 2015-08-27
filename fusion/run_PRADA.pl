#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Cwd;

my $usage = "\n\tusage: $0 reads.left.fq reads.right.fq output_directory\n\n";

my $left_fq = $ARGV[0] or die $usage;
my $right_fq = $ARGV[1] or die $usage;
my $output_dir = $ARGV[2] or die $usage;



my $PRADA_DIR = $ENV{PRADA_DIR}; #"/seq/regev_genome_portal/SOFTWARE/PRADA/pyPRADA_1.2/";
unless ($PRADA_DIR) {
    die "Error, must specify env var PRADA_DIR to point to PRADA installation directory";
}


main: {

    $left_fq = &ensure_full_path($left_fq);
    $right_fq = &ensure_full_path($right_fq);
    $output_dir = &ensure_full_path($output_dir);

    unless (-d $output_dir) {
        &process_cmd("mkdir -p $output_dir");
    }
    
    chdir $output_dir or die "Error, cannot cd to $output_dir";
    
    my $left_fq_use = "illumina.end1.fastq";
    unless (-e $left_fq_use) {
        &make_fastq_file($left_fq, $left_fq_use);
    }
    
    my $right_fq_use = "illumina.end2.fastq";
    unless (-e $right_fq_use) {
        &make_fastq_file($right_fq, $right_fq_use);
    }
    
    my $read_length = &get_read_length($left_fq_use);
    my $pct_80_readlength = int(0.8 * $read_length);
    
    unless (-e "illumina.bam") {
        &process_cmd("echo null > illumina.bam");
    }
    
    ## pre-process fastqs, align to genome 
    my $cmd = "$PRADA_DIR/prada-preprocess-bi2 -ref $PRADA_DIR/conf.txt -tag prada -step 2_e1_1 -intermediate no -docstr ladeda -outdir ./prada_preprocess -sample illumina  -pbs run_prada_preprocess ";
    unless (-s "illumina.withRG.GATKRecalibrated.flagged.bam") {
        &process_cmd($cmd);
        
    
        $cmd = "sh prada_preprocess/run_prada_preprocess.pbs";
        &process_cmd($cmd);
    }

    ## run prada-fusion
    $cmd = "$PRADA_DIR/prada-fusion -bam illumina.withRG.GATKRecalibrated.flagged.bam -ref $PRADA_DIR/conf.txt -tag prada -mm 2 -junL $pct_80_readlength -outdir prada_fusions";
    &process_cmd($cmd);

    
    unlink($left_fq_use, $right_fq_use);
    

  update_summary_file:
    ## reformat:
    my $prada_summary_file = "$output_dir/prada_fusions/prada.fus.summary.txt";
    unless (-s $prada_summary_file) {
        die "Error, cannot locate output file: $prada_summary_file";
    }
    my $adj_output_file = "$output_dir/prada.fusions.summary";
    open (my $fh, $prada_summary_file) or die "Error, cannot open file $prada_summary_file";
    open (my $ofh, ">$adj_output_file") or die "Error, cannot write to file $adj_output_file";
    my $header = <$fh>;
    print $ofh "#fusion_name\t$header";
    while (<$fh>) {
        my @x = split(/\t/);
        my ($geneA, $geneB) = ($x[0], $x[1]);
        unshift (@x, "$geneA--$geneB");
        print $ofh join("\t", @x);
    }
    close $fh;
    close $ofh;
    

    exit(0);
    
}


####
sub get_read_length {
    my ($fastq) = @_;

    my $num_reads = 10000;

    my $cmd = "sed -n '2~4p' $fastq | head -n$num_reads";
    my @results = `$cmd`;
    if ($?) {
        die "Error, cmd: $cmd died with ret $?";
    }
    chomp @results;

    my $sum_len = 0;
    foreach my $seq (@results) {
        $seq =~ s/\s//g;
        $sum_len += length($seq);
    }
    
    my $avg_len = int($sum_len/$num_reads);

    if ($avg_len < 30) {
        die "Error, avg read length found to be $avg_len .... reads seem to be too short for this.";
    }
    return($avg_len);
}


####
sub make_fastq_file {
    my ($in_file, $out_file) = @_;
    
    if ($in_file =~ /\.gz$/) {
        my $cmd = "gunzip -c $in_file > $out_file";
        &process_cmd($cmd);
    }
    else {
        # just symlink
        my $cmd = "ln -sf $in_file $out_file";
        &process_cmd($cmd);
    }
}


####
sub ensure_full_path {
    my ($path) = @_;

    if ($path !~ m|^/|) {
        $path = cwd() . "/$path";
    }

    return($path);
}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}


__END__


## BWA alignment

# align end1
/Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/bwa-0.5.7-mh/bwa aln -t 12 /RIS/home/verhaakgroup/PRADA/hg19broad/Ensembl64.transcriptome.plus.genome.fasta /Users/bhaas/BioIfx/PRADA/tmp/illumina.end1.fastq > illumina.end1.sai
/Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/bwa-0.5.7-mh/bwa samse -s -n 100 /RIS/home/verhaakgroup/PRADA/hg19broad/Ensembl64.transcriptome.plus.genome.fasta illumina.end1.sai /Users/bhaas/BioIfx/PRADA/tmp/illumina.end1.fastq > illumina.end1.sam
/Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/samtools-0.1.16/samtools view -bS -o illumina.end1.bam illumina.end1.sam
/Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/samtools-0.1.16/samtools sort -n -m 1000000000 illumina.end1.bam illumina.end1.sorted

# align end2
/Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/bwa-0.5.7-mh/bwa aln -t 12 /RIS/home/verhaakgroup/PRADA/hg19broad/Ensembl64.transcriptome.plus.genome.fasta /Users/bhaas/BioIfx/PRADA/tmp/illumina.end2.fastq > illumina.end2.sai
/Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/bwa-0.5.7-mh/bwa samse -s -n 100 /RIS/home/verhaakgroup/PRADA/hg19broad/Ensembl64.transcriptome.plus.genome.fasta illumina.end2.sai /Users/bhaas/BioIfx/PRADA/tmp/illumina.end2.fastq > illumina.end2.sam
/Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/samtools-0.1.16/samtools view -bS -o illumina.end2.bam illumina.end2.sam
/Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/samtools-0.1.16/samtools sort -n -m 1000000000 illumina.end2.bam illumina.end2.sorted

# remap alignments end1
java -Djava.io.tmpdir=tmp/ -cp /Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/GATK//RemapAlignments.jar -Xmx8g org.broadinstitute.cga.tools.gatk.rna.RemapAlignments M=/RIS/home/verhaakgroup/PRADA/hg19broad/Ensembl64.transcriptome.plus.genome.map IN=illumina.end1.sorted.bam OUT=illumina.end1.remapped.bam R=/RIS/home/verhaakgroup/PRADA/hg19broad/Homo_sapiens_assembly19.fasta REDUCE=TRUE
/Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/samtools-0.1.16/samtools sort -n -m 1000000000 illumina.end1.remapped.bam illumina.end1.remapped.sorted

# remap alignments end2
java -Djava.io.tmpdir=tmp/ -cp /Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/GATK//RemapAlignments.jar -Xmx8g org.broadinstitute.cga.tools.gatk.rna.RemapAlignments M=/RIS/home/verhaakgroup/PRADA/hg19broad/Ensembl64.transcriptome.plus.genome.map IN=illumina.end2.sorted.bam OUT=illumina.end2.remapped.bam R=/RIS/home/verhaakgroup/PRADA/hg19broad/Homo_sapiens_assembly19.fasta REDUCE=TRUE
/Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/samtools-0.1.16/samtools sort -n -m 1000000000 illumina.end2.remapped.bam illumina.end2.remapped.sorted

# make paired-end alignments
java -Djava.io.tmpdir=tmp/ -Xmx8g -jar /Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/GATK//PairMaker.jar IN1=illumina.end1.remapped.sorted.bam IN2=illumina.end2.remapped.sorted.bam OUTPUT=illumina.paired.bam TMP_DIR=tmp/
/Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/samtools-0.1.16/samtools sort -m 1000000000 illumina.paired.bam illumina.paired.sorted


# quality score recalibration        
java -Xmx8g -jar /Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/Picard//AddOrReplaceReadGroups.jar I=illumina.paired.sorted.bam O=illumina.withRG.paired.sorted.bam RGLB=prada RGPL=illumina RGPU=prada RGSM=prada
/Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/samtools-0.1.16/samtools index illumina.withRG.paired.sorted.bam
java -Xmx8g -jar /Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/GATK//GenomeAnalysisTK.jar -l INFO -R /RIS/home/verhaakgroup/PRADA/hg19broad/Homo_sapiens_assembly19.fasta --default_platform illumina --knownSites /RIS/home/verhaakgroup/PRADA/hg19broad/dbsnp_135.b37.vcf -I illumina.withRG.paired.sorted.bam --downsample_to_coverage 10000 -T CountCovariates -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -nt 12 -recalFile illumina.orig.csv
java -Xmx8g -jar /Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/GATK//GenomeAnalysisTK.jar -l INFO -R /RIS/home/verhaakgroup/PRADA/hg19broad/Homo_sapiens_assembly19.fasta --default_platform illumina -I illumina.withRG.paired.sorted.bam -T TableRecalibration --out illumina.withRG.GATKRecalibrated.bam -recalFile illumina.orig.csv

# mark duplicates
java -Xmx8g -jar /Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/Picard//MarkDuplicates.jar I=illumina.withRG.GATKRecalibrated.bam O=illumina.withRG.GATKRecalibrated.flagged.bam METRICS_FILE=illumina.Duplicates_metrics.txt VALIDATION_STRINGENCY=SILENT TMP_DIR=tmp/
/Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/samtools-0.1.16/samtools index illumina.withRG.GATKRecalibrated.flagged.bam


java -Xmx8g -jar /Users/bhaas/BioIfx/PRADA/pyPRADA_1.2/tools/RNA-SeQC_v1.1.7.jar -ttype 2 -t /RIS/home/verhaakgroup/PRADA/hg19broad/Homo_sapiens.GRCh37.64.gtf -r /RIS/home/verhaakgroup/PRADA/hg19broad/Homo_sapiens_assembly19.fasta -s 'illumina|illumina.withRG.GATKRecalibrated.flagged.bam|Disc' -o illumina/
