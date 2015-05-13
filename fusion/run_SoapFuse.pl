#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use Cwd;
use FindBin;


my $SOAPFUSE_DIR = $ENV{SOAPFUSE_DIR}; # "/seq/regev_genome_portal/SOFTWARE/SoapFuse/SOAPfuse-v1.26";
#my $UTIL_DIR = "$FindBin::Bin/util";

unless ($SOAPFUSE_DIR) {
    die "Error, must set SOAPFUSE_DIR env var to path where SOAPFUSE is installed";
}


my $usage = "usage: $0 reads.left.fq reads.right.fq output_directory/ [NO_CLEANUP_FLAG]\n\n";
my $left_fq = $ARGV[0] or die $usage;
my $right_fq = $ARGV[1] or die $usage;
my $output_directory = $ARGV[2] or die $usage;
my $NO_CLEANUP_FLAG = $ARGV[3] || 0;


my $MIN_SEQ_LEN = 25;

my $soap_fuse_config_file = "$SOAPFUSE_DIR/config/gencode_v19.config";


$left_fq = &ensure_full_path($left_fq, 1);
$right_fq = &ensure_full_path($right_fq, 1);
$soap_fuse_config_file = &ensure_full_path($soap_fuse_config_file, 1);
$output_directory = &ensure_full_path($output_directory, 0);



main: {

    my $checkpoint = "$output_directory/SoapFuse.ok";
    if (-e $checkpoint) {
        print STDERR "Already ran SoapFuse here: $checkpoint.\n";
        exit(0);
    }
    

    unless (-d $output_directory) {
        mkdir $output_directory or die "Error, cannot mkdir $output_directory";
    }

    chdir $output_directory or die "Error, cannot cd to $output_directory";
    
    ## Prep reads in SoapFuse-expected file structure

    my $soap_sample_dir = "sample/lib";
    if (! -d $soap_sample_dir) {
        &process_cmd("mkdir -p $soap_sample_dir");
    }
    chdir $soap_sample_dir or die "Error, cannot chdir $soap_sample_dir";
    

    
    &ensure_gzip_fq_file($left_fq, "run_1.fq.gz");
    &ensure_gzip_fq_file($right_fq, "run_2.fq.gz");    


    
    ## Run SoapFuse
    chdir $output_directory or die "Error, cannot cd to $output_directory";

    my $read_length = &get_read_length($left_fq);

    my $samples_file = "samples_list.txt";
    open (my $ofh, ">$samples_file") or die "Error, cannot write to file $samples_file";
    print $ofh join("\t", "sample", "lib", "run", $read_length);
    close $ofh;

    my $cmd = "$SOAPFUSE_DIR/SOAPfuse-RUN.pl -c $SOAPFUSE_DIR/config/gencode_v19.config -fd . -l samples_list.txt -o soap_fuse_outdir ";
    
    if (! -e $checkpoint) {
        
        &process_cmd($cmd);
        
        
        ## extract the fusion reads:
        #my $cmd = "$UTIL_DIR/soapFuse_junc_read_extractor.pl ./soap_fuse_outdir/final_fusion_genes/sample/sample.final.Span_reads > ./soap_fuse_outdir/final_fusion_genes/sample/sample.final.Span_reads.fa";
        #&process_cmd($cmd);
        
        #$cmd = "$UTIL_DIR/soapFuse_junc_read_extractor.pl ./soap_fuse_outdir/final_fusion_genes/sample/sample.final.Junc_reads > ./soap_fuse_outdir/final_fusion_genes/sample/sample.final.Junc_reads.fa";
        #&process_cmd($cmd);
        
        
        #############
        ## Clean-up
        #############
        
        unless ($NO_CLEANUP_FLAG) {
            
            &process_cmd("mv ./soap_fuse_outdir/junction_seq .") if (-d "./soap_fuse_outdir/junction_seq");
            &process_cmd("mv ./soap_fuse_outdir/final_fusion_genes .") if (-d "./soap_fuse_outdir/final_fusion_genes");
            
            &process_cmd("rm -rf ./soap_fuse_outdir") if (-d "./soap_fuse_outdir");
            
            
        }
    
        &process_cmd("touch $checkpoint"); 

    }
        
    
    exit(0);
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


####
sub ensure_full_path {
    my ($path, $require_exists) = @_;

    unless ($path =~ m|^/|) {
        $path = cwd() . "/$path";
    }

    if ($require_exists && ! -e $path) {
        confess ("Error, cannot locate file: $path");
    }

    return($path);
}

####
sub get_read_length {
    my ($fq_file) = @_;

    my $cmd = ($fq_file =~ /\.gz$/) ? "gunzip -c $fq_file | head -n4 | tail -n1" : "head -n4 $fq_file | tail -n1";

    my $seq_line = `$cmd`;
    chomp $seq_line;

    unless (length($seq_line) > $MIN_SEQ_LEN) {
        confess "Error, sequence length of $seq_line is < min=$MIN_SEQ_LEN";
    }

    return(length($seq_line));
}



####
sub ensure_gzip_fq_file {
    my ($start_fq, $target_fq) = @_;
    
    if ($start_fq =~ /\.gz$/) {
        
        ## just symlink it:
        &process_cmd("ln -sf $start_fq $target_fq");
    }
    elsif (-e "$start_fq.gz") {
        # use the gzip-version
        &process_cmd("ln -sf $start_fq.gz $target_fq");
    }
    else {
        ## create a gzipped version, since soapfuse seems to want it.
        &process_cmd("gzip -c $start_fq > $target_fq");
    }

    return;
}
