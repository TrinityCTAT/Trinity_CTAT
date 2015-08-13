#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";



my @files_to_keep = qw (cleanme.pl 
                        reads.left.simPE.fq.gz
                        reads.right.simPE.fq.gz
                        runMe.pl
                        test_fusions.list
                        test_fusions.list2
                        test_fusions.list3

test_FusionInspector.pl
__runMe_kcopipe.sh
test_TrinityFusion.pl
test_STAR-fusion.pl
test_PRADA.pl
test_SoapFuse.pl
test_STAR-FGene.pl
                        );


my %keep = map { + $_ => 1 } @files_to_keep;


`rm -rf ./Gsnap_Fusion` if (-d "Gsnap_Fusion");
`rm -rf ./Star_Fusion` if (-d "Star_Fusion");
`rm -rf ./Trinity_Fusion` if (-d "Trinity_Fusion");
`rm -rf ./Fusion_Inspector` if (-d "Fusion_Inspector");
`rm -rf ./trinity_out_dir` if (-d "trinity_out_dir");
`rm -rf ./prada_outdir` if (-d "prada_outdir");
`rm -rf ./soapfuse_outdir` if (-d "soapfuse_outdir");
`rm -rf ./Star_FGene` if (-d "Star_FGene");
`rm -rf ./_STAR*`;


foreach my $file (<*>) {
	
	if (-f $file && ! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}


exit(0);
