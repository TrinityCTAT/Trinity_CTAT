#!/usr/bin/env perl

use strict;
use warnings;
use CGI::Carp qw(fatalsToBrowser);
use CGI;
use File::Basename;


my $scriptname = basename($0);
my $indexed_flag = "false";


$|++;

main: {

    my $cgi = new CGI();

    print $cgi->header();

    print $cgi->start_html(-title => "Fusion Transcript Portal",
        -head => q{
<head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="description" content="">
    <meta name="author" content="">
    <link rel="shortcut icon" href="//igv.org/web/img/favicon.ico">
    <title>Integrative Genomics Viewer</title>

    <!-- Bootstrap CSS - NOT REQUIRED FOR IGV -->
    <link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css">

    <!-- jQuery UI CSS -->
    <link rel="stylesheet" type="text/css"
          href="//ajax.googleapis.com/ajax/libs/jqueryui/1.11.2/themes/smoothness/jquery-ui.css"/>

    <!-- Font Awesome CSS -->
    <link rel="stylesheet" type="text/css" href="//maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css">

    <!-- IGV CSS -->
    <link rel="stylesheet" type="text/css" href="//igv.org/web/beta/igv.css">

    <!-- jQuery JS -->
    <script type="text/javascript" src="//ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
    <script type="text/javascript" src="//ajax.googleapis.com/ajax/libs/jqueryui/1.11.2/jquery-ui.min.js"></script>

    <!-- Bootstrap JS -  NOT REQUIRED FOR IGV -->
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/js/bootstrap.min.js"></script>

    <!-- IGV JS-->
    <script type="text/javascript" src="//igv.org/web/beta/igv-all.min.js"></script>


    
    <style type="text/css">

        #trackList {

            font-family: 'PT Sans', sans-serif;
            font-size: small;
            font-weight: 400;
        }

        div#trackList > div, div#trackList > h3 {

            color: #444;
            margin-left: 48px;
            margin-top: 4px;
            margin-bottom: 4px;

            padding-left: 32px;

            width: 300px;
        }

        div#trackList > div:hover,
        div#trackList > div:focus,
        div#trackList > div:active {
            cursor: pointer;
            color: white;
            background-color: rgba(49, 106, 246, 1);
        }


       .igv-viewport-div {
           overflow-y: auto;
       }

   </style>



</head>    




});
    
    my $params_href = $cgi->Vars();
    
    my $dataset = $params_href->{dataset} or die "Error, must specify dataset param";
    
    &js_IGV("CTAT_GALAXY_DATA/$dataset");
    
    print $cgi->end_html();


    exit(0);
}


####
sub js_IGV {
    my ($project_dir) = @_;
    
    
    print "<h2>Fusion Inspector Web<h2>\n";
    
    
    &add_fusion_list("$project_dir/FInspector.fa.fai", "$project_dir/FInspector.fusion_predictions.txt");
    
    
    print <<__EOIGV;

    <div id="fusion_descr"></div>

    <div id="myDiv" style="padding-top: 50px;padding-bottom: 20px;"></div>


    <script type="text/javascript">

    \$(document).ready\(function () {

        var div = \$("#myDiv")\[0],

                      options = {
                        showKaryo: false,
                        showNavigation: true,
                        fastaURL: "$project_dir/FInspector.fa",
                        cytobandURL: "$project_dir/cytoBand.txt",
                        tracks: [
                            {
                              type: "sequence",
                              order: 1
                            },

                            {
                              label: "ref_annot",
                              url: "$project_dir/FInspector.bed.sorted.bed.gz",
                              displayMode: "SQUISHED",
                              indexed: $indexed_flag,
                              order: 10
                                  
                            },

                            

                            {
                              label: "JuncSpanView",
                              type: "FusionJuncSpan",
                              url: "$project_dir/FInspector.igv.FusionJuncSpan",
                              displayMode: "SQUISHED",
                              indexed: false,
                              order: 20
                            },
                 



                            {
                              label: "Trinity_fusion",
                              url: "$project_dir/FInspector.gmap_trinity_GG.fusions.gff3.bed.sorted.bed.gz",
                              displayMode: "EXPANDED",
                              color: "green",
                              indexed: $indexed_flag,
                              order: 30

                            },
            

                            {
                              url: "$project_dir/FInspector.junction_reads.bam.bed.sorted.bed.gz",
                              label: "FInspector_junction_reads",
                              displayMode: "SQUISHED",
                              minHeight: 30,
                              maxHeight: 100,
                              color: "cyan",
                              indexed: $indexed_flag,
                              order: 40
                            
                            },

                            {
                              url: "$project_dir/FInspector.junction_reads.bam",
                              label: "FInspector_junction_reads_bam",
                              height: 100,
                              indexed: true,
                              visibilityWindow: 2000000,
                              order: 45
                            
                            },


                            
                            {
                              url: "$project_dir/FInspector.spanning_reads.bam.bed.sorted.bed.gz",
                              label: "FInspector_spanning_reads",
                              displayMode: "SQUISHED",
                              minHeight: 30,
                              maxHeight: 100,
                              color: "purple",
                              indexed: $indexed_flag,
                              order: 50
                            },

{
                              url: "$project_dir/FInspector.spanning_reads.bam",
                              label: "FInspector_spanning_reads_bam",
                              height: 100,
                              indexed: true,
                              visibilityWindow: 2000000,
                              order: 55
                            
                            }


                            
        ]
                      };

                      igv.createBrowser(div, options);


    


    });

</script>


__EOIGV

    ;

    return;


}




sub add_fusion_list {
    my ($file, $prog_info) = @_;

    my %fusion_to_prog_list;
    {
        
        my %fusion_to_support;

        open (my $fh, $prog_info) or die "Error, cannot open file $prog_info";
        while (<$fh>) {
            if (/^\#/) { next; }
            chomp;
            my @x = split(/\t/);
            my ($fusion_name, $score, $geneA, $chr_brkpt_A, $geneB, $chr_brkpt_B, $splice_type, $junc_support, $span_support, @rest) = @x;
            if ($junc_support) {
                push (@{$fusion_to_support{"$geneA--$geneB"}}, "$chr_brkpt_A-$chr_brkpt_B\[J-$junc_support,S-$span_support,$splice_type]");
            }
        }
        close $fh;
        
        foreach my $fusion (keys %fusion_to_support) {
            my @entries = sort (@{$fusion_to_support{$fusion}});
            $fusion_to_prog_list{$fusion} = join("; ", @entries);
        }
    }

    ## list the contigs.
    my @fusion_list;
    open (my $fh, "$file") or die $!; 
    while (<$fh>) {
        chomp;
        my ($fusion, @rest) = split(/\t/);
        push (@fusion_list, $fusion);
    }
    close $fh;
    
    @fusion_list = sort @fusion_list;
    
    my $html = "<select id='fusion_list'>\n";
    $html .= "    <option>-- fusion --</option>\n";
    foreach my $fusion (sort keys %fusion_to_prog_list) {
        #print "<div id=\'$fusion\' class='fusion_type' onclick=\"selected_fusion(\'$fusion\');\">$fusion</div>\n";
        my $fusion_support = $fusion_to_prog_list{$fusion};
        $html .= "   <option value=\'$fusion\'>$fusion $fusion_support</option>\n";
    }
    $html .= "</select>\n";

    $html .= <<__EOSCR;

    <script>
        document.getElementById("fusion_list").onchange = function() {
           var fusion = this.value;
           igv.browser.search(fusion);
          
        }
    </script>


__EOSCR

        ;
    
    print $html;

    return;
}


