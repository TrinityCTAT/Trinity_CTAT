#!/usr/bin/env Rscript

library( ggplot2 )
library( optparse )

C_I_PRIMARY_POS = 1
C_I_PRIMARY_REF = 2
C_I_PRIMARY_GT = 3
C_I_PRIMARY_DEPTH = 4
C_I_SECONDARY_POS = 5
C_I_SECONDARY_REF = 6
C_I_SECONDARY_GT = 7
C_I_SECONDARY_DEPTH = 8

C_STR_DETAIL_FILE = "detail_validation.pdf"

# Argument parsing
pArgs = OptionParser( usage ="%prog files1.tab files2.tab" )
pArgs = add_option( pArgs, c("-k","--title_key"),type="character",action="store",dest="str_title_key",default="Primary vs Secondary",help="Key identifying the contrast being visualized (eg \"DNA vs RNA\")")
pArgs = add_option( pArgs, c("-o","--group_output_dir"),type="character",action="store",dest="str_output_dir",default=NA,help="Output directory (required).")
lsArgs = parse_args( pArgs, positional_arguments=TRUE )

# Handle arguments
v_str_files = lsArgs$args
str_title_key = lsArgs$options$str_title_key
str_output_dir = lsArgs$options$str_output_dir
if( is.na(str_output_dir) )
{
  stop("Please include the output directory. Use -o.")
}
dir.create( str_output_dir )

# Information about the files as a group
# Holds file information for the global images at the end.
print( "Number of files reading:" )
print( length( v_str_files ) )
ls_group_sensitivity = c()
ls_group_specificity = c()
ls_group_files = c()

# Process each file
for( str_file in v_str_files )
{
  print( "Processing file:" )
  print( str_file )

  # id \t chr_position \t Ref_bas \t call_genotype \t depth \t chr_position \t ref_bas \t call_genotype \t depth
  # 1-4 is for the reference sample
  # 5-8 is for the secondary assay sample
  df_tab = read.table( str_file )

  # Make sure that the file does not include calls that do not get called in both primary and secondary sources
  vi_empty_calls = intersect( which( is.na( df_tab[[ C_I_PRIMARY_GT ]]  ) ), which( is.na( df_tab[[ C_I_SECONDARY_GT ]] ) ) )
  if( length( vi_empty_calls ) > 0 )
  {
    df_tab = df_tab[-1*vi_empty_calls,]
  }

  # Read depth can vary between evidence.
  # To try to normalize this, feature depth are divided by the sum of all the depth in the evidence
  df_tab[[ C_I_PRIMARY_DEPTH ]] = df_tab[[ C_I_PRIMARY_DEPTH ]] / sum( df_tab[[ C_I_PRIMARY_DEPTH ]], na.rm = TRUE )
  df_tab[[ C_I_SECONDARY_DEPTH ]] = df_tab[[ C_I_SECONDARY_DEPTH ]] / sum( df_tab[[ C_I_SECONDARY_DEPTH ]], na.rm = TRUE )

  # Different measurements used later on
  ## Primary
  ### Locations with no depth
  vi_primary_no_depth = which( is.na( df_tab[[ C_I_PRIMARY_DEPTH ]] ))
  ### Locations with depth
  vi_primary_read_evidence = which( !is.na( df_tab[[ C_I_PRIMARY_DEPTH ]] ))
  ### Primary calls
  vi_primary_calls = which( !is.na(  df_tab[[ C_I_PRIMARY_GT ]] ) )
  ## Secondary
  ### Locations with no depth
  vi_secondary_no_depth = which( is.na( df_tab[[ C_I_SECONDARY_DEPTH ]] ))
  ### Locations with depth
  vi_secondary_read_evidence = which( !is.na( df_tab[[ C_I_SECONDARY_DEPTH ]] ))
  i_count_secondary_read_evidence = length( vi_secondary_read_evidence )
  ### Calls with or without depth
  vi_secondary_calls = which( !is.na( df_tab[[ C_I_SECONDARY_GT ]] ) )
  ### Calls with depth
  vi_secondary_calls_with_evidence = intersect( vi_secondary_read_evidence, vi_secondary_calls )
  i_count_secondary_calls_with_evidence = length( vi_secondary_calls_with_evidence )
  ## Both
  vi_agreeing_calls = intersect( vi_primary_calls, vi_secondary_calls )
  ### Correct secondary calls (defined as secondary agreeing with primary)
  vi_secondary_calls_agreeing_primary = intersect( vi_secondary_calls_with_evidence, vi_agreeing_calls )
  i_count_secondary_calls_agreeing_primary = length( vi_secondary_calls_agreeing_primary )
  
  # False / True negative and postive rates for primary and secondary calls
  ## False negative rates
  ### Secondary: Calls not called in the secondary that had read depth and were called in the primary
  vi_secondary_false_negatives = setdiff( intersect( vi_secondary_read_evidence, vi_primary_calls ), vi_secondary_calls )
#  vi_secondary_false_negatives_depth = df_tab[[ C_I_PRIMARY_DEPTH ]][ vi_secondary_false_negatives ]
#  vi_order_secondary_false_negatives = order( vi_secondary_false_negatives_depth )
  ## False positive rates
  ### Secondary: Calls in the secondary that are not in the primary calls
  vi_secondary_false_positives = setdiff( vi_secondary_calls, vi_primary_calls )
#  vi_order_secondary_false_positives = order( df_tab[[ C_I_SECONDARY_DEPTH ]][ vi_secondary_false_positives ] )
  ## True negative rates
  ### Secondary: Not recorded
  ## True postive rates
  ### Secondary: Calls in the secondary with read depth that are also in the primary
  vi_secondary_true_positives = intersect( vi_secondary_calls, vi_primary_calls )
#  vi_order_secondary_true_positives = order( df_tab[[ C_I_SECONDARY_DEPTH ]][ vi_secondary_true_positives ] )

  # How many calls occured with no read depth to support them
  print( "Just a check: How many primary calls had no read depth" )
  print( length( intersect( vi_primary_calls, vi_primary_no_depth ) ) )
  print( "Just a check: How many secondary calls had no read depth" )
  print( length( intersect( vi_secondary_calls, vi_secondary_no_depth ) ) )

  ## Primary are MAF or DNA
  ## Secondary are DNA or RNA
  # How many primary calls are not supported by the secondary assay
  # How many secondary calls are not supported by the primary assay
  vi_primary_exclusive = intersect( vi_primary_calls, vi_secondary_no_depth )
  vi_secondary_exclusive = setdiff( vi_secondary_calls, vi_primary_no_depth )
  print( "Percent of the primary calls with no evidence in the secondary calls" )
  print( length( vi_primary_exclusive )/length( vi_primary_calls ) * 100 )
  print( "Percent of the secondary calls with no evidence in the primary calls." )
  print( length( vi_secondary_exclusive )/length( vi_secondary_calls ) * 100 )
  print( "Remember this is a problematic view for maf files given they do not have thier own depth files or it is the same as the DNA files." )

  # How does read depth associate with secondary assays being called
  # What are the distributions of read depth for good and bad calls
  # v_sec_calls_bad = setdiff( vi_secondary_calls_with_evidence , vi_secondary_calls_agreeing_primary )
#  v_histogram = c( df_tab[[ C_I_SECONDARY_DEPTH ]][ vi_secondary_true_positives ],
#                   df_tab[[ C_I_SECONDARY_DEPTH ]][ vi_secondary_false_positives ],        
#                   df_tab[[ C_I_PRIMARY_DEPTH ]][ vi_secondary_false_negatives ] )
#  i_bin = floor( max( v_histogram )/200 )
#  df_histogram = data.frame( values = v_histogram,
#                           groups = c(rep("True_positives",length( vi_secondary_true_positives )),
#                                      rep("False_Positives", length( vi_secondary_false_positives )),
#                                      rep("False_Negatives",length( vi_secondary_false_negatives ))))
#
#  print("Creating detail ROC")
#  print( file.path( str_output_dir, paste( basename( str_file ), "rates_hist", C_STR_DETAIL_FILE, sep = "_"  )) )
#  ggplot( df_histogram, aes(x=values,fill=groups) ) +
#    geom_histogram( data=subset(df_histogram,groups=="True_positives"), fill = "green", alpha = .4, binwidth=i_bin ) +
#    geom_histogram( data=subset(df_histogram,groups=="False_Positives"), fill = "orange", alpha = .4, binwidth=i_bin ) +
#    geom_histogram( data=subset(df_histogram,groups=="False_Negatives"), fill = "purple", alpha = .4, binwidth=i_bin ) +
#    geom_vline(aes(xintercept=20,colour="red"), linetype="dashed", size=.25) +
#    ggtitle( paste( "Depth distribution in secondary calls (",str_title_key,")") )
#  ggsave(file=file.path( str_output_dir, paste( basename( str_file ), "rates_hist", C_STR_DETAIL_FILE, sep = "_"  )))
  
  # Sensitivity (y-axis): Exome Calls with min rna-seq coverage, what percent did rna-seq call
  # sensitivity (y-axis) = given exome snp calls having min alpha rna-seq coverage, what percent did rna-seq call as snps?
  # Specificity: Given RNA-seq calls having min rna-seq coverage, what percent are correct
  # specificity (x-axis) = given rna-seq snp calls having min alpha rna-seq coverage, what percent are correct (exome-called)?
  # Alpha: Min depth of rna-seq coverage

  # Really high sensitivity
  Sensitivity = c()
  Specificity = c()
  # alpha=min depth of rna-seq coverage
  # Minimum coverage defines sample space in primary and secondary evidence (depth)
  # What is the depth in both rna and dna?
  # So at a depth both evidence must meet the minium depth, if not they are dropped (completely)
  # At a certain level the RNASEQ depth will drownd out the DNASeq.
  # Stop, at max( depth( primary ), depth( secondary ) )
 
  # The depth that was used in the sensitivity or specificity
  vi_roc_depth_used = c() 
  vi_roc_depth = c()
  # Which depths has the lowest max
  # This is the range that is used in the ROCs.
  # Looking at depths that are no na
  if( median( df_tab[[ C_I_PRIMARY_DEPTH ]][ vi_primary_calls ], na.rm = TRUE ) < median( df_tab[[ C_I_SECONDARY_DEPTH ]][ vi_secondary_calls ], na.rm = TRUE ) )
  {
    vi_roc_depth = df_tab[[ C_I_PRIMARY_DEPTH ]][ vi_primary_calls ]
  } else {
    vi_roc_depth = df_tab[[ C_I_SECONDARY_DEPTH ]][ vi_secondary_calls ]
  }
  vi_roc_depth = sort( unique( vi_roc_depth ), decreasing=FALSE )
  for( i_alpha_depth in vi_roc_depth )
  {
    # Primary locations with a minimum read evidence and no nas
    vi_primary_with_evidence_at_depth = intersect( which( df_tab[ C_I_PRIMARY_DEPTH ] >= i_alpha_depth ), vi_primary_read_evidence )

    # Secondary locations with a minimum read evidence and no nas
    vi_secondary_with_evidence_at_depth = intersect( which( df_tab[ C_I_SECONDARY_DEPTH ] >= i_alpha_depth ), vi_secondary_read_evidence )

    # Take the locations that are primary calls
    vi_primary_calls_at_depth = intersect( vi_primary_calls, vi_primary_with_evidence_at_depth )
    # Take the locations that are secondary calls
    vi_secondary_calls_at_depth = intersect( vi_secondary_calls, vi_secondary_with_evidence_at_depth )

    # Of the primary calls, which have the same minimum coverage in the secondary assay
    vi_primary_matched_calls_at_depth = intersect( vi_primary_calls_at_depth, vi_secondary_with_evidence_at_depth )
    # Of the secondary calls, which have the same minimum coverage in the primary assay
    vi_secondary_matched_calls_at_depth = intersect( vi_secondary_calls_at_depth, vi_primary_with_evidence_at_depth )

    if( ( length( vi_primary_matched_calls_at_depth ) > 0 ) && ( length( vi_secondary_matched_calls_at_depth ) > 0 ) )
    {
      # Sensitivity is: How many of the primary calls at a depth had a call in the secondary given a minimum coverage in both
      # Primary calls at depth which had the same depth and were called in the secondary / primary calls at depth
      Sensitivity = c( Sensitivity, length( intersect( vi_primary_matched_calls_at_depth, vi_secondary_calls_at_depth ) ) / length( vi_primary_matched_calls_at_depth ) )
      # Specificity (x-axis) = given rna-seq snp calls having min alpha rna-seq coverage, what percent are correct (exome-called)?
      # Secondary calls at depth which had the same depth and were called in the primary / secondary calls at depth
      Specificity = c( Specificity, length( intersect( vi_secondary_matched_calls_at_depth, vi_primary_calls_at_depth ) ) / length( vi_secondary_matched_calls_at_depth ) )
      # Depth for sensitivity and specificity measurements
      vi_roc_depth_used = c( vi_roc_depth_used, i_alpha_depth)
    }
  }
  df_cur = data.frame( sensitivity=Sensitivity, specificity=Specificity )

  # New plot
  # At a given minimum, x = depth, (y) sensitivity vs ( x ) minimum read threshold
  pdf( file.path( str_output_dir, paste( basename( str_file ), "sensitivity_min_read_coverage.pdf", sep = "_" ) ), useDingbats = FALSE )
  plot( vi_roc_depth_used, Sensitivity, main = "Positive Predictive Value vs Min Read Covereage", xlab = "Minimum Read Coverage", ylab = "Positive Predictive Value" ) 
  dev.off()

  pdf( file.path( str_output_dir, paste( basename( str_file ), "ROC2", C_STR_DETAIL_FILE, sep = "_" )), useDingbats = FALSE )
  plot( df_cur$specificity, df_cur$sensitivity )
  plot( df_cur$specificity, df_cur$sensitivity, xlim =c(0,1), ylim=c(0,1) )
  dev.off()

#  ggplot( data = df_cur, aes( x=specificity, y=sensitivity )) + geom_line() + ggtitle("Specificity vs Sensitivity by Depth")
#  gg = ggplot()
#  gg = gg + geom_line( data=df_cur, aes( x=Specificity, y=Sensitivity ) )
#  gg = gg + ggtitle( paste( "Specificity vs Sensitivity at depth (", str_title_key, ")" ) )+ylim(c(0,1.1))+xlim(c(0,1.1))
#  ggsave(file=file.path( str_output_dir, paste( basename( str_file ), "ROC2", C_STR_DETAIL_FILE, sep = "_" )))

  # New plot
  # At a given minimum, x = depth, (y) sensitivity vs ( x ) minimum read threshold

  # Update group infomation
  print("Updating group information")
  ls_group_sensitivity = c( ls_group_sensitivity, Sensitivity )
  ls_group_specificity = c( ls_group_specificity, Specificity )
  ls_group_files = c( ls_group_files, rep( str_file, length( Sensitivity ) ) )
}

# Group ROCs
ggplot( data = data.frame( Sensitivity=ls_group_sensitivity, Specificity=ls_group_specificity, File = ls_group_files ),
        aes( x=Specificity, y=Sensitivity, group=File )) + geom_line() + ggtitle("Specificity vs Sensitivity by Depth")
ggsave( file=file.path( str_output_dir, paste( str_title_key, "group.pdf", sep = "_" ) ) )

ggplot( data = data.frame( Sensitivity=ls_group_sensitivity, Specificity=ls_group_specificity, File = ls_group_files ),
        aes( x=Specificity, y=Sensitivity, group=File )) + geom_line() + ggtitle("Specificity vs Sensitivity by Depth") + ylim(c(0,1))+ xlim(c(0,1))
ggsave( file=file.path( str_output_dir, paste( str_title_key, "group_2.pdf", sep = "_" ) ) )
