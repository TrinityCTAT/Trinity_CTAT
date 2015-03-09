#!/usr/bin/env Rscript

library( RColorBrewer )
library( ggplot2 )
library( optparse )
library( vioplot )

C_I_PRIMARY_POS = 1
C_I_PRIMARY_REF = 2
C_I_PRIMARY_GT = 3
C_I_PRIMARY_DEPTH = 4
C_I_SECONDARY_POS = 5
C_I_SECONDARY_REF = 6
C_I_SECONDARY_GT = 7
C_I_SECONDARY_DEPTH = 8

C_STR_DETAIL_FILE = "detail_validation.pdf"

C_I_MIN_PERCENT_FEATURES = .1
C_I_INDIVIDUALS_PER_BIN = 5

# Argument parsing
pArgs = OptionParser( usage ="%prog files1.tab files2.tab" )
pArgs = add_option( pArgs, c("-k","--title_key"),type="character",action="store",dest="str_title_key",default="Primary vs Secondary",help="Key identifying the contrast being visualized (eg \"DNA vs RNA\")")
pArgs = add_option( pArgs, c("-o","--group_output_dir"),type="character",action="store",dest="str_output_dir",default=NA,help="Output directory (required).")
lsArgs = parse_args( pArgs, positional_arguments=TRUE )

# Calculate membership in different error classes
# i_alpha_depth : the minimum depth that a feature must have to be looked at
# vi_primary_calls : indices of features that are calls in the primary ( truth ) data set
# vi_secondary_calls : indices of features that are calls in the secondary ( evaluate ) data set
# vi_depth : The depths of the features
func_calculate_categories = function( i_alpha_depth, vi_primary_calls, vi_secondary_calls, vi_depth )
{
  # Locations of interest for this depth
  # Matched at the depth given a min coverage in both evidence
  vi_features = intersect( which( vi_depth >= i_alpha_depth ), which( !is.na( vi_depth ) ) )
  vi_primary_calls_at_depth = intersect( vi_features, vi_primary_calls )
  vi_secondary_calls_at_depth = intersect( vi_features, vi_secondary_calls )

  # Calculate error classes
  # TP
  vi_TP = intersect( vi_primary_calls_at_depth, vi_secondary_calls_at_depth )
  # FP
  vi_FP = setdiff( vi_secondary_calls_at_depth, vi_primary_calls_at_depth )
  # FN
  vi_FN = setdiff( vi_primary_calls_at_depth, vi_secondary_calls_at_depth )

  # Count class membership
  # Calculate depth distributions per error class
  # Add to return list
  return( list( TP = length( vi_TP ), 
                FP = length( vi_FP ),
                FN = length( vi_FN ),
                Truth_depths = vi_depth[ vi_primary_calls_at_depth ],
                TP_depths = vi_depth[ vi_TP ],
                FP_depths = vi_depth[ vi_FP ],
                FN_depths = vi_depth[ vi_FN ],
                Calls_primary_indices = vi_primary_calls_at_depth,
                Calls_secondary_indices = vi_secondary_calls_at_depth,
                Feature_count = length( vi_features )))
}

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

  # Different measurements used later on
  ## Primary
  ### Locations with no depth
  vi_primary_no_depth = which( is.na( df_tab[[ C_I_PRIMARY_DEPTH ]] ))
  ### Locations with depth
#  vi_primary_read_evidence = which( !is.na( df_tab[[ C_I_PRIMARY_DEPTH ]] ))
  ### Primary calls
  vi_primary_calls = which( !is.na(  df_tab[[ C_I_PRIMARY_GT ]] ) )
  ## Secondary
  ### Locations with no depth
  vi_secondary_no_depth = which( is.na( df_tab[[ C_I_SECONDARY_DEPTH ]] ))
  ### Locations with depth
#  vi_secondary_read_evidence = which( !is.na( df_tab[[ C_I_SECONDARY_DEPTH ]] ))
#  i_count_secondary_read_evidence = length( vi_secondary_read_evidence )
  ### Calls with or without depth
  vi_secondary_calls = which( !is.na( df_tab[[ C_I_SECONDARY_GT ]] ) )
  ### Calls with depth
#  vi_secondary_calls_with_evidence = intersect( vi_secondary_read_evidence, vi_secondary_calls )
#  i_count_secondary_calls_with_evidence = length( vi_secondary_calls_with_evidence )
  ## Both
#  vi_agreeing_calls = intersect( vi_primary_calls, vi_secondary_calls )
  ### Correct secondary calls (defined as secondary agreeing with primary)
#  vi_secondary_calls_agreeing_primary = intersect( vi_secondary_calls_with_evidence, vi_agreeing_calls )
#  i_count_secondary_calls_agreeing_primary = length( vi_secondary_calls_agreeing_primary )
  
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

  # Which depths has the lowest max per feature
  # This is the range that is used in the ROCs.
  # This can have NAs in it if an evidence does not have a read depth
  vi_min_depth = apply( df_tab[ c( C_I_PRIMARY_DEPTH, C_I_SECONDARY_DEPTH ) ], MARGIN = 1, min )
  vi_roc_depth_used = c()

  # Plot the distribution of the data and get the central tendency
  # Remove NAs and 0 entries from the list.
  vi_min_depth_no_na = vi_min_depth[ intersect( which( !is.na( vi_min_depth ) ), which( vi_min_depth != 0 ) ) ]
  vi_truth_depth_no_na = df_tab[[ C_I_PRIMARY_DEPTH ]][ which( !is.na( df_tab[[ C_I_PRIMARY_DEPTH ]] )) ]
  vi_evaluated_depth_no_na = df_tab[[ C_I_SECONDARY_DEPTH ]][ which( !is.na( df_tab[[ C_I_SECONDARY_DEPTH ]] )) ]
  i_min_depth_no_na = length( vi_min_depth_no_na )

  pdf( file.path( str_output_dir, paste( basename( str_file ), "depth_distributions.pdf", sep = "_" ) ), useDingbats = FALSE )

  # Get central value for min coverage
  i_center = mean( log( vi_min_depth_no_na, 2 ) )
  i_center_raw = round( 2^i_center )

  hist( log( vi_min_depth_no_na, 2 ), main = "Distribution of Min Coverage (Log2)", xlab = "Min Coverage", breaks = round( length( vi_min_depth_no_na ) / C_I_INDIVIDUALS_PER_BIN ) )
  abline( v = i_center, col = "RED", lwd = 2 )

  hist( log( vi_truth_depth_no_na, 2 ), main = "Distribution of Depths for the Truth Data Set (Log2)", xlab = "Truth Depth", breaks = round( length( vi_min_depth_no_na ) / C_I_INDIVIDUALS_PER_BIN ) )

  hist( log( vi_evaluated_depth_no_na, 2 ), main = "Distribution of Depths for the Evaluated Data Set (Log2)", xlab = "Evaluated Depth", breaks = round(length( vi_evaluated_depth_no_na ) / C_I_INDIVIDUALS_PER_BIN ))
  dev.off()

  # At any depth what are the TP / FP / FN
  vstr_bar_names = c("Truth","TP","FP","FN")
  vstr_bar_colors = c("Grey","aquamarine","darkorchid1","darkgoldenrod1")
  list_categories_all = func_calculate_categories( i_alpha_depth = -1, vi_primary_calls = vi_primary_calls, 
                                                   vi_secondary_calls = vi_secondary_calls, 
                                                   vi_depth = vi_min_depth_no_na )
  vi_bar_values = c( length( list_categories_all$Truth_depths ), list_categories_all$TP, list_categories_all$FP, list_categories_all$FN )
  pdf( file.path( str_output_dir, paste( basename( str_file ), "raw_class_distributions", C_STR_DETAIL_FILE, sep = "_" )), useDingbats = FALSE )
  barplot( vi_bar_values, names.arg = paste( vstr_bar_names," (", vi_bar_values ,")" ), main = "Proportions of Call Classes", xlab = "Classes",
           col = vstr_bar_colors, border = "grey", legend = vstr_bar_names, beside = TRUE )

  # How does read depth associate with secondary assays being called
  # What are the distributions of read depth for good and bad calls
  vioplot( list_categories_all$Truth_depths, list_categories_all$TP_depths, list_categories_all$FP_depths, list_categories_all$FN_depths,
           names = paste( vstr_bar_names," (", vi_bar_values ,")" ), col = vstr_bar_colors)
  title( "Min Coverage by Class" )

  dev.off()

  Sensitivity = c()
  Specificity = c()
  TP = c()
  FP = c()
  FN = c()

  i_sensitivity_at_center = NA
  i_specificity_at_center = NA
  # Loop through each depths and calculate error rates
  for( i_alpha_depth in 1:max( vi_min_depth_no_na, na.rm = TRUE ) )
  {
    list_categories = func_calculate_categories( i_alpha_depth=i_alpha_depth, vi_primary_calls=vi_primary_calls,
                                            vi_secondary_calls=vi_secondary_calls, vi_min_depth_no_na )

    # Check to make sure we have at least the min percent of features or break
    if( ( list_categories$Feature_count / i_min_depth_no_na ) < C_I_MIN_PERCENT_FEATURES )
    {
      break
    }
 
    i_TP = list_categories$TP
    i_FP = list_categories$FP
    i_FN = list_categories$FN

    if( ( i_TP + i_FN ) > 0 && ( i_TP + i_FP ) > 0 )
    {
      # Sensitivity
      # TP / ( TP + FN )
      i_sensitivity = i_TP / ( i_TP + i_FN )
      Sensitivity = c( Sensitivity, i_sensitivity )
      # TP / ( TP + FP )
      # Positive Predictive value
      i_specificity = 1 - ( i_TP / ( i_TP + i_FP ))
      Specificity = c( Specificity, i_specificity )
      # Record the depth measurement
      vi_roc_depth_used = c( vi_roc_depth_used, i_alpha_depth )
      # Record error classes at each depth
      TP = c( TP, i_TP )
      FP = c( FP, i_FP )
      FN = c( FN, i_FN )

      if( i_alpha_depth == i_center_raw )
      {
        i_sensitivity_at_center = i_sensitivity
        i_specificity_at_center = i_specificity
      }
    }
  }
  df_cur = data.frame( sensitivity=Sensitivity, specificity=Specificity, depth=vi_roc_depth_used, tp=TP, fp=FP, fn=FN )

  # Colors used in plotting
  ## Fill
  i_depth_used = length( vi_roc_depth_used )
  i_depth_half = floor( i_depth_used /2 )
  plt_blues = adjustcolor( colorRampPalette(brewer.pal( 9, "Blues" ))( i_depth_used ), alpha.f = 0.5 )
  ## Borders
  str_border_color = "#000000"
  plt_border = sapply( c(round( ( i_depth_half:1 ) / i_depth_half, 2 ),rep(0, i_depth_used-i_depth_half)), function( x ) adjustcolor( str_border_color, alpha.f = x ) )

  # At a given minimum, x = depth, (y) sensitivity vs ( x ) minimum read threshold
  pdf( file.path( str_output_dir, paste( basename( str_file ), "sensitivity_min_read_coverage_norm.pdf", sep = "_" ) ), useDingbats = FALSE )
  plot( df_cur$depth, df_cur$sensitivity, main = "Sensitivity vs Min Read Coverage", xlab = paste("Minimum Read Coverage (Requiring ",C_I_MIN_PERCENT_FEATURES * 100,"% data)",sep=""), ylab = "Sensitivity", pch = 21, col = plt_border, bg = plt_blues ) 
  abline( v = i_center_raw, col = "darkgoldenrod1" )
  dev.off()

  # At a given minimum, x = depth, (y) postivive predictive value vs ( x ) minimum read threshold
  pdf( file.path( str_output_dir, paste( basename( str_file ), "fdr_min_read_coverage_norm.pdf", sep = "_" ) ), useDingbats = FALSE )
  plot( df_cur$depth, df_cur$specificity, main = "False Discovery Rate vs Min Read Coverage", xlab = paste("Minimum Read Coverage (Requiring ",C_I_MIN_PERCENT_FEATURES * 100,"% data)",sep=""), ylab = "False Discovery Rate", pch = 21, col = plt_border, bg = plt_blues ) 
  abline( v = i_center_raw, col = "darkgoldenrod1" )
  dev.off()

  # Optimization plot
  pdf( file.path( str_output_dir, paste( basename( str_file ), "optimize", C_STR_DETAIL_FILE, sep = "_" )), useDingbats = FALSE )
  plot( df_cur$specificity, df_cur$sensitivity, main = paste("Sensitivity vs False Discovery Rate (Requiring ",C_I_MIN_PERCENT_FEATURES * 100,"% data)",sep=""), xlab = "False Discovery Rate 1 - ( TP/TP+FP )", ylab = "Sensitivity (TP/TP+FN)", pch = 21, col = plt_border, bg = plt_blues )
  lines( x = c( i_specificity_at_center, 1 ), y = c( i_sensitivity_at_center, i_sensitivity_at_center), col = "darkgoldenrod1" )
  lines( x = c( i_specificity_at_center, i_specificity_at_center ), y = c( 0, i_sensitivity_at_center), col = "darkgoldenrod1" )
  plot( df_cur$specificity, df_cur$sensitivity, main = "Sensitivity vs False Discovery Rate", xlab = "False Discovery Rate 1 - ( TP/TP+FP )", ylab = "Sensitivity (TP/TP+FN)", pch = 21, 
                                                xlim = c( 0, max( df_cur$specificity ) ), ylim = c( 0, max( df_cur$sensitivity ) ), col = plt_border, bg = plt_blues )
  lines( x = c( i_specificity_at_center, 1 ), y = c( i_sensitivity_at_center, i_sensitivity_at_center), col = "darkgoldenrod1" )
  lines( x = c( i_specificity_at_center, i_specificity_at_center ), y = c( 0, i_sensitivity_at_center), col = "darkgoldenrod1" )
  dev.off()

  # Write measurements to file for troubleshooting
  write.table( df_cur, file = file.path( str_output_dir, paste( basename( str_file ), "data.txt", sep = "_" ) ) )
}
