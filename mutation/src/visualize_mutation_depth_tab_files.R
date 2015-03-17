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
C_I_HOLD_TRUTH_AT_DEPTH = 2

I_ROC_LINES = 6

# Argument parsing
pArgs = OptionParser( usage ="%prog files1.tab files2.tab" )
pArgs = add_option( pArgs, c("-k","--title_key"),type="character",action="store",dest="str_title_key",default="Primary vs Secondary",help="Key identifying the contrast being visualized (eg \"DNA vs RNA\")")
pArgs = add_option( pArgs, c("-o","--group_output_dir"),type="character",action="store",dest="str_output_dir",default=NA,help="Output directory (required).")
lsArgs = parse_args( pArgs, positional_arguments=TRUE )

func_select_even_roc_depths = function( vi_depths, vi_depths_indicies, i_number_of_selections = length( vi_depths_indicies ) )
{
  # Given a list of depths and indicies indicating which depths are of interest.
  # Return depths that are in order, and even selected throughout the sequence of depths (even in order not value)
  # Will return the number of selections or the max number of unique depths, which ever is smaller
  # Also select representative depths from the data
  # Order the indices by the depth and then evenly subselect a group of indices
  vi_unique_indicies = vi_depths_indicies[ which( ! duplicated( vi_depths[ vi_depths_indicies ] ) ) ]
  vi_ordered_unique_indicies = vi_unique_indicies[ order( vi_depths[ vi_unique_indicies ], decreasing = FALSE )]
  i_max_lines = min( i_number_of_selections, length( vi_ordered_unique_indicies ) )
  vi_selected_depth_indicies = vi_ordered_unique_indicies[ seq(1, length( vi_ordered_unique_indicies ), floor( length( vi_ordered_unique_indicies ) / i_max_lines))]
  vi_selected_depths = vi_depths[ vi_selected_depth_indicies ]
  return( vi_selected_depths )
}

# Calculate membership in different error classes
# i_alpha_depth : the minimum depth that a feature must have to be looked at
# vi_primary_calls : indices of features that are calls in the primary ( truth ) data set
# vi_secondary_calls : indices of features that are calls in the secondary ( evaluate ) data set
# vi_depth : The depths of the features
func_calculate_categories = function( i_alpha_depth, vi_primary_calls, vi_secondary_calls, vi_depth_indicies, vi_depth, i_hold_primary_at = NA, vi_primary_calls_hold = NA )
{
  # Locations of interest for this depth
  # Matched at the depth given a min coverage in both evidence
  vi_features = intersect( which( vi_depth[ vi_depth_indicies] >= i_alpha_depth ), which( !is.na( vi_depth[ vi_depth_indicies] ) ) )
  vi_features = vi_depth_indicies[ vi_features ]
  # The primary calls will either undr the same constraints as the secondary calls (requiring a certain depth) or
  # Will just have a floor to be held at (must atleast be a minimal depth)
  vi_primary_calls_at_depth = intersect( vi_features, vi_primary_calls )
  if( !is.na( vi_primary_calls_hold ) && !is.na( vi_primary_calls_hold))
  {
    vi_primary_calls_at_depth = vi_primary_calls_hold
  }
  vi_secondary_calls_at_depth = intersect( vi_features, vi_secondary_calls )

  # Calculate error classes
  # TP
  vi_TP = intersect( vi_primary_calls_at_depth, vi_secondary_calls_at_depth )
  # FP
  vi_FP = setdiff( vi_secondary_calls_at_depth, vi_primary_calls_at_depth )
  # FN
  vi_FN = setdiff( vi_primary_calls_at_depth, vi_secondary_calls_at_depth )
  if( !is.na( vi_primary_calls_hold ) && !is.na( i_hold_primary_at ) )
  {
    # FN accumulate TP that are less than i_alpha_depth
    # FN stay FN
    # FP are ignored once they are below i_alpha_depth
    # FN that are being dropped are not na, are alteast the hold, less than the alpha, and otherwise a true positive
    vi_FN_by_depth = intersect( vi_secondary_calls, which( !is.na( vi_depth )) )
    vi_FN_by_depth = intersect( vi_FN_by_depth, which( vi_depth >= i_hold_primary_at ))
    vi_FN_by_depth = intersect( vi_FN_by_depth, which( vi_depth < i_alpha_depth ) )
    vi_global_TP = intersect( vi_primary_calls, vi_secondary_calls )
    vi_FN_by_depth = intersect( vi_global_TP, vi_FN_by_depth )
    vi_FN = union( setdiff( vi_primary_calls_at_depth, vi_secondary_calls_at_depth ), vi_FN_by_depth )
  }

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

# Calculate rocs
# vi_primary_calls: Indicies of the truth calls
# vi_secondary_calls: Indicies of calls for the calls being investigated
# vi_depth_indicies: Indicies of depth that will be used
# vi_depth:  The original depth values.
func_calculate_roc_values = function( vi_primary_calls, vi_secondary_calls, vi_depth_indicies, vi_depth )
{
  vi_depth_used = c()
  vi_tp = c()
  vi_fp = c()
  vi_fn = c()
  vi_cumulative_tp_r = c()
  i_cumulative_tp = 0
  vi_cumulative_fp_r = c()
  i_cumulative_fp = 0

  print( "Primary")
  print( vi_primary_calls )
  print( "Secondary" )
  print( vi_secondary_calls )
  print( "vi_depth_indicies" )
  print( vi_depth_indicies )

  # Get the calls given the depth restrictions
  vi_positives = intersect( vi_primary_calls, vi_depth_indicies )
  print("positives")
  print( vi_positives )
  vi_negatives = setdiff( intersect( vi_secondary_calls, vi_depth_indicies ),intersect( vi_primary_calls, vi_depth_indicies ) )
  print("negatives")
  print( vi_negatives )

  # Total positives and negative to use when calculating cumulative rates 
  i_total_positive = length( vi_positives )
  i_total_negative = length( vi_negatives )
  print( "Total pos / neg" )
  print( i_total_positive )
  print( i_total_negative )

  # Calculate per depth
  print( "D")
  print( vi_depth[ vi_depth_indicies ] )
  print( sort( vi_depth[ vi_depth_indicies ], decreasing = FALSE ) )
  for( i_cur_depth in sort( unique( vi_depth[ vi_depth_indicies ] ), decreasing = FALSE ) )
  {
    print("Depth")
    print( i_cur_depth )
    # Select the indicies for the depth
    vi_cur_features = intersect( which( vi_depth == i_cur_depth ), vi_depth_indicies )
    print( "features")
    print( vi_cur_features )
    vi_primary_features = intersect( vi_cur_features, vi_primary_calls )
    vi_secondary_features = intersect( vi_cur_features, vi_secondary_calls )
    print( "Pri / sec features" )
    print( vi_primary_features )
    print( vi_secondary_features )

    # Of these indices which are called
    # In both (TP)
    vi_cur_tp = intersect( vi_primary_features, vi_secondary_features )
    i_cur_tp = length( vi_cur_tp )

    # In Primary only (FN)
    vi_cur_fn = setdiff( vi_primary_features, vi_secondary_features )
    i_cur_fn = length( vi_cur_fn )

    # In Secondary only (FP)
    vi_cur_fp = setdiff( vi_secondary_features, vi_primary_features )
    i_cur_fp = length( vi_cur_fp )

    # Accumulate
    vi_depth_used = c( vi_depth_used, i_cur_depth )
    vi_tp = c( vi_tp, i_cur_tp )
    vi_fp = c( vi_fp, i_cur_fp )
    vi_fn = c( vi_fn, i_cur_fn )
    print( i_cur_tp / i_total_positive )
    print( i_cur_fp / i_total_negative )
    i_cumulative_tp = i_cumulative_tp + ( i_cur_tp / i_total_positive )
    i_cumulative_fp = i_cumulative_fp + ( i_cur_fp / i_total_negative )
    vi_cumulative_tp_r = c( vi_cumulative_tp_r, round( i_cumulative_tp, 3 ) )
    vi_cumulative_fp_r = c( vi_cumulative_fp_r, round( i_cumulative_fp, 3 ) )
    print( vi_cumulative_tp_r )
    print( vi_cumulative_fp_r )
  }
  # Return Depth, TP , FP, FN and cumulative sums
  return( data.frame( Depth=vi_depth_used, TP=vi_tp, FP=vi_fp, FN=vi_fn, TPR=vi_cumulative_tp_r, FPR=vi_cumulative_fp_r ) )
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
  vi_primary_no_depth = union( which( is.na( df_tab[[ C_I_PRIMARY_DEPTH ]] )), which( df_tab[[ C_I_PRIMARY_DEPTH ]] < 1 ))
  ### Primary calls
  vi_primary_calls = which( !is.na(  df_tab[[ C_I_PRIMARY_GT ]] ) )
  ## Secondary
  ### Locations with no depth
  vi_secondary_no_depth = union( which( is.na( df_tab[[ C_I_SECONDARY_DEPTH ]])), which( df_tab[[ C_I_SECONDARY_DEPTH ]] < 1 ))
  ### Calls with or without depth
  vi_secondary_calls = which( !is.na( df_tab[[ C_I_SECONDARY_GT ]] ) )

  ## Primary are MAF or DNA
  ## Secondary are DNA or RNA
  # How many primary calls are not supported by the secondary assay
  # How many secondary calls are not supported by the primary assay
  vi_primary_exclusive = intersect( vi_primary_calls, vi_secondary_no_depth )
  vi_secondary_exclusive = intersect( vi_secondary_calls, vi_primary_no_depth )
  print( "Percent of the primary calls with no evidence in the secondary calls" )
  print( length( vi_primary_exclusive )/length( vi_primary_calls ) * 100 )
  print( "Percent of the secondary calls with no evidence in the primary calls." )
  print( length( vi_secondary_exclusive )/length( vi_secondary_calls ) * 100 )
  print( "Remember this is a problematic view for maf files given they do not have thier own depth files or it is the same as the DNA files." )

  # Which depths has the lowest max per feature
  # This is the range that is used in the ROCs.
  # This can have NAs or 0s in it if an evidence does not have a read depth
  vi_min_depth = apply( df_tab[ c( C_I_PRIMARY_DEPTH, C_I_SECONDARY_DEPTH ) ], MARGIN = 1, min )
  vi_roc_depth_used = c()
  vi_roc_depth_used_at_2 = c()

  # Plot the distribution of the data and get the central tendency
  # Remove NAs and 0 entries from the list.
  vi_min_depth_no_na_indicies = intersect( which( !is.na( vi_min_depth ) ), which( vi_min_depth != 0 ) )
  vi_min_depth_no_na = vi_min_depth[ vi_min_depth_no_na_indicies ]
  i_min_depth_no_na = length( vi_min_depth_no_na )

  # Also select representative depths from the data
  # Order the indices by the depth and then evenly subselect a group of indices
  vi_selected_depths = func_select_even_roc_depths( vi_depths=vi_min_depth, vi_depths_indicies=vi_min_depth_no_na_indicies, i_number_of_selections=I_ROC_LINES )

  # Get one individual distributions without NA and 0 depth features
  vi_truth_depth_no_na = df_tab[[ C_I_PRIMARY_DEPTH ]][ intersect( which( !is.na( df_tab[[ C_I_PRIMARY_DEPTH ]] )), which( df_tab[[ C_I_PRIMARY_DEPTH ]] != 0 ) ) ]
  vi_evaluated_depth_no_na = df_tab[[ C_I_SECONDARY_DEPTH ]][ intersect( which( !is.na( df_tab[[ C_I_SECONDARY_DEPTH ]] )), which( df_tab[[ C_I_SECONDARY_DEPTH ]] != 0 ) ) ]

  # Plot distributions
  pdf( file.path( str_output_dir, paste( basename( str_file ), "depth_distributions.pdf", sep = "_" ) ), useDingbats = FALSE )

  # Get central value for min coverage
  i_center = mean( log( vi_min_depth_no_na, 2 ) )
  i_center_raw = round( 2^i_center )

  # Plot depth distributions no matter the call state
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
                                                   vi_depth_indicies = vi_min_depth_no_na_indicies,
                                                   vi_depth = vi_min_depth )
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

  # Changing bounds
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
    # For optimization plot
    list_categories = func_calculate_categories( i_alpha_depth=i_alpha_depth, vi_primary_calls=vi_primary_calls,
                                                 vi_secondary_calls=vi_secondary_calls, vi_depth_indicies=vi_min_depth_no_na_indicies,
                                                 vi_depth=vi_min_depth )
    # Check to make sure we have at least the min percent of features or break
    if( ( list_categories$Feature_count / i_min_depth_no_na ) < C_I_MIN_PERCENT_FEATURES )
    {
      break
    }

    # Get error rates 
    i_TP = list_categories$TP
    i_FP = list_categories$FP
    i_FN = list_categories$FN
    # Calculating Sensitivity / Specificity Measurements
    if( ( i_TP + i_FN ) > 0 && ( i_TP + i_FP ) > 0 )
    {
      # Sensitivity
      # TP / ( TP + FN )
      i_sensitivity = i_TP / ( i_TP + i_FN )
      Sensitivity = c( Sensitivity, i_sensitivity )
      # TP / ( TP + FP )
      # FDR
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

  # Go through indices for each ROC line
  pdf( file.path( str_output_dir, paste( basename( str_file ), "roc.pdf", sep = "_" ) ), useDingbats = FALSE )
  plot.new()
  vi_roc_ticks = c(0.0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0)
  lines( x = vi_roc_ticks, y = vi_roc_ticks, col = "grey" )
  i_roc_index = 1
  print("Primary calls" )
  print( vi_primary_calls )
  print("Secondary calls" )
  print( vi_secondary_calls )
  vstr_roc_colors = rainbow( length( vi_selected_depths ) )
  for( i_cur_roc_depth in vi_selected_depths )
  {
    df_roc_results = func_calculate_roc_values( vi_primary_calls = vi_primary_calls,
                               vi_secondary_calls = vi_secondary_calls, 
                               vi_depth_indicies = which( vi_min_depth >= i_cur_roc_depth ),
                               vi_depth = vi_min_depth )
    print( "df_roc_results" )
    print( df_roc_results )
    lines( x = c(0,df_roc_results[[ "FPR" ]],1), y = c(0,df_roc_results[[ "TPR" ]],1), col = vstr_roc_colors[ i_roc_index ]  )
    i_roc_index = i_roc_index + 1
    write.table( df_roc_results, file = file.path( str_output_dir, paste( basename( str_file ), "data_roc", i_cur_roc_depth,".txt", sep = "_" ) ) )
#    return( list( Depth=vi_depth_used, TP=vi_tp, FP=vi_fp, FN=vi_fn, TPR=vi_cumulative_tp_r, FPR=vi_cumulative_fp_r ) )
  }
  title( main="TPR vs FPR varying by depth *remeber remove 10%" )
  legend( "bottomright", legend= c( vi_selected_depths, "Random" ), fill = c( vstr_roc_colors, "grey" ), border = "black", title = "Min Depth" )
  vstr_roc_ticks = paste( vi_roc_ticks )
  axis( 1, at=vi_roc_ticks, labels=vstr_roc_ticks )
  mtext( side = 1, "False Positive Rate", line = 2 )
  axis( 2, at=vi_roc_ticks, labels=vstr_roc_ticks ) 
  mtext( side = 2, "True Positive Rate", line = 2 )
  dev.off()

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
  lines( x = c( i_specificity_at_center, 1 ), y = c( i_sensitivity_at_center, i_sensitivity_at_center ), col = "darkgoldenrod1" )
  lines( x = c( i_specificity_at_center, i_specificity_at_center ), y = c( 0, i_sensitivity_at_center ), col = "darkgoldenrod1" )
  plot( df_cur$specificity, df_cur$sensitivity, main = "Sensitivity vs False Discovery Rate", xlab = "False Discovery Rate 1 - ( TP/TP+FP )", ylab = "Sensitivity (TP/TP+FN)", pch = 21, 
                                                xlim = c( 0, max( df_cur$specificity ) ), ylim = c( 0, max( df_cur$sensitivity ) ), col = plt_border, bg = plt_blues )
  lines( x = c( i_specificity_at_center, 1 ), y = c( i_sensitivity_at_center, i_sensitivity_at_center ), col = "darkgoldenrod1" )
  lines( x = c( i_specificity_at_center, i_specificity_at_center ), y = c( 0, i_sensitivity_at_center ), col = "darkgoldenrod1" )
  dev.off()

  # Write measurements to file for troubleshooting
  write.table( df_cur, file = file.path( str_output_dir, paste( basename( str_file ), "data.txt", sep = "_" ) ) )
}
