#!/usr/bin/env Rscript

library(colorRamps )
library( RColorBrewer )
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

C_I_MIN_PERCENT_FEATURES = .1
C_I_INDIVIDUALS_PER_BIN = 5

VI_ROC_DEPTHS = c( 1:10, 20, 30, 40, 50, 60 )
VI_ROC_RNA_DEPTHS = c( 1, 2, 3, 4, 5, 10, 15, 20 )
VI_ROC_TRUTH_MIN_DEPTH = c( 1, 2, 3, 4, 5, 10, 15, 20 )
I_SELECTED_EVAL_MIN_COV = 1
I_SELECTED_TRUTH_MIN_COV = 10

# Argument parsing
pArgs = OptionParser( usage ="%prog files1.tab files2.tab" )
pArgs = add_option( pArgs, c("-k","--title_key"),type="character",action="store",dest="str_title_key",default="Primary vs Secondary",help="Key identifying the contrast being visualized (eg \"DNA vs RNA\")")
pArgs = add_option( pArgs, c("-o","--group_output_dir"),type="character",action="store",dest="str_output_dir",default=NA,help="Output directory (required).")
lsArgs = parse_args( pArgs, positional_arguments=TRUE )

func_plot_roc = function( list_TPR, list_FDR, vi_depths, i_mean_depth, str_pdf_file_name, str_title, str_legend_vary_title  )
{
  pdf( str_pdf_file_name, useDingbats = FALSE )
  plot.new()
  # This is the ROC that cutts off calls less than the depth
  vi_roc_ticks = c(0.0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0)
  lines( x = vi_roc_ticks, y = vi_roc_ticks, col = "grey" )
  # These ROCs select on the secondary data depth
  vstr_roc_colors = rainbow( length( vi_depths ))
  # Plot lines
  vstr_exome_cov = names( list_TPR )
  # Exome cov colors
  vstr_exome_colors = primary.colors( length( vstr_exome_cov ) )
  for( i_exome_cov_i in 1:length( vstr_exome_cov ))
  {
    vi_FDR = list_FDR[[ vstr_exome_cov[ i_exome_cov_i ] ]]
    vi_TPR = list_TPR[[ vstr_exome_cov[ i_exome_cov_i ] ]]
    lines( x=vi_FDR, y=vi_TPR, col=vstr_exome_colors[ i_exome_cov_i ] )
    points( x=vi_FDR, y=vi_TPR, pch=21, col=vstr_exome_colors[ i_exome_cov_i ], bg= vstr_roc_colors )
  }
  title( str_title )
  # Legend Parameters
  vstr_legend_labels_exome = c( vstr_exome_cov, "Random" )
  vstr_legend_labels_rna = c( vi_depths, paste("Mean (",i_mean_depth,")",sep="" ))
  vstr_legend_col_exome = c( vstr_exome_colors, "grey" )
  vstr_legend_fill_rna = c( vstr_roc_colors, "white" )
  vstr_legend_border_rna = c( rep("black", length( vstr_roc_colors )),"white")
  vstr_legend_line_exome = c( rep(1,length(vstr_exome_cov)),1 )
  legend( "bottomright", legend=vstr_legend_labels_rna, fill=vstr_legend_fill_rna, border=vstr_legend_border_rna, title="Vary RNA Cov", cex=.75 )
  legend( "topleft", legend=vstr_legend_labels_exome, col=vstr_legend_col_exome, title=str_legend_vary_title, lty=vstr_legend_line_exome, lwd=2, cex=.75 )
  vstr_roc_ticks = paste( vi_roc_ticks )
  axis( 2, at=vi_roc_ticks, labels=vstr_roc_ticks )
  mtext( side = 2, "True Positive Rate (TP/FP+FN)", line = 2 )
  axis( 1, at=vi_roc_ticks, labels=vstr_roc_ticks ) 
  mtext( side = 1, "False Discovery Rate (FP/FP+TP)", line = 2 )
  dev.off()
}

func_plot_seperate_metrics = function( df_measurements )
{
  # Cut off at 10% features
  i_number_features = df_measurements[[ "TP_min" ]][ 1 ] + df_measurements[[ "FP_min" ]][ 1 ] + df_measurements[[ "FN_min" ]][ 1 ]
  i_cut_off = floor( i_number_features * .1 )
  i_cut_off_index = nrow( df_measurements )
  for( i in 1:nrow( df_measurements ))
  {
    if( df_measurements[[ "TP_min" ]][ i ] + df_measurements[[ "FP_min" ]][ i ] + df_measurements[[ "FN_min" ]][ i ] < i_cut_off )
    {
      i_cut_off_index = i - 1
      break
    }
  }
  # Colors used in plotting
  ## Fill
  i_depth_used = length( df_measurements$Depth[ 1:i_cut_off_index ] )
  i_depth_half = floor( i_depth_used /2 )
  plt_blues = adjustcolor( colorRampPalette(brewer.pal( 9, "Blues" ))( i_depth_used ), alpha.f = 0.5 )

  ## Borders
  str_border_color = "#000000"
  plt_border = sapply( c(round( ( i_depth_half:1 ) / i_depth_half, 2 ),rep(0, i_depth_used-i_depth_half)), function( x ) adjustcolor( str_border_color, alpha.f = x ) )

  # At a given minimum, x = depth, (y) sensitivity vs ( x ) minimum read threshold
  plot( df_measurements$Depth[ 1:i_cut_off_index ], df_measurements[[ "TPR_min" ]][ 1:i_cut_off_index ], main = "TPR vs Min Read Coverage", xlab = paste("Min Coverage (Requiring ",C_I_MIN_PERCENT_FEATURES * 100,"% data)",sep=""), ylab = "TPR", pch = 21, col = plt_border, bg = plt_blues ) 
#  abline( v = i_center_raw, col = "darkgoldenrod1" )

  # At a given minimum, x = depth, (y) postivive predictive value vs ( x ) minimum read threshold
  plot( df_measurements$Depth[ 1:i_cut_off_index ], df_measurements[[ "FDR_min" ]][ 1:i_cut_off_index ], main = "FDR vs Min Read Coverage", xlab = paste("Min Read Coverage (Requiring ",C_I_MIN_PERCENT_FEATURES * 100,"% data)",sep=""), ylab = "FDR", pch = 21, col = plt_border, bg = plt_blues ) 

  # Optimization plot
  plot( df_measurements[[ "FDR_min" ]][ 1:i_cut_off_index ], df_measurements[[ "TPR_min" ]][ 1:i_cut_off_index ], main = paste("TPR vs FDR (Requiring ",C_I_MIN_PERCENT_FEATURES * 100,"% data)",sep=""), xlab = "False Discovery Rate (FP/FP+TP)", ylab = "Sensitivity (TP/TP+FN)", pch = 21, col = plt_border, bg = plt_blues )
  plot( df_measurements[[ "FDR_min" ]][ 1:i_cut_off_index ], df_measurements[[ "TPR_min" ]][ 1:i_cut_off_index ], main = "TPR vs FDR", xlab = "False Discovery Rate (FP/TP+FP)", ylab = "TPR (TP/TP+FN)", pch = 21, col = plt_border, bg = plt_blues )
}

# Calculate rocs
# vi_primary_calls: Indicies of the truth calls
# vi_secondary_calls: Indicies of calls for the calls being investigated
func_calculate_roc_values = function( vi_primary_calls, vi_secondary_calls )
{
  # Total positives and negative to use when calculating cumulative rates 
  i_total_tp = length( intersect( vi_secondary_calls, vi_primary_calls ))
  i_total_fp = length( setdiff( vi_secondary_calls, vi_primary_calls ))
  i_total_fn = length( setdiff( vi_primary_calls, vi_secondary_calls ))
  return( data.frame( TP=i_total_tp, FP=i_total_fp, FN=i_total_fn,
                      TPR=i_total_tp / ( i_total_tp + i_total_fn ),
                      FDR = i_total_fp / ( i_total_tp + i_total_fp ) ) )
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
  df_orig = read.table( str_file )

  # Make sure that the file does not include calls that do not get called in both primary and secondary sources
  # Make sure the coverage is not 0 or NA
  vi_remove_calls = intersect( which( is.na( df_orig[[ C_I_PRIMARY_GT ]]  ) ), which( is.na( df_orig[[ C_I_SECONDARY_GT ]] ) ) )
  vi_remove_calls = union( vi_remove_calls, which( is.na( df_orig[[ C_I_PRIMARY_DEPTH ]] ) ) )
  vi_remove_calls = union( vi_remove_calls, which( is.na( df_orig[[ C_I_SECONDARY_DEPTH ]] ) ) )
  vi_remove_calls = union( vi_remove_calls, which( df_orig[[ C_I_SECONDARY_DEPTH ]] < 1 ) )
  vi_remove_calls = union( vi_remove_calls, which( df_orig[[ C_I_PRIMARY_DEPTH ]] < 1 ) )
  vi_remove_calls = unique( vi_remove_calls )
  if( length( vi_remove_calls ) > 0 )
  {
    df_orig = df_orig[-1*vi_remove_calls,]
  }

  # Mean after minimal filtering
  i_mean_depth=round( mean( df_orig[[ C_I_SECONDARY_DEPTH ]] ),2 )

  # Go through Exome min coverage settings
  # Hold the TPR and FDR lines for each setting in a list
  list_TPR = list()
  list_FDR = list()
  for( i_cur_exome_min_coverage in VI_ROC_TRUTH_MIN_DEPTH )
  {
    print("Exome min coverage")
    print( i_cur_exome_min_coverage )
    # Set the initial problem state.
    # Min Exome coverage
    df_tab = df_orig
    vi_remove_calls = which( df_orig[[ C_I_PRIMARY_DEPTH ]] < i_cur_exome_min_coverage )
    vi_remove_calls = which( df_orig[[ C_I_SECONDARY_DEPTH ]] < I_SELECTED_EVAL_MIN_COV )
    if( length( vi_remove_calls ) > 0 )
    {
      df_tab = df_orig[-1*vi_remove_calls,]
    }
    # Different measurements used later on
    ### Primary calls
    vi_primary_calls = which( !is.na(  df_tab[[ C_I_PRIMARY_GT ]] ) )
    ## Secondary
    ### Calls with or without depth
    vi_secondary_calls = which( !is.na( df_tab[[ C_I_SECONDARY_GT ]] ) )
    # Update the secondary calls, if the genotypes do no match the primary calls
    # Remove from the secondary calls
    vi_common_calls = intersect( vi_secondary_calls, vi_primary_calls )
    vi_incorrect_common_calls = c()
    for( i in vi_common_calls )
    {
      if( df_tab[[ C_I_PRIMARY_GT ]][ i ] != df_tab[[ C_I_SECONDARY_GT ]][ i ] )
      {
        vi_incorrect_common_calls = c( vi_incorrect_common_calls, i )
      }
    }
    print( paste( "Percent of the commonly calls features which do NOT agree in genotype. ",length( vi_incorrect_common_calls )," / ",length( vi_common_calls )," = ", round( length( vi_incorrect_common_calls ) / length( vi_common_calls ), 4 ) * 100, sep="" ) )
    vi_secondary_calls = setdiff( vi_secondary_calls, vi_incorrect_common_calls )
    # Classify calls
    vi_TPR = c()
    vi_FDR = c()
    vi_TP = c()
    vi_FP = c()
    vi_FN = c()
    vi_roc_secondary_calls = vi_secondary_calls
    vi_roc_depths = sort( VI_ROC_DEPTHS, decreasing=FALSE )
    for( i_cur_roc_depth in vi_roc_depths )
    {
      # Secondary calls that d not meet the depth criteria change
      # TP -> FN, FN -> FN, FP -> TN
      # Easiest thing is to remove them from the secondary calls, then this happens automagically
      vi_roc_secondary_calls = intersect( vi_roc_secondary_calls, which( df_tab[[ C_I_SECONDARY_DEPTH ]] >= i_cur_roc_depth ) )
      # Ignore the case of the secondary not at the depth level
      df_roc_results = func_calculate_roc_values( vi_primary_calls = vi_primary_calls, vi_secondary_calls = vi_roc_secondary_calls )
      # Store infor to write.
      vi_TPR = c( vi_TPR, df_roc_results[[ "TPR" ]] )
      vi_FDR = c( vi_FDR, df_roc_results[[ "FDR" ]] )
      vi_TP = c( vi_TP, df_roc_results[[ "TP" ]] )
      vi_FP = c( vi_FP, df_roc_results[[ "FP" ]] )
      vi_FN = c( vi_FN, df_roc_results[[ "FN" ]] )
    }
    # Update list of TPR and FDR
    list_TPR[[paste(i_cur_exome_min_coverage)]] = vi_TPR
    list_FDR[[paste(i_cur_exome_min_coverage)]] = vi_FDR
    # Write data per exome depth
    df_roc = data.frame( TPR=vi_TPR, FDR=vi_FDR, TP=vi_TP, FP=vi_FP, FN=vi_FN )
    rownames( df_roc ) = vi_roc_depths
    write.table( df_roc, file = file.path( str_output_dir, paste( basename( str_file ), paste("data_roc_exome_",i_cur_exome_min_coverage,".txt",sep=""), sep = "_" ) ) )
  }
  # Go through indices for each ROC line, filter at depth in RNASEQ
  str_file_roc_rnaseq = file.path( str_output_dir, paste( basename( str_file ), "roc_rnaseq_depth.pdf", sep = "_" ) )
  func_plot_roc( list_TPR=list_TPR, list_FDR=list_FDR, vi_depths=vi_roc_depths,
                 i_mean_depth=i_mean_depth,
                 str_pdf_file_name=str_file_roc_rnaseq,
                 str_title=paste("TPR vs FDR varying Exome Min Cov ( RNA Seq Min Cov ",I_SELECTED_EVAL_MIN_COV,")",sep=""),
                 str_legend_vary_title = "Min Exome Cov" )

  # Go through RNASEQ min coverage settings
  # Hold the TPR and FDR lines for each setting in a list
  list_TPR = list()
  list_FDR = list()
  for( i_cur_rnaseq_min_coverage in VI_ROC_RNA_DEPTHS )
  {
    print("RNASEQ min coverage")
    print( i_cur_rnaseq_min_coverage)
    # Set the initial problem state.
    # Min Exome coverage
    df_tab = df_orig
    vi_remove_calls = which( df_orig[[ C_I_PRIMARY_DEPTH ]] < I_SELECTED_TRUTH_MIN_COV )
    vi_remove_calls = unique( union( vi_remove_calls, which( df_orig[[ C_I_SECONDARY_DEPTH ]] < i_cur_rnaseq_min_coverage ) ) )
    if( length( vi_remove_calls ) > 0 )
    {
      df_tab = df_orig[-1*vi_remove_calls,]
    }
    # Different measurements used later on
    ### Primary calls
    vi_primary_calls = which( !is.na(  df_tab[[ C_I_PRIMARY_GT ]] ) )
    ## Secondary
    ### Calls with or without depth
    vi_secondary_calls = which( !is.na( df_tab[[ C_I_SECONDARY_GT ]] ) )
    # Update the secondary calls, if the genotypes do no match the primary calls
    # Remove from the secondary calls
    vi_common_calls = intersect( vi_secondary_calls, vi_primary_calls )
    vi_incorrect_common_calls = c()
    for( i in vi_common_calls )
    {
      if( df_tab[[ C_I_PRIMARY_GT ]][ i ] != df_tab[[ C_I_SECONDARY_GT ]][ i ] )
      {
        vi_incorrect_common_calls = c( vi_incorrect_common_calls, i )
      }
    }
    print( paste( "Percent of the commonly calls features which do NOT agree in genotype. ",length( vi_incorrect_common_calls )," / ",length( vi_common_calls )," = ", round( length( vi_incorrect_common_calls ) / length( vi_common_calls ), 4 ) * 100, sep="" ) )
    vi_secondary_calls = setdiff( vi_secondary_calls, vi_incorrect_common_calls )
    # Store calls
    vi_TPR = c()
    vi_FDR = c()
    vi_TP = c()
    vi_FP = c()
    vi_FN = c()
    vi_roc_secondary_calls = vi_secondary_calls
    vi_roc_depths = sort( VI_ROC_DEPTHS, decreasing=FALSE )
    for( i_cur_roc_depth in vi_roc_depths )
    {
      if( i_cur_roc_depth < i_cur_rnaseq_min_coverage )
      {
        # Store infor to write.
        vi_TPR = c( vi_TPR, NA )
        vi_FDR = c( vi_FDR, NA )
        vi_TP = c( vi_TP, NA )
        vi_FP = c( vi_FP, NA )
        vi_FN = c( vi_FN, NA )
      } else {
        # Secondary calls that d not meet the depth criteria change
        # TP -> FN, FN -> FN, FP -> TN
        # Easiest thing is to remove them from the secondary calls, then this happens automagically
        vi_roc_secondary_calls = intersect( vi_roc_secondary_calls, which( df_tab[[ C_I_SECONDARY_DEPTH ]] >= i_cur_roc_depth ) )
        # Ignore the case of the secondary not at the depth level
        df_roc_results = func_calculate_roc_values( vi_primary_calls = vi_primary_calls, vi_secondary_calls = vi_roc_secondary_calls )
        # Store infor to write.
        vi_TPR = c( vi_TPR, df_roc_results[[ "TPR" ]] )
        vi_FDR = c( vi_FDR, df_roc_results[[ "FDR" ]] )
        vi_TP = c( vi_TP, df_roc_results[[ "TP" ]] )
        vi_FP = c( vi_FP, df_roc_results[[ "FP" ]] )
        vi_FN = c( vi_FN, df_roc_results[[ "FN" ]] )
      }
    }
    # Update list of TPR and FDR
    list_TPR[[paste(i_cur_rnaseq_min_coverage)]] = vi_TPR
    list_FDR[[paste(i_cur_rnaseq_min_coverage)]] = vi_FDR
    # Write data per exome depth
    df_roc = data.frame( TPR=vi_TPR, FDR=vi_FDR, TP=vi_TP, FP=vi_FP, FN=vi_FN )
    rownames( df_roc ) = vi_roc_depths
    write.table( df_roc, file = file.path( str_output_dir, paste( basename( str_file ), paste("data_roc_rnaseq_",i_cur_exome_min_coverage,".txt",sep=""), sep = "_" ) ) )
  }
  # Go through indices for each ROC line, filter at depth in RNASEQ
  str_file_roc_rnaseq = file.path( str_output_dir, paste( basename( str_file ), "roc_rnaseq_depth_2_vary_rna.pdf", sep = "_" ) )
  func_plot_roc( list_TPR=list_TPR, list_FDR=list_FDR, vi_depths=vi_roc_depths,
                 i_mean_depth=i_mean_depth,
                 str_pdf_file_name=str_file_roc_rnaseq,
                 str_title = paste("TPR vs FDR (Exome Min Cov ",I_SELECTED_TRUTH_MIN_COV,"/ RNA Seq Varying Cov)",sep=" " ),
                 str_legend_vary_title = "Min RNA Cov")
}
