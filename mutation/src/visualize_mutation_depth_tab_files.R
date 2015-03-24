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
VI_ROC_PRED_MIN_DEPTH = c( 1, 2, 3, 4, 5, 10, 15, 20 )
VI_ROC_TRUTH_MIN_DEPTH = c( 1, 2, 3, 4, 5, 10, 15, 20 )
VI_ROC_TICKS = c(0.0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0)
I_SELECTED_PRED_MIN_COV = 1
I_SELECTED_TRUTH_MIN_COV = 10

# Argument parsing
pArgs = OptionParser( usage ="%prog -o output_dir file1.tab file2.tab" )
pArgs = add_option( pArgs, c("-k","--title_key"),type="character",action="store",dest="str_title_key",default="Primary vs Secondary",help="Key identifying the contrast being visualized (eg \"DNA vs RNA\")")
pArgs = add_option( pArgs, c("-o","--group_output_dir"),type="character",action="store",dest="str_output_dir",default=NA,help="Output directory (required).")
pArgs = add_option( pArgs, c("-t","--measure_transitions" ), type="logical", action="store_true",dest="f_calculate_transitions",default=FALSE,help="Turns on calculating nucleotide transitions, can take time to calculate.")
lsArgs = parse_args( pArgs, positional_arguments=TRUE )

func_plot_roc = function( list_TPR, list_FDR, vi_depths, i_mean_depth, str_pdf_file_name, str_title, str_legend_vary_title  )
{
  # Generate basis for plot
  pdf( str_pdf_file_name, useDingbats = FALSE )
  plot.new()
  lines( x=VI_ROC_TICKS, y=VI_ROC_TICKS, col="grey" )

  # These ROCs select on the secondary data depth
  vstr_roc_colors = rainbow( length( vi_depths ))
  # Plot lines
  vstr_line_names = names( list_TPR )
  # Line colors
  vstr_exome_colors = primary.colors( length( vstr_line_names ) )
  # Plot each line (TPR vs FDR)
  for( i_name in 1:length( vstr_line_names ))
  {
    vi_FDR = list_FDR[[ vstr_line_names[ i_name ] ]]
    vi_TPR = list_TPR[[ vstr_line_names[ i_name ] ]]
    lines( x=vi_FDR, y=vi_TPR, col=vstr_exome_colors[ i_name ] )
    points( x=vi_FDR, y=vi_TPR, pch=21, col=vstr_exome_colors[ i_name ], bg=vstr_roc_colors )
  }
  title( str_title )

  # Legend Parameters
  vstr_legend_labels_exome = c( vstr_line_names, "Random" )
  vstr_legend_labels_rna = c( vi_depths, paste("Mean (",i_mean_depth,")",sep="" ))
  vstr_legend_col_exome = c( vstr_exome_colors, "grey" )
  vstr_legend_fill_rna = c( vstr_roc_colors, "white" )
  vstr_legend_border_rna = c( rep("black", length( vstr_roc_colors )),"white")
  vstr_legend_line_exome = c( rep(1,length(vstr_line_names)),1 )
  legend( "bottomright", legend=vstr_legend_labels_rna, fill=vstr_legend_fill_rna, border=vstr_legend_border_rna, title="Vary RNA Cov", cex=.75 )
  legend( "topleft", legend=vstr_legend_labels_exome, col=vstr_legend_col_exome, title=str_legend_vary_title, lty=vstr_legend_line_exome, lwd=2, cex=.75 )
  vstr_roc_ticks=paste( VI_ROC_TICKS )
  axis( 2, at=VI_ROC_TICKS, labels=vstr_roc_ticks )
  mtext( side=2, "True Positive Rate (TP/FP+FN)", line = 2 )
  axis( 1, at=VI_ROC_TICKS, labels=vstr_roc_ticks ) 
  mtext( side=1, "False Discovery Rate (FP/FP+TP)", line = 2 )
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
  list_return = list()

  # Total positives and negative to use when calculating cumulative rates 
  vi_tp = intersect( vi_secondary_calls, vi_primary_calls )
  vi_fp = setdiff( vi_secondary_calls, vi_primary_calls )
  vi_fn = setdiff( vi_primary_calls, vi_secondary_calls )

  i_total_tp = length( vi_tp )
  i_total_fp = length( vi_fp )
  i_total_fn = length( vi_fn )
  
  list_return[["TP"]] = i_total_tp
  list_return[["FP"]] = i_total_fp
  list_return[["FN"]] = i_total_fn
  list_return[["TPR"]] = ( i_total_tp / ( i_total_tp + i_total_fn ) )
  list_return[["FDR"]] = ( i_total_fp / ( i_total_tp + i_total_fp ) )
  list_return[["TP_indices"]] = vi_tp
  list_return[["FP_indices"]] = vi_fp
  list_return[["FN_indices"]] = vi_fn
  return( list_return )
}

func_filter_unequal_genotypes = function( df_data, vi_indicies )
{
  vi_incorrect_calls = c()
  for( i in vi_indicies )
  {
    if( df_data[[ C_I_PRIMARY_GT ]][ i ] != df_data[[ C_I_SECONDARY_GT ]][ i ] )
    {
      vi_incorrect_calls = c( vi_incorrect_calls, i )
    }
  }
  return( vi_incorrect_calls )
}

func_vary_coverage_and_measure_classes = function( vi_vary_truth, vi_vary_prediction, str_file_base_name, str_plot_title, str_legend_title, f_vary_truth, f_calculate_transitions = TRUE ) 
{
  # Go through Exome min coverage settings
  # Hold the TPR and FDR lines for each setting in a list
  list_TPR = list()
  list_FDR = list()
  for( i_cur_index in 1:length(vi_vary_truth) )
  {
    i_cur_truth_coverage = vi_vary_truth[ i_cur_index ]
    i_cur_pred_coverage = vi_vary_prediction[ i_cur_index ]
    # Set the initial problem state.
    # Min Exome coverage
    df_tab = df_orig
    vi_remove_calls = which( df_tab[[ C_I_PRIMARY_DEPTH ]] < i_cur_truth_coverage )
    vi_remove_calls = unique( union( vi_remove_calls, which( df_tab[[ C_I_SECONDARY_DEPTH ]] < i_cur_pred_coverage ) ) )
    if( length( vi_remove_calls ) > 0 )
    {
      df_tab = df_tab[-1*vi_remove_calls,]
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
    vi_incorrect_common_calls = func_filter_unequal_genotypes( df_data=df_tab, vi_indicies=vi_common_calls )
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
    vNA_init = rep(NA,length( vi_roc_depths ))
    df_roc = data.frame( TPR=vNA_init, FDR=vNA_init, TP=vNA_init, FP=vNA_init, FN=vNA_init,
                           FP_AG=vNA_init, FP_AC=vNA_init, FP_AT=vNA_init, FP_GC=vNA_init, FP_GA=vNA_init, FP_GT=vNA_init, FP_CG=vNA_init, FP_CT=vNA_init, FP_CA=vNA_init,
                           TP_AG=vNA_init, TP_AC=vNA_init, TP_AT=vNA_init, TP_GC=vNA_init, TP_GA=vNA_init, TP_GT=vNA_init, TP_CG=vNA_init, TP_CT=vNA_init, TP_CA=vNA_init )
    rownames( df_roc ) = vi_roc_depths
    f_measure_transitions = f_calculate_transitions
    for( i_cur_roc_depth in vi_roc_depths )
    {
      print( paste( "Depth ",i_cur_roc_depth, sep = " " ) )
      if( i_cur_roc_depth < i_cur_pred_coverage )
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
        ls_roc_results = func_calculate_roc_values( vi_primary_calls=vi_primary_calls, vi_secondary_calls=vi_roc_secondary_calls )
        # Store error class to write.
        vi_TPR = c( vi_TPR, ls_roc_results[[ "TPR" ]] )
        vi_FDR = c( vi_FDR, ls_roc_results[[ "FDR" ]] )
        vi_TP = c( vi_TP, ls_roc_results[[ "TP" ]] )
        vi_FP = c( vi_FP, ls_roc_results[[ "FP" ]] )
        vi_FN = c( vi_FN, ls_roc_results[[ "FN" ]] )
        if( f_measure_transitions )
        {
          # Used to count nucleotide transitions in errors
          ltrans_fp = list( "A->G"=0, "A->C"=0,"A->T"=0, "G->C"=0, "G->A"=0, "G->T"=0, "C->G"=0, "C->T"=0, "C->A"=0 )
          ltrans_tp = list( "A->G"=0, "A->C"=0,"A->T"=0, "G->C"=0, "G->A"=0, "G->T"=0, "C->G"=0, "C->T"=0, "C->A"=0 )
          # Store transitions
          # Get transition in the prediction
          for( i_TP_index in ls_roc_results[[ "TP_indices" ]] )
          {
            df_row = df_tab[ i_TP_index, ]
            str_ref = as.character( df_row[[ C_I_SECONDARY_REF ]] )
            vstr_genotype = unlist( strsplit( as.character( df_row[[ C_I_SECONDARY_GT ]] ), "/" ))
            for( str_gt in vstr_genotype )
            {
              if( ! str_gt == str_ref )
              {
                str_trans = paste( str_ref, str_gt, sep = "->" )
                ltrans_tp[[ str_trans ]] = ltrans_tp[[ str_trans ]] + 1
              }
            }
          }
          for( i_FP_index in ls_roc_results[[ "FP_indices" ]] )
          {
            df_row = df_tab[ i_FP_index, ]
            str_ref = as.character( df_row[[ C_I_SECONDARY_REF ]] )
            vstr_genotype = unlist( strsplit( as.character(df_row[[ C_I_SECONDARY_GT ]]), "/" ) )
            for( str_gt in vstr_genotype )
            {
              if( ! str_gt == str_ref )
              {
                str_trans = paste( str_ref, str_gt, sep = "->" )
                ltrans_fp[[ str_trans ]] = ltrans_fp[[ str_trans ]] + 1
              }
            }
          }
          df_roc[[ i_cur_roc_depth, "FP_AC" ]] = ltrans_fp[[ "A->C" ]]
          df_roc[[ i_cur_roc_depth, "FP_AG" ]] = ltrans_fp[[ "A->G" ]]
          df_roc[[ i_cur_roc_depth, "FP_AT" ]] = ltrans_fp[[ "A->T" ]]
          df_roc[[ i_cur_roc_depth, "FP_GC" ]] = ltrans_fp[[ "G->C" ]]
          df_roc[[ i_cur_roc_depth, "FP_GT" ]] = ltrans_fp[[ "G->T" ]]
          df_roc[[ i_cur_roc_depth, "FP_GA" ]] = ltrans_fp[[ "G->A" ]]
          df_roc[[ i_cur_roc_depth, "FP_CT" ]] = ltrans_fp[[ "C->T" ]]
          df_roc[[ i_cur_roc_depth, "FP_CA" ]] = ltrans_fp[[ "C->A" ]]
          df_roc[[ i_cur_roc_depth, "FP_CG" ]] = ltrans_fp[[ "C->G" ]]
          df_roc[[ i_cur_roc_depth, "TP_AC" ]] = ltrans_tp[[ "A->C" ]]
          df_roc[[ i_cur_roc_depth, "TP_AG" ]] = ltrans_tp[[ "A->G" ]]
          df_roc[[ i_cur_roc_depth, "TP_AT" ]] = ltrans_tp[[ "A->T" ]]
          df_roc[[ i_cur_roc_depth, "TP_GC" ]] = ltrans_tp[[ "G->C" ]]
          df_roc[[ i_cur_roc_depth, "TP_GT" ]] = ltrans_tp[[ "G->T" ]]
          df_roc[[ i_cur_roc_depth, "TP_GA" ]] = ltrans_tp[[ "G->A" ]]
          df_roc[[ i_cur_roc_depth, "TP_CT" ]] = ltrans_tp[[ "C->T" ]]
          df_roc[[ i_cur_roc_depth, "TP_CA" ]] = ltrans_tp[[ "C->A" ]]
          df_roc[[ i_cur_roc_depth, "TP_CG" ]] = ltrans_tp[[ "C->G" ]]
          f_measure_transitions = FALSE
        }
      }
    }
    # Update list of TPR and FDR
    if( f_vary_truth )
    {
      list_TPR[[ paste(i_cur_truth_coverage) ]] = vi_TPR
      list_FDR[[ paste(i_cur_truth_coverage) ]] = vi_FDR
    } else {
      list_TPR[[ paste(i_cur_pred_coverage) ]] = vi_TPR
      list_FDR[[ paste(i_cur_pred_coverage) ]] = vi_FDR
    }
    # Write data per exome depth
    df_roc[[ "TPR" ]] = vi_TPR
    df_roc[[ "FDR" ]] = vi_FDR
    df_roc[[ "TP" ]] = vi_TP
    df_roc[[ "FP" ]] = vi_FP
    df_roc[[ "FN" ]] = vi_FN
    write.table( df_roc, file = file.path( str_output_dir, paste( str_file_base_name, paste("data_roc_vary_truth_",i_cur_truth_coverage,"_pred_",i_cur_pred_coverage,".txt",sep=""), sep = "_" ) ) )
  }

  # Go through indices for each ROC line, filter at depth in RNASEQ
  str_file_roc_rnaseq = file.path( str_output_dir, paste( str_file_base_name, paste("data_roc_vary_truth_",i_cur_truth_coverage,"_pred_",i_cur_pred_coverage,".pdf",sep=""), sep = "_" ) )
  func_plot_roc( list_TPR=list_TPR, list_FDR=list_FDR, vi_depths=vi_roc_depths,
                 i_mean_depth=i_mean_depth,
                 str_pdf_file_name=str_file_roc_rnaseq,
                 str_title=str_plot_title,
                 str_legend_vary_title=str_legend_title )
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
  vi_remove_calls = intersect( which( is.na( df_orig[[ C_I_PRIMARY_GT ]] ) ), which( is.na( df_orig[[ C_I_SECONDARY_GT ]] ) ) )
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

  # Vary the problem space holding the pred min coverage at 1 and then varying the truth sample space
  # Then look at min coverage for the pred in the resulting sample space
  func_vary_coverage_and_measure_classes( vi_vary_truth=VI_ROC_TRUTH_MIN_DEPTH,
                                          vi_vary_prediction=rep(I_SELECTED_PRED_MIN_COV,length(VI_ROC_TRUTH_MIN_DEPTH)),
                                          str_file_base_name=basename( str_file ), 
                                          str_plot_title=paste("TPR vs FDR varying Min Cov ( RNA Seq Min Cov ",I_SELECTED_PRED_MIN_COV," )",sep=""),
                                          str_legend_title="Min Exome Cov", f_vary_truth=TRUE, f_calculate_transitions=lsArgs$options$f_calculate_transitions ) 

  # Vary the problem space holding truth at 10 min coverage and vary pred min coverage
  # then investigate each pred min coverage in the resulting sample space
  func_vary_coverage_and_measure_classes( vi_vary_truth=rep(I_SELECTED_TRUTH_MIN_COV,length(VI_ROC_PRED_MIN_DEPTH)),
                                          vi_vary_prediction=VI_ROC_PRED_MIN_DEPTH,
                                          str_file_base_name=basename( str_file ),
                                          str_plot_title=paste("TP vs FDR varying Min Cov ( DNA Seq Min Cov ",I_SELECTED_TRUTH_MIN_COV," )",sep=""),
                                          str_legend_title="Min RNA Seq Cov", f_vary_truth=FALSE, f_calculate_transitions = lsArgs$options$f_calculate_transitions )
}
