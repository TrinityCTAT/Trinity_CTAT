#!/usr/bin/env Rscript

library( colorRamps )
library( RColorBrewer )
library( ggplot2 )
library( optparse )

# TAB file columns indices
C_I_PRIMARY_POS = 1
C_I_PRIMARY_REF = 2
C_I_PRIMARY_GT = 3
C_I_PRIMARY_DEPTH = 4
C_I_SECONDARY_POS = 5
C_I_SECONDARY_REF = 6
C_I_SECONDARY_GT = 7
C_I_SECONDARY_DEPTH = 8

#C_STR_DETAIL_FILE = "detail_validation.pdf"

#C_I_MIN_PERCENT_FEATURES = .1
#C_I_INDIVIDUALS_PER_BIN = 5

# Settings for plotting ROCs
VI_ROC_DEPTHS = c( 1:10, 20, 30, 40, 50, 60 )
VI_ROC_PRED_MIN_DEPTH = c( 1, 3, 5, 10, 20 )
VI_ROC_TRUTH_MIN_DEPTH = c( 1, 3, 5, 10, 20 )
VI_ROC_TICKS = c( 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 )
I_SELECTED_PRED_MIN_COV = 1
I_SELECTED_TRUTH_MIN_COV = 10

# Argument parsing
pArgs = OptionParser( usage ="%prog -o output_dir file1.tab file2.tab" )
pArgs = add_option( pArgs, c( "-c","--compare" ), type="character", action="store", dest="str_compare_file", default=NULL, help="Compares each of the given output files with the tab field given here as a reference. Used in side-by-side plots of two methods.")
pArgs = add_option( pArgs, c( "-k","--title_key" ), type="character", action="store", dest="str_title_key", default="Primary vs Secondary", help="Key identifying the contrast being visualized (eg \"DNA vs RNA\")")
pArgs = add_option( pArgs, c( "-o","--group_output_dir" ), type="character", action="store", dest="str_output_dir", default=NULL, help="Output directory (required).")
pArgs = add_option( pArgs, c( "--method" ), type="character", action="store", dest="str_method_name", default=NULL, help="The name of the method being evaluated (Should match the input file.")
pArgs = add_option( pArgs, c( "--method_compare" ), type="character", action="store", dest="str_method_name_compare", default=NULL, help="The name of the method that is used for comparison (should match the --compare file).")
pArgs = add_option( pArgs, c( "-s", "--sample_tag" ), type="character", action="store", dest="str_sample_key", default=NULL, help="Adds this tag to the output pdf name." )
pArgs = add_option( pArgs, c( "--serial_plots" ), type="logical", action="store_true", dest="f_make_serial_plots", default=FALSE, help="After the sample space is defined, additionally plots each depth as a seperate plot.")
lsArgs = parse_args( pArgs, positional_arguments=TRUE )

func_plot_roc = function( list_TPR, list_FDR, vi_depths, i_mean_depth, str_pdf_file_name, str_title, str_legend_vary_title,
                          list_TPR_compare=NULL, list_FDR_compare=NULL, vi_depths_compare=NULL, str_method_name="",
                          str_method_name_compare="" )
{
  # Plots both ROCs and Optimization plots (coverage by metric)

  # Check parameters for comparison mode.
  if( ! is.null( vi_depths_compare ) )
  {
    if( (! length( vi_depths ) == length( vi_depths_compare ) ) || (! all( vi_depths==vi_depths_compare ) ) )
    {
      print( "Depths are not the same and so the two methods can not be compared" )
      print( vi_depths )
      print( vi_depths_compare )
      return( FALSE )
    }
  }
  if( length( list_TPR_compare ) > 0 || length( list_FDR_compare ) > 0 )
  {
    if( is.null( str_method_name ) || is.null( str_method_name_compare ) )
    {
      print( "When comparing please make sure to include a method name and comparison method name" )
      return( FALSE )
    }
  }

  # Make plot for 
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

  # Plot comparison reference first if given
  if( ! is.null( vi_depths_compare ) )
  {
    for( i_name_compare in 1:length( vstr_line_names ) )
    {
      vi_FDR_compare = list_FDR_compare[[ vstr_line_names[ i_name_compare ] ]]
      vi_TPR_compare = list_TPR_compare[[ vstr_line_names[ i_name_compare ] ]]
      lines( x=vi_FDR_compare, y=vi_TPR_compare, col=vstr_exome_colors[ i_name_compare ], lty = 3 )
      points( x=vi_FDR_compare, y=vi_TPR_compare, pch=24, col=vstr_exome_colors[ i_name_compare ], bg=vstr_roc_colors )
    }
  }
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
  vstr_legend_shape_exome = c( rep( NA, length( vstr_line_names )), NA )
  # Update legend for compare information
  if( ! is.null( vi_depths_compare ) )
  {
    vstr_legend_labels_exome = c( vstr_legend_labels_exome, str_method_name, str_method_name_compare )
    vstr_legend_col_exome = c( vstr_legend_col_exome, "black", "black" )
    vstr_legend_line_exome = c( vstr_legend_line_exome, 1, 3 )
    vstr_legend_shape_exome = c( vstr_legend_shape_exome, 16, 17 )
  }
  legend( "bottomright", legend=vstr_legend_labels_rna, fill=vstr_legend_fill_rna, border=vstr_legend_border_rna, title="Vary RNA Cov", cex=.75 )
  legend( "topright", legend=vstr_legend_labels_exome, col=vstr_legend_col_exome, title=str_legend_vary_title, lty=vstr_legend_line_exome, lwd=2, cex=.75, pch=vstr_legend_shape_exome )
  vstr_roc_ticks=paste( VI_ROC_TICKS )
  axis( 2, at=VI_ROC_TICKS, labels=vstr_roc_ticks )
  mtext( side=2, "True Positive Rate (TP/TP+FN)", line = 2 )
  axis( 1, at=VI_ROC_TICKS, labels=vstr_roc_ticks ) 
  mtext( side=1, "False Discovery Rate (FP/FP+TP)", line = 2 )
  dev.off()

  # Plot the measurements over depth
  # TPR over Depth
  # FDR over Depth
  i_max_coverage = ceiling( max( vi_depths, na.rm=TRUE ) * 1.1 )
  vi_serial_ticks= seq( 0, i_max_coverage, ceiling( max( vi_depths,na.rm=TRUE ) / 10.0 ) ) 
  vstr_serial_ticks=paste( vi_serial_ticks )
  for( i_name_single_plots in 1:length( vstr_line_names ) )
  {
    str_file_path=dirname( str_pdf_file_name )
    str_file_base=basename( str_pdf_file_name )
    pdf( file.path( str_file_path, paste( "FDR", str_file_base,sep="_")), useDingbats=FALSE)
    # Plot comparison line if given (FDR)
    str_cur_depth=vstr_line_names[ i_name_single_plots ]
    # Plot predictor line (FDR)
    plot(x=vi_depths, y=list_FDR[[ vstr_line_names[ i_name_single_plots ]]], col=vstr_exome_colors[ i_name_single_plots ], ylim=c(0,1), xlab=NA, ylab=NA, xaxt="n" )
    lines(x=vi_depths, y=list_FDR[[ vstr_line_names[ i_name_single_plots ]]], col=vstr_exome_colors[ i_name_single_plots ], ylim=c(0,1) )
    points(x=vi_depths, y=list_FDR[[ vstr_line_names[ i_name_single_plots ]]], pch=24, col=vstr_exome_colors[ i_name_single_plots ], bg=vstr_roc_colors )
    if( ! is.null( vi_depths_compare ) )
    {
      lines(x=vi_depths_compare, y=list_FDR_compare[[ vstr_line_names[ i_name_single_plots ]]], col=vstr_exome_colors[ i_name_single_plots ], lty=3 )
      points(x=vi_depths_compare, y=list_FDR_compare[[ vstr_line_names[ i_name_single_plots ]]], pch=24, col=vstr_exome_colors[ i_name_single_plots ], bg=vstr_roc_colors )
    }
    title( "FDR vs Min RNA-Seq Coverage" )
    # Legend
    legend( "topright", legend=vstr_legend_labels_rna, fill=vstr_legend_fill_rna, border=vstr_legend_border_rna, title="Min RNA-Seq Cov", cex=.75 )
    legend( "bottomleft", legend=vstr_legend_labels_exome, col=vstr_legend_col_exome, title=str_legend_vary_title, lty=vstr_legend_line_exome, lwd=2, cex=.75, pch=vstr_legend_shape_exome )
    axis( 2, at=vi_serial_ticks, labels=vstr_serial_ticks )
    mtext( side=1, "Min RNA-Seq Coverage", line = 2 )
    axis( 1, at=vi_serial_ticks, labels=vstr_serial_ticks ) 
    mtext( side=2, "False Discovery Rate (FP/FP+TP)", line = 2 )
    dev.off()

    pdf( file.path( str_file_path, paste( "TPR", str_file_base,sep="_")), useDingbats=FALSE)
    # Plot comparison line if given (TPR)
    str_cur_depth=vstr_line_names[ i_name_single_plots ]
    # Plot predictor line (TPR)
    plot(x=vi_depths, y=list_TPR[[ vstr_line_names[ i_name_single_plots ]]], col=vstr_exome_colors[ i_name_single_plots ], ylim=c(0,1), xlab=NA, ylab=NA, xaxt="n" )
    lines(x=vi_depths, y=list_TPR[[ vstr_line_names[ i_name_single_plots ]]], col=vstr_exome_colors[ i_name_single_plots ] )
    points(x=vi_depths, y=list_TPR[[ vstr_line_names[ i_name_single_plots ]]], pch=24, col=vstr_exome_colors[ i_name_single_plots ], bg=vstr_roc_colors )
    if( ! is.null( vi_depths_compare ) )
    {
      lines(x=vi_depths_compare, y=list_TPR_compare[[ vstr_line_names[ i_name_single_plots ]]], col=vstr_exome_colors[ i_name_single_plots ], lty=3 )
      points(x=vi_depths_compare, y=list_TPR_compare[[ vstr_line_names[ i_name_single_plots ]]], pch=24, col=vstr_exome_colors[ i_name_single_plots ], bg=vstr_roc_colors )
    }
    title( "TPR vs Min RNA-Seq Coverage" )
    # Legend
    legend( "bottomleft", legend=vstr_legend_labels_rna, fill=vstr_legend_fill_rna, border=vstr_legend_border_rna, title="Vary RNA Cov", cex=.75 )
    legend( "topright", legend=vstr_legend_labels_exome, col=vstr_legend_col_exome, title=str_legend_vary_title, lty=vstr_legend_line_exome, lwd=2, cex=.75, pch=vstr_legend_shape_exome )
    vstr_roc_ticks=paste( VI_ROC_TICKS )
    axis( 2, at=vi_serial_ticks, labels=vstr_serial_ticks )
    mtext( side=1, "Min RNA-Seq Coverage", line = 2 )
    axis( 1, at=vi_serial_ticks, labels=vstr_serial_ticks ) 
    mtext( side=2, "True Positive Rate (TP/TP+FN)", line = 2 )
    dev.off()
  }
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

# Return features that did not call the same genotype
# df_data: data frame with 8 entries
# vi_indices, indices that are looked at
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

func_vary_coverage_and_measure_classes = function( df_data, vi_vary_truth, vi_vary_prediction, str_file_base_name, f_vary_truth ) 
{
  # Go through Exome min coverage settings
  # Hold the TPR and FDR lines for each setting in a list
  list_TPR = list()
  list_FDR = list()
  for( i_cur_index in 1:length(vi_vary_truth) )
  {
    i_cur_truth_coverage = vi_vary_truth[ i_cur_index ]
    i_cur_pred_coverage = vi_vary_prediction[ i_cur_index ]

    # Reduce the full data frame by requiring the DNA side down to be a min coverage
    # Set the initial problem state.
    # Min Exome coverage
    df_tab = df_data
    vi_remove_calls = which( df_tab[[ C_I_PRIMARY_DEPTH ]] < i_cur_truth_coverage )
    vi_remove_calls = unique( union( vi_remove_calls, which( df_tab[[ C_I_SECONDARY_DEPTH ]] < i_cur_pred_coverage ) ) )
    if( length( vi_remove_calls ) > 0 )
    {
      df_tab = df_tab[-1*vi_remove_calls,]
    }

    # Different measurements used later on
    ### Primary calls ( not Na and atleast a minimim coverage given above )
    vi_primary_calls = which( !is.na(  df_tab[[ C_I_PRIMARY_GT ]] ) )
    ## Secondary
    ### Calls with or without depth
    vi_secondary_calls = which( !is.na( df_tab[[ C_I_SECONDARY_GT ]] ) )
    # Update the secondary calls, if the genotypes do no match the primary calls
    # Remove from the secondary calls
    vi_common_calls = intersect( vi_secondary_calls, vi_primary_calls )
    # NOTE:
    # Calls are removed from the secondary calls that do not have the same variant call
    # This is because downstream the secondary calls are assumed to be not only called but correctly (matching the DNA).
    # Removing the calls from scondary calls shuffles the false calls into:
    # if called in the DNA: False Negatives
    # if not called in the DNA: True Negatives, which are not used in calculations so far.
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
    
    # Do this for each depth of interest (given in constants).
    vi_roc_depths = sort( VI_ROC_DEPTHS, decreasing=FALSE )
    vNA_init = rep(NA,length( vi_roc_depths ))
    df_roc = data.frame( TPR=vNA_init, FDR=vNA_init, TP=vNA_init, FP=vNA_init, FN=vNA_init )
    rownames( df_roc ) = vi_roc_depths
    for( i_cur_roc_depth in vi_roc_depths )
    {
      if( i_cur_roc_depth < i_cur_pred_coverage )
      {
        # Store infor to write.
        vi_TPR = c( vi_TPR, NA )
        vi_FDR = c( vi_FDR, NA )
        vi_TP = c( vi_TP, NA )
        vi_FP = c( vi_FP, NA )
        vi_FN = c( vi_FN, NA )
      } else {
        # NOTE:
        # Secondary calls that do not meet the depth criteria change
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
    write.table( df_roc, file = file.path( str_output_dir, paste( str_file_base_name, paste("roc_truth_",i_cur_truth_coverage,"_pred_",i_cur_pred_coverage,".txt",sep=""), sep = "_" ) ) )
  }
  return( list( TPR=list_TPR, FDR=list_FDR, DEPTHS=vi_roc_depths ) )
}

# Read in tab data
# Remove features that have NA in both DNA and RNA
# Remove features that have no coverage in DNA and RNA
# Fix the data modes just incase
func_read_and_filter = function( str_input_file )
{
  print( "Reading file:" )
  print( str_input_file )

  # id \t chr_position \t Ref_bas \t call_genotype \t depth \t chr_position \t ref_bas \t call_genotype \t depth
  # 1-4 is for the reference sample
  # 5-8 is for the secondary assay sample
  df_orig = read.table( str_input_file )
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

  # Set type
  df_orig[[ C_I_PRIMARY_POS ]] = as.character( df_orig[[ C_I_PRIMARY_POS ]] )
  df_orig[[ C_I_PRIMARY_REF ]] = as.character( df_orig[[ C_I_PRIMARY_REF ]] )
  df_orig[[ C_I_PRIMARY_GT ]] = as.character( df_orig[[ C_I_PRIMARY_GT ]] )
  df_orig[[ C_I_PRIMARY_DEPTH ]] = as.numeric( df_orig[[ C_I_PRIMARY_DEPTH ]] )
  df_orig[[ C_I_SECONDARY_POS ]] = as.character( df_orig[[ C_I_SECONDARY_POS ]] )
  df_orig[[ C_I_SECONDARY_REF ]] = as.character( df_orig[[ C_I_SECONDARY_REF ]] )
  df_orig[[ C_I_SECONDARY_GT ]] = as.character( df_orig[[ C_I_SECONDARY_GT ]] )
  df_orig[[ C_I_SECONDARY_DEPTH ]] = as.numeric( df_orig[[ C_I_SECONDARY_DEPTH ]] )
  return( df_orig )
}

# Handle arguments
# Require the output directory
v_str_files = lsArgs$args
str_title_key = lsArgs$options$str_title_key
str_output_dir = lsArgs$options$str_output_dir
if( is.null(str_output_dir) )
{
  stop("Please include the output directory. Use -o.")
}
dir.create( str_output_dir )

# Information about the files as a group
# Holds file information for the global images at the end.
print( "STARTING visualize_mutation_depth_tab_files.R" )
print( "Number of files reading:" )
print( length( v_str_files ) )

ls_classes_compare_vary_min_truth = list( "TPR"=NULL, "FDR"=NULL, "DEPTHS"=NULL )
ls_classes_compare_vary_min_pred = list( "TPR"=NULL, "FDR"=NULL, "DEPTHS"=NULL )
if( ! is.null( lsArgs$options$str_compare_file ) )
{
  # Read in a reference file to compare against if given
  df_compare = func_read_and_filter( str_input_file=lsArgs$options$str_compare_file )

  # Vary the problem space holding the pred min coverage at 1 and then varying the truth sample space
  # Then look at min coverage for the pred in the resulting sample space
  # Add a base to the txt file that can have the sample and title keys for uniqueness
  str_cur_base_name = gsub( "\\.", "_", basename( lsArgs$options$str_compare_file ) )
  if( !is.null( lsArgs$options$str_sample_key ))
  {
    str_cur_base_name = paste( str_cur_base_name, lsArgs$options$str_sample_key, sep = "_" )
  }
  if( !is.null( str_title_key ))
  {
    str_cur_base_name = paste( str_cur_base_name, str_title_key, sep = "_" )
  }
  ls_classes_compare_vary_min_truth = func_vary_coverage_and_measure_classes( df_data=df_compare, vi_vary_truth=VI_ROC_TRUTH_MIN_DEPTH,
                                          vi_vary_prediction=rep(I_SELECTED_PRED_MIN_COV,length(VI_ROC_TRUTH_MIN_DEPTH)),
                                          str_file_base_name=str_cur_base_name )

  # Vary the problem space holding truth at 10 min coverage and vary pred min coverage
  # then investigate each pred min coverage in the resulting sample space
  ls_classes_compare_vary_min_pred = func_vary_coverage_and_measure_classes( df_data=df_compare, vi_vary_truth=rep(I_SELECTED_TRUTH_MIN_COV,length(VI_ROC_PRED_MIN_DEPTH)),
                                          vi_vary_prediction=VI_ROC_PRED_MIN_DEPTH,
                                          str_file_base_name=str_cur_base_name )
}

# Process each file
for( str_file in v_str_files )
{
  # Read in and filter perform initial filtering
  df_orig = func_read_and_filter( str_input_file=str_file )
  # Mean after minimal filtering
  i_mean_depth=round( mean( df_orig[[ C_I_SECONDARY_DEPTH ]] ),2 )

  # Add a base to the txt file that can have the sample and title keys for uniqueness
  str_cur_base_name = gsub( "\\.", "_", basename( str_file ) )
  if( !is.null( lsArgs$options$str_sample_key ))
  {
    str_cur_base_name = paste( str_cur_base_name, lsArgs$options$str_sample_key, sep = "_" )
  }
  if( !is.null( str_title_key ))
  {
    str_cur_base_name = paste( str_cur_base_name, str_title_key, sep = "_" )
  }

  # Vary the problem space holding the pred min coverage at 1 and then varying the truth sample space
  # Then look at min coverage for the pred in the resulting sample space
  ls_classes_vary_min_truth = func_vary_coverage_and_measure_classes( df_data=df_orig, vi_vary_truth=VI_ROC_TRUTH_MIN_DEPTH,
                                          vi_vary_prediction=rep(I_SELECTED_PRED_MIN_COV,length(VI_ROC_TRUTH_MIN_DEPTH)),
                                          str_file_base_name=str_cur_base_name,
                                          f_vary_truth=TRUE ) 

  # Go through indices for each ROC line, filter at depth in RNASEQ
  # Create file name with optional sample tag if given
  str_file_roc_rnaseq = file.path( str_output_dir, paste( str_title_key, "_roc_truth_vary_pred_", I_SELECTED_PRED_MIN_COV, ".pdf", sep="") )
  if( !is.null( lsArgs$options$str_sample_key ) )
  {
      str_file_roc_rnaseq = file.path( str_output_dir, paste( lsArgs$options$str_sample_key, "_", str_title_key, "_roc_truth_vary_pred_", I_SELECTED_PRED_MIN_COV, ".pdf", sep="") )
  }
  func_plot_roc( list_TPR=ls_classes_vary_min_truth[[ "TPR" ]], list_FDR=ls_classes_vary_min_truth[[ "FDR" ]], vi_depths=ls_classes_vary_min_truth[[ "DEPTHS" ]],
                 list_TPR_compare=ls_classes_compare_vary_min_truth[["TPR"]], list_FDR_compare=ls_classes_compare_vary_min_truth[["FDR"]],
                 vi_depths_compare=ls_classes_compare_vary_min_truth[["DEPTHS"]],
                 i_mean_depth=i_mean_depth, str_method_name=lsArgs$options$str_method_name, str_method_name_compare=lsArgs$options$str_method_name_compare,
                 str_pdf_file_name=str_file_roc_rnaseq,
                 str_title=paste("TPR vs FDR varying Min Cov ( RNA-Seq Min Cov ",I_SELECTED_PRED_MIN_COV," )",sep=""),
                 str_legend_vary_title="Min Exome Cov" )

  # Make serial plots if requested
  if( lsArgs$options$f_make_serial_plots )
  {
    # Plot for each setting 
    for( str_name in names( ls_classes_vary_min_truth[[ "TPR" ]] ) )
    {
      lvf_serial_TPR=list()
      lvf_serial_TPR_compare=list()
      lvf_serial_TPR[[ str_name ]] = ls_classes_vary_min_truth[["TPR"]][[ str_name ]]
      lvf_serial_TPR_compare[[ str_name ]] = ls_classes_compare_vary_min_truth[["TPR"]][[ str_name ]]
      lvf_serial_FDR=list()
      lvf_serial_FDR_compare=list()
      lvf_serial_FDR[[ str_name ]] = ls_classes_vary_min_truth[["FDR"]][[ str_name ]]
      lvf_serial_FDR_compare[[ str_name ]] = ls_classes_compare_vary_min_truth[["FDR"]][[ str_name ]]
      # Create file name with optional sample tag if given
      str_file_roc_rnaseq = file.path( str_output_dir, paste( str_title_key,  "_roc_truth_", str_name,"_pred_",I_SELECTED_PRED_MIN_COV,".pdf",sep="") )
      if( !is.null( lsArgs$options$str_sample_key ) )
      {
          str_file_roc_rnaseq = file.path( str_output_dir, paste( lsArgs$options$str_sample_key, "_", str_title_key, "_roc_truth_", str_name,"_pred_",I_SELECTED_PRED_MIN_COV,".pdf", sep="") )
      }
      func_plot_roc( list_TPR=lvf_serial_TPR, list_FDR=lvf_serial_FDR, vi_depths=ls_classes_vary_min_truth[[ "DEPTHS" ]],
                 list_TPR_compare=lvf_serial_TPR_compare, list_FDR_compare=lvf_serial_FDR_compare, vi_depths_compare=ls_classes_compare_vary_min_truth[["DEPTHS"]],
                 i_mean_depth=i_mean_depth, str_method_name=lsArgs$options$str_method_name, str_method_name_compare=lsArgs$options$str_method_name_compare,
                 str_pdf_file_name=str_file_roc_rnaseq,
                 str_title=paste("TPR vs FDR ( DNA-Seq Min ",str_name,"; RNA-Seq Min Cov ",I_SELECTED_PRED_MIN_COV," )",sep=""),
                 str_legend_vary_title="Min Exome Cov" )
    }
  }

  # Vary the problem space holding truth at 10 min coverage and vary pred min coverage
  # then investigate each pred min coverage in the resulting sample space
  ls_classes_vary_min_pred = func_vary_coverage_and_measure_classes( df_data=df_orig, vi_vary_truth=rep(I_SELECTED_TRUTH_MIN_COV,length(VI_ROC_PRED_MIN_DEPTH)),
                                          vi_vary_prediction=VI_ROC_PRED_MIN_DEPTH,
                                          str_file_base_name=str_cur_base_name,
                                          f_vary_truth=FALSE )

  # Go through indices for each ROC line, filter at depth in RNASEQ
  # Create file name with optional sample tag if given
  str_file_roc_rnaseq = file.path( str_output_dir, paste( str_title_key,"_roc_truth_",I_SELECTED_TRUTH_MIN_COV,"_pred_vary.pdf",sep="") )
  if( !is.null( lsArgs$options$str_sample_key ) )
  {
      str_file_roc_rnaseq = file.path( str_output_dir, paste( lsArgs$options$str_sample_key, "_", str_title_key, "_roc_truth_",I_SELECTED_TRUTH_MIN_COV,"_pred_vary.pdf", sep="") )
  }
  func_plot_roc( list_TPR=ls_classes_vary_min_pred[[ "TPR" ]], list_FDR=ls_classes_vary_min_pred[[ "FDR" ]], vi_depths=ls_classes_vary_min_pred[[ "DEPTHS" ]],
                 list_TPR_compare=ls_classes_compare_vary_min_pred[["TPR"]], list_FDR_compare=ls_classes_compare_vary_min_pred[["FDR"]],
                 vi_depths_compare=ls_classes_compare_vary_min_pred[["DEPTHS"]],
                 i_mean_depth=i_mean_depth, str_method_name=lsArgs$options$str_method_name, str_method_name_compare=lsArgs$options$str_method_name_compare,
                 str_pdf_file_name=str_file_roc_rnaseq,
                 str_title=paste("TP vs FDR varying Min Cov ( DNA-Seq Min Cov ",I_SELECTED_TRUTH_MIN_COV," )",sep=""),
                 str_legend_vary_title="Min RNA-Seq Cov" )

  # Make serial plots if requested
  if( lsArgs$options$f_make_serial_plots )
  {
    # Plot for each setting (vary truth)
    for( str_name in names( ls_classes_vary_min_truth[[ "TPR" ]] ) )
    {
      lvf_serial_TPR=list()
      lvf_serial_TPR_compare=list()
      lvf_serial_TPR[[ str_name ]] = ls_classes_vary_min_truth[["TPR"]][[ str_name ]]
      lvf_serial_TPR_compare[[ str_name ]] = ls_classes_compare_vary_min_truth[["TPR"]][[ str_name ]]
      lvf_serial_FDR=list()
      lvf_serial_FDR_compare=list()
      lvf_serial_FDR[[ str_name ]] = ls_classes_vary_min_truth[["FDR"]][[ str_name ]]
      lvf_serial_FDR_compare[[ str_name ]] = ls_classes_compare_vary_min_truth[["FDR"]][[ str_name ]]
      # Create file name with optional sample tag if given
      str_file_roc_rnaseq = file.path( str_output_dir, paste( basename( str_file ), "_", str_title_key, "_roc_truth_", str_name,"_pred_",I_SELECTED_PRED_MIN_COV,".pdf",sep="") )
      if( !is.null( lsArgs$options$str_sample_key ) )
      {
          str_file_roc_rnaseq = file.path( str_output_dir, paste( lsArgs$options$str_sample_key, "_", basename( str_file ), "_", str_title_key, "_roc_truth_", str_name,"_pred_",I_SELECTED_PRED_MIN_COV,".pdf", sep="") )
      }
      func_plot_roc( list_TPR=lvf_serial_TPR, list_FDR=lvf_serial_FDR, vi_depths=ls_classes_vary_min_truth[[ "DEPTHS" ]],
                 list_TPR_compare=lvf_serial_TPR_compare, list_FDR_compare=lvf_serial_FDR_compare, vi_depths_compare=ls_classes_compare_vary_min_truth[["DEPTHS"]],
                 i_mean_depth=i_mean_depth, str_method_name=lsArgs$options$str_method_name, str_method_name_compare=lsArgs$options$str_method_name_compare,
                 str_pdf_file_name=str_file_roc_rnaseq,
                 str_title=paste("TPR vs FDR ( DNA-Seq Min ",str_name,"; RNA-Seq Min Cov ",I_SELECTED_PRED_MIN_COV," )",sep=""),
                 str_legend_vary_title="Min Exome Cov" )


      # Plot for each setting (vary pred)
      for( str_name_pred in names( ls_classes_vary_min_pred[[ "TPR" ]] ))
      {
          lvf_serial_TPR=list()
          lvf_serial_TPR_compare=list()
          lvf_serial_TPR[[ str_name_pred ]]=ls_classes_vary_min_pred[["TPR"]][[ str_name_pred ]]
          lvf_serial_TPR_compare[[ str_name_pred ]]=ls_classes_compare_vary_min_pred[["TPR"]][[ str_name_pred ]]
          lvf_serial_FDR=list()
          lvf_serial_FDR_compare=list()
          lvf_serial_FDR[[ str_name_pred ]]=ls_classes_vary_min_pred[["FDR"]][[ str_name_pred ]]
          lvf_serial_FDR_compare[[ str_name_pred ]]=ls_classes_compare_vary_min_pred[["FDR"]][[ str_name_pred ]]
          # Create file name with optional sample tag if given
          str_file_roc_rnaseq = file.path( str_output_dir, paste( basename( str_file ), "_", str_title_key, "_roc_truth_",I_SELECTED_TRUTH_MIN_COV,"_pred_",str_name_pred,".pdf",sep="") )
          if( !is.null( lsArgs$options$str_sample_key ) )
          {
              str_file_roc_rnaseq = file.path( str_output_dir, paste( lsArgs$options$str_sample_key, "_", basename( str_file ), "_", str_title_key, "_roc_truth_",I_SELECTED_TRUTH_MIN_COV,"_pred_",str_name_pred,".pdf", sep="") )
          }
          func_plot_roc( list_TPR=lvf_serial_TPR, list_FDR=lvf_serial_FDR, vi_depths=ls_classes_vary_min_pred[[ "DEPTHS" ]],
                 list_TPR_compare=lvf_serial_TPR_compare, list_FDR_compare=lvf_serial_FDR_compare,
                 vi_depths_compare=ls_classes_compare_vary_min_pred[["DEPTHS"]],
                 i_mean_depth=i_mean_depth, str_method_name=lsArgs$options$str_method_name, str_method_name_compare=lsArgs$options$str_method_name_compare,
                 str_pdf_file_name=str_file_roc_rnaseq,
                 str_title=paste("TP vs FDR ( DNA-Seq Min ",I_SELECTED_TRUTH_MIN_COV,"; RNA-Seq Min Cov ",str_name_pred,")",sep=""),
                 str_legend_vary_title="Min RNA-Seq Cov" )
      }
    }
  }
}
