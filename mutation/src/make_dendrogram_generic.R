#!/usr/bin/env Rscript

library( optparse )

# Command line arguments
pArgs <- OptionParser( usage = "%prog [options] input.matrix output.pdf" )
pArgs <- add_option( pArgs, c("-i","--input_matrix"), type="character", action="store", dest="str_input_matrix", metavar="Input_matrix", help="Input genotype matrix to visualize." )
pArgs <- add_option( pArgs, c("-p","--output_png"), type="character", action="store", dest="str_output_png", metavar="Output_png", help="An output png of the hclusted samples." )
pArgs <- add_option( pArgs, c("-o","--output_distance_matrix"), type="character", action="store", dest="str_distance_matrix", metavar="Output_distance", help="An output file of distance matrices." )
pArgs <- add_option( pArgs, c("-f","--distance_function"), type="character", action="store", dest="str_distance_function", metavar="Input_function", help="A function used for measuring distance." )
args <- parse_args( pArgs )
print(args)
# Source in the distance function
source( args$str_distance_function )

# Read in input matrix
data = read.table( args$str_input_matrix )

# Make distance matrix
d = dist_matrix(data) #!!!

# Write distance matrix to directory
write.table(d, args$str_distance_matrix, quote=F, sep="\t")

# Hclust and write to png
h = hclust(as.dist(d))
png( args$str_output_png, width=20)
plot(h, cex=0.3)
dev.off()
