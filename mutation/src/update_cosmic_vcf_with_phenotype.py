#!/usr/bin/env python

# Libraries
import argparse
import csv
import gzip
import os

# Constants
## COSMIC VCF
c_VCF_COMMENT = "#"
c_VCF_DELIM = "\t"
c_VCF_COSMIC_ID_INDEX = 2
c_VCF_ANN_DELIM = ";"
c_VCF_HEADER = "#CHROM"
c_VCF_ANN_INDEX = 7

## COSMIC
c_PHENOTYPE_DELIM = "\t"
c_COSMIC_NO_ENTRY = [ "NS" ]
i_COSMIC_ID_INDEX = 12
str_COSMIC_ID = "COSMIC_ID"
str_COSMIC_ID_DES = "COSMIC mutation id (unique)."
#i_COSMIC_GENE_NAME_INDEX = 0
#str_COSMIC_GENE_NAME = "GENE_NAME"
#str_COSMIC_GENE_NAME_DES = "The gene name for which the data has been curated in COSMIC. In most cases this is the accepted HGNC identifier."
i_COSMIC_TISSUE_INDEX = 7
i_COSMIC_SUBTISSUE_INDEX = 8
str_COSMIC_TISSUE = "TISSUE"
str_COSMIC_TISSUE_DES = "The primary tissue/cancer and subtype from which the sample originated."
i_COSMIC_TUMOR_INDEX = 9
i_COSMIC_SUBTUMOR_INDEX = 10
str_COSMIC_TUMOR = "TUMOR"
str_COSMIC_TUMOR_DES = "The histological classification of the sample."
i_COSMIC_FATHMM_INDEX = 21
str_COSMIC_FATHMM = "FATHMM"
str_COSMIC_FATHMM_DES = "FATHMM (Functional Analysis through Hidden Markov Models). 'Pathogenic'=Cancer or damaging, 'Neutral'=Passanger or Tolerated."
i_COSMIC_SOMATIC_INDEX = 22
str_COSMIC_SOMATIC = "SOMATIC"
str_COSMIC_SOMATIC_DES = "Information on whether the sample was reported to be Confirmed Somatic. 'Confirmed somatic'=if the mutation has been confimed to be somatic in the experiment by sequencing both the tumour and a matched normal from the same patient, 'Previously Observed'=when the mutation has been reported as somatic previously but not in current paper, 'variant of unknown origin'=when the mutation is known to be somatic but the tumour was sequenced without a matched normal"
i_COSMIC_PUBMED_INDEX = 23
str_COSMIC_PUBMED = "PUBMED_COSMIC"
str_COSMIC_PUBMED_DES = "The PUBMED ID for the paper that the sample was noted in COSMIC."

# Arg parse
prsr_arguments = argparse.ArgumentParser( prog = "update_cosmic_vcf_with_phenotype.py", description = "Update a cosmic vcf with phenotype information.", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "--in_vcf", required=True, dest="str_vcf_file_in", action="store", help="Input vcf file." )
prsr_arguments.add_argument( "--in_pheno", required=True, dest="str_pheno_file_in", action="store", help="Input Phenotype file from cosmic." )
prsr_arguments.add_argument( dest="str_output_vcf", action="store", help="Output VCF file updated with phenotype annotations." )
args = prsr_arguments.parse_args()


# Indicate files
print "VCF File: " + args.str_vcf_file_in
print "Phenotype File: " + args.str_pheno_file_in
print "Writing to: " + args.str_output_vcf


# Holds the phenotype information
## { cosmic_id: c_VCF_ANN_DELIM.join( [ tissue, tumor, FATHMM_prediction, somatic_status, pubmed_id ) }
dict_pheno = {}

# Skips the header
f_skip_header = True

# Read in phenotype information
print "Phenotype file..."
with gzip.open( args.str_pheno_file_in, "rb" ) as hndl_phenotype_in:
  for lstr_line in csv.reader( hndl_phenotype_in, delimiter=c_PHENOTYPE_DELIM ):
    # First line is a header
    if f_skip_header:
      f_skip_header = False
      continue

    # Tissue
    str_tissue = lstr_line[ i_COSMIC_TISSUE_INDEX ]
    if not lstr_line[ i_COSMIC_SUBTISSUE_INDEX ] in c_COSMIC_NO_ENTRY:
      str_tissue = str_tissue + " -- " + lstr_line[ i_COSMIC_SUBTISSUE_INDEX ]

    # Tumor
    str_tumor = lstr_line[ i_COSMIC_TUMOR_INDEX ]
    if not lstr_line[ i_COSMIC_SUBTUMOR_INDEX ] in c_COSMIC_NO_ENTRY:
      str_tumor = str_tumor + " -- " + lstr_line[ i_COSMIC_SUBTUMOR_INDEX ]

    # Store data of interest
    dict_pheno[ lstr_line[ i_COSMIC_ID_INDEX ] ] = c_VCF_ANN_DELIM.join( [ "=".join( [ str_COSMIC_ID, lstr_line[ i_COSMIC_ID_INDEX ] ] ), 
                                                     "=".join( [ str_COSMIC_TISSUE, str_tissue ] ),
                                                     "=".join( [ str_COSMIC_TUMOR, str_tumor ] ),
                                                     "=".join( [ str_COSMIC_FATHMM, lstr_line[ i_COSMIC_FATHMM_INDEX ] ] ),
                                                     "=".join( [ str_COSMIC_SOMATIC, lstr_line[ i_COSMIC_SOMATIC_INDEX ] ] ),
                                                     "=".join( [ str_COSMIC_PUBMED, lstr_line[ i_COSMIC_PUBMED_INDEX ] ] ) ] )
print "Read in " + str( len( dict_pheno ) )  + " entries."

# Open handle to vcf
print "Reading VCF file..."
with gzip.open( args.str_vcf_file_in, "rb" ) as hndl_vcf_in:
  with gzip.open( args.str_output_vcf, "wb" ) as hndl_vcf_out:

    # Read in VCF file
    for str_vcf_line in hndl_vcf_in:
      # Split line, CSV reader is no good, some fields are too long for it.
      lstr_vcf_line = str_vcf_line.split( c_VCF_DELIM )

      # Store comments
      if lstr_vcf_line[0][0] == c_VCF_COMMENT:
        if lstr_vcf_line[ 0 ] == c_VCF_HEADER:
          # Write new comments
          for str_id, str_des in zip( [ str_COSMIC_ID, str_COSMIC_TISSUE,
                                        str_COSMIC_TUMOR, str_COSMIC_FATHMM, str_COSMIC_SOMATIC, str_COSMIC_PUBMED ],
                                      [ str_COSMIC_ID_DES, str_COSMIC_TISSUE_DES,
                                        str_COSMIC_TUMOR_DES, str_COSMIC_FATHMM_DES, str_COSMIC_SOMATIC_DES, str_COSMIC_PUBMED_DES ] ):
            hndl_vcf_out.write( "##INFO=<ID=" + str_id + ",Type=String,Description=\"" + str_des + "\">\n" )
        # Write line
        hndl_vcf_out.write( c_VCF_DELIM.join( lstr_vcf_line ) )
        continue

      # Update line if the cosmic id matches
      if lstr_vcf_line[ c_VCF_COSMIC_ID_INDEX ] in dict_pheno:
        lstr_vcf_line[ c_VCF_ANN_INDEX ] = c_VCF_ANN_DELIM.join( [ dict_pheno[ lstr_vcf_line[ c_VCF_COSMIC_ID_INDEX ] ], lstr_vcf_line[ c_VCF_ANN_INDEX ] ] )

      # Write VCF file
      hndl_vcf_out.write( c_VCF_DELIM.join( lstr_vcf_line ) )
