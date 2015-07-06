#!/usr/bin/env python

import argparse
import json
import os
import requests
import time
import urllib

# Constants
str_response_job_id = "jobid"
str_response_log = "logurl"
str_response_status = "status"
str_response_url = "resultfileurl"
str_response_pass = u'Success'
lstr_response_fail = [ u'Error', u'submissionfailed' ]

# VCF
str_vcf_comment = "#"
str_cravat_header = "\t".join( [ "# UID", "Chr.", "Position", "Strand", "Ref. base", "Alt. base" ] )
c_VCF_delimiter = "\t"
i_vcf_chr = 0
i_vcf_position = 1
i_vcf_ref = 3
i_vcf_alt = 4

# Parse arguments
prsr_arguments = argparse.ArgumentParser( prog = "annotate_with_cravat.py", description = "Annotate VCF file with cravat.", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "--classifier", dest = "str_classifier", action = "store", required=True, help = "Tissue type." )
prsr_arguments.add_argument( "--is_hg18", dest = "f_hg_18", action = "store", help = "Indicates if Hg18 (True=Hg18; False=Hg19)." )
prsr_arguments.add_argument( "--email", dest = "str_email", action = "store", required=True, help = "Email contact for job." )
prsr_arguments.add_argument( "--max_attempts", dest = "i_max_attempts", default=100, action = "store", type=int, help = "Max attempts of querying response before timing out." )
prsr_arguments.add_argument( "--wait", dest = "i_wait", action = "store", default=10, type=int, help = "Wait in seconds before querying the response." )
prsr_arguments.add_argument( dest = "str_input_file", action = "store", help = "VCF file to update." )
prsr_arguments.add_argument( dest = "str_output_dir", action = "store", help = "Output Zip file (please use the extension .zip) Can be opened with 'unzip'." )
args_call = prsr_arguments.parse_args()

def func_check_job( str_json_id ):
  """
  Uses the CRAVAT status web service to check if a job is complete.
  """

  # Send file over to service
  response_check_job = requests.get( "http://www.cravat.us/rest/service/status", params={ str_response_job_id: str_json_id } )
  return response_check_job.json()

def func_download_cravat_result( str_cravat_url, str_download_location ):
  """
  Using the download URL from CRAVAT, downloads to the given location.
  """

  # Download Cravat result
  try:
    urllib.urlretrieve( str_cravat_url, str_download_location )
    return True
  except IOError:
    print "annotate_with_cravat::Error in retrieving results. " + str( IOError )
    return False

def func_get_cravat_response( str_json_id, i_max_attempts, i_wait ):
  """
  Checks CRAVAT for a period of time until the success, failure, or time out are given for a job id.
  """

  for x_attempt in xrange( i_max_attempts ):
    print "annotate_with_cravat::Attempting to get CRAVAT response. Attempt " + str( x_attempt )
    time.sleep( i_wait )
    response_json = func_check_job( str_json_id )
    if not response_json:
      print "annotate_with_cravat::Error did not get a resposne from CRAVAT"
      return( None )
    str_success = str( response_json.get( str_response_status, None) )
    if ( not str_success ) or ( str_success in lstr_response_fail ):
      print " ".join( [ "annotate_with_cravat::Error was not successful in getting info from CRAVAT.",
                        "Job id =" + str( str_json_id ) + ".",
                        "Success =" + str( str_success ) + ".",
                        response_json.get( str_response_log, "" ) ] )
      return( None )
    if str_success == str_response_pass:
      return response_json.get( str_response_url, None )
  print " ".join( [ "annotate_with_cravat::Error gave up on CRAVAT after " + str( i_max_attempts * i_wait ) + " seconds.",
                    "Job id =" + str( str_json_id ) + ".",
                    "Success =" + str( str_success ) + "." ] )
  return( None )

def func_vcf_to_cravat_mutations( str_path ):
  """
  Read in a VCF file and reduce it in size for CRAVAT
  """

  # Read in file into a string
  lstr_reduced_vcf = [ str_cravat_header ]
  str_sample = os.path.basename( str_path )
  i_counter = 0
  with open( str_path, "rb" ) as hndl_in:
    for str_line in hndl_in:
      if str_line[ 0 ] == str_vcf_comment:
        continue
    lstr_cur_line = str_line.split( c_VCF_delimiter )
    lstr_reduced_vcf.append( "\t".join( [ str( i_counter ), lstr_cur_line[ i_vcf_chr ], lstr_cur_line[ i_vcf_position ],
                                          "+", lstr_cur_line[ i_vcf_ref ], lstr_cur_line[ i_vcf_alt ] ] ) )
    i_counter = i_counter + 1

  return "\n".join( lstr_reduced_vcf )

def func_request_cravat_service( str_vcf_path, str_classifier, f_hg_18, str_email ):
  """
  Request a job to occur with CRAVAT.
  """

#  # Read in file into a string
#  str_initial_vcf = func_vcf_to_cravat_mutations( str_vcf_path )
#  if not str_initial_vcf:
#    print " ".join( [ "annotate_with_cravat::Error did not read VCF file.",
#                      "Path =" + str( str_vcf_path ) + "." ] )
#    return( None )

  # Encode request info
  pyld_request = { "chasmclassifier": str( str_classifier ),
                   "mupitinput": "on",
                   "hg18": "on" if f_hg_18 else "off",
                   "tsvreport": "on",
                   "analyses": "CHASM;SnvGet",
                   "functionalannotation": "on",
                   "analysistype": "driver",
                   "email": str_email }

  # Send file over to service
  response_cravat = requests.post( "http://www.cravat.us/rest/service/submit", files={"inputfile": open( str_vcf_path )}, data=pyld_request )
  json_response = response_cravat.json()
  # Get return (job id)
  return json_response.get( str_response_job_id, None )

# Runs
# Request CRAVAT service
str_job_id = func_request_cravat_service( args_call.str_input_file,
                                          args_call.str_classifier,
                                          args_call.f_hg_18,
                                          args_call.str_email )
if not str_job_id:
  print " ".join( [ "annotate_with_cravat::Error Job id not found.",
                    "Job id =" + str( str_job_id ) ] )
  exit( 100 )

# Wait for result
str_download_url = func_get_cravat_response( str_job_id, args_call.i_max_attempts, args_call.i_wait )
if not str_download_url:
  print " ".join( [ "annotate_with_cravat::Error did not recieve valid URL download.",
                    "Job id =" + str( str_job_id ) + ".",
                    "URL =" + str( str_download_url ) + "." ] )
  exit( 101 )

# Get zip file
f_success = func_download_cravat_result( str_download_url, args_call.str_output_dir )
if not f_success:
  print " ".join( [ "annotate_with_cravat::Error, could not download data in URL.",
                    "Job id =" + str( str_job_id ) + ".",
                    "URL =" + str( str_download_url ) + ".",
                    "Download =" + str( args_call.str_output_dir ) + "." ] )
  exit( 102 )
