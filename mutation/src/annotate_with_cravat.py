#!/usr/bin/env python

import argparse
import json
import requests
import time
import urllib

# Constants
str_response_job_id = "jobid"
str_response_status = "status"
str_response_success = "status"
str_response_url = "resultfileurl"
str_response_pass = u'Success'
str_response_fail = u'submissionfailed'

# for testing delete
str_email = "ttickle@broadinstitute.org"
str_classifier = "Blood-Lymphocyte"
f_hg_18 = True

# Wait time in seconds
i_wait = 10
i_max_attempts = 100

# Parse arguments
prsr_arguments = argparse.ArgumentParser( prog = "annotate_with_cravat.py", description = "Annotate VCF file with cravat.", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( c( "--classifier" ), dest = "str_classifier", action = "store", required=True, help = "Tissue type." )
prsr_arguments.add_argument( c( "--is_hg18" ), dest = "f_hg_18", action = "store", require=True help = "Indicates if Hg18 (True=Hg18; False=Hg19)." )
prsr_arguments.add_argument( c( "--email" ), dest = "str_email", action = "store", required=True, help = "Email contact for job." )
prsr_arguments.add_argument( dest = "str_input_file", action = "store", help = "VCF file to update." )
prsr_arguments.add_argument( dest = "str_output_dir", action = "store", help = "Location to download directory." )
args_call = prsr_arguments.parse_args()

# Runs
def func_check_job( str_json_id ):
  # Send file over to service
  response_check_job = requests.get( "http://www.cravat.us/rest/service/status", params={ str_response_job_id: str_json_id } )
  return response_check_job.json()

# Runs
def func_download_cravat_result( str_cravat_url, str_download_location ):
  # Download Cravat result
  try:
    urllib.urlretrieve( str_cravat_url, str_download_location )
    return True
  except IOError:
    print "annotate_with_cravat::Error in retrieving results. " + str( IOError )
    return False

# Runs
def func_get_cravat_response( str_json_id ):
  for x_attempt in xrange( i_max_attempts ):
    print "annotate_with_cravat::Attempting to get CRAVAT response. Attempt " + str( x_attempt )
    time.sleep( i_wait )
    response_json = func_check_job( str_json_id )
    if not response_json:
      print "annotate_with_cravat::Error did not get a resposne from CRAVAT"
      return( None )
    str_success = str( response_json.get( str_response_status, None) )
    if ( not str_success ) or ( str_success == str_response_fail ):
      print " ".join( [ "annotate_with_cravat::Error was not successful in getting info from CRAVAT.",
                        "Job id =" + str( str_json_id ) + ".",
                        "Success =" + str( str_success ) + "." ] )
      return( None )
    if str_success == str_response_pass:
      return response_json.get( str_response_url, None )
  print " ".join( [ "annotate_with_cravat::Error gave up on CRAVAT after " + str( i_max_attempts * i_wait ) + " seconds.",
                    "Job id =" + str( str_json_id ) + ".",
                    "Success =" + str( str_success ) + "." ] )
  return( None )

# Runs 
def func_read_vcf( str_path ):
  # Read in file into a string
  str_initial_vcf = ""
  with open( str_path, "rb" ) as hndl_in:
    str_initial_vcf = hndl_in.read()
  return str_initial_vcf

# Runs
def func_request_cravat_service( str_vcf_path, str_classifier, f_hg_18, str_email ):
  # Read in file into a string
  str_initial_vcf = func_read_vcf( str_vcf_path )
  if not str_initial_vcf:
    print " ".join( [ "annotate_with_cravat::Error did not read VCF file.",
                      "Path =" + str( str_vcf_path ) + "." ] )
    return( None )

  # Encode request info
  pyld_request = { "chasmclassifier": str( str_classifier ),
                   "mupitinput": "on",
                   "hg18": "on" if f_hg_18 else "off",
                   "tsvreport": "on",
                   "analyses": "CHASM%3BSnvGet",
                   "functionalannotation": "on",
                   "analysistype": "driver",
                   "email": str_email,
                   "mutations": str_initial_vcf }

  # Send file over to service
  response_cravat = requests.get( "http://www.cravat.us/rest/service/submit", params=pyld_request )
  print( response_cravat.url )
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
str_download_url = func_get_cravat_response( str_job_id )
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
  
# Update VCF file
# ?
