#!/usr/bin/env python

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2015"
__credits__ = ["Timothy Tickle", "Brian Haas"]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@broadinstitute.org"
__status__ = "Development"

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
lstr_response_fail = [u'Error', u'submissionfailed']

# VCF
str_vcf_comment = "#"
str_cravat_header = "\t".join(["# UID", "Chr.", "Position",
                               "Strand", "Ref. base", "Alt. base"])
c_VCF_delimiter = "\t"
i_vcf_chr = 0
i_vcf_position = 1
i_vcf_ref = 3
i_vcf_alt = 4

# Parse arguments
prsr_arguments = argparse.ArgumentParser(prog = "annotate_with_cravat.py", description = "Retrieve annotations associated with variants in a VCF file using CRAVAT.", formatter_class = argparse.ArgumentDefaultsHelpFormatter)
prsr_arguments.add_argument("--analysis", dest = "str_analysis", action = "store", default="CHASM;VEST;", help = "Analysis type parameter. For options see cravat.us/help.jsp")
prsr_arguments.add_argument("--classifier", dest = "str_classifier", action = "store", required=True, help = "Tissue type for the classifier. For options see cravat.us/help.jsp.")
prsr_arguments.add_argument("--email", dest = "str_email", action = "store", required=True, help = "Email contact for job.")
prsr_arguments.add_argument("--is_hg19", dest = "f_hg_19", action = "store_true", default=False, help = "Indicates the reference is Hg19 (By default assumed to be Hg38).")
prsr_arguments.add_argument("--max_attempts", dest = "i_max_attempts", default=100, action = "store", type=int, help = "Max attempts of querying response before timing out.")
prsr_arguments.add_argument("--wait", dest = "i_wait", action = "store", default=10, type=int, help = "Wait in seconds before querying the response.")
prsr_arguments.add_argument(dest = "str_input_file", action = "store", help = "Path to VCF file containing variants.")
prsr_arguments.add_argument(dest = "str_output_dir", action = "store", help = "Output Zip file (please use the extension .zip) Can be opened with 'unzip'.")
args_call = prsr_arguments.parse_args()

def func_check_job(str_json_id):
    """
    Uses the CRAVAT status web service to check if a job is complete.

    * str_json_id : CRAVAT job id.
                  : string
    """

    # Send file over to service
    response_check_job = requests.get("http://www.cravat.us/CRAVAT/rest/service/status", params={ str_response_job_id: str_json_id })
    return response_check_job.json()

def func_download_cravat_result(str_cravat_url, str_download_location):
    """
    Using the download URL from CRAVAT, downloads to the given location.

    * str_cravat_url : Url provided by CRAVAT fro results download.
                     : string
    * str_download_location : Path (including file name) to write the download.
                            : string
    """

    # Download Cravat result
    try:
        urllib.urlretrieve(str_cravat_url, str_download_location)
        return True
    except IOError:
        print "annotate_with_cravat::Error in retrieving results. " + str(IOError)
        return False

def func_get_cravat_response(str_json_id, i_max_attempts, i_wait):
    """
    Checks CRAVAT for a period of time until the success, failure, or time out are given for a job id.

    *  str_json_id : CRAVAT job id.
                   : string
    * i_max_attempts : Max times to check CRAVAT for job completion.
                     : int
    * i_wait : Time in seconds to wait between checking CRAVAT for job completion.
             : int
    """

    for x_attempt in xrange(i_max_attempts):
        print "annotate_with_cravat::Attempting to get CRAVAT response. Attempt " + str(x_attempt + 1)
        time.sleep(i_wait)
        response_json = func_check_job(str_json_id)
        if not response_json:
            print "annotate_with_cravat::Error did not get a response from CRAVAT"
            return(None)
        str_success = str(response_json.get(str_response_status, None))
        if (not str_success) or (str_success in lstr_response_fail):
            print " ".join(["annotate_with_cravat::Error was not successful in getting info from CRAVAT.",
                              "Job id =" + str(str_json_id) + ".",
                              "Success =" + str(str_success) + ".",
                              response_json.get(str_response_log, "")])
            return(None)
        if str_success == str_response_pass:
            print "annotate_with_cravat:: Reponse URL on success: " + str(response_json)
            return response_json.get(str_response_url, None)
    print " ".join(["annotate_with_cravat::Error gave up on CRAVAT after " + str(i_max_attempts * i_wait) + " seconds.",
                      "Job id =" + str(str_json_id) + ".",
                      "Success =" + str(str_success) + "."])
    return(None)


def func_request_cravat_service(str_vcf_path, str_analysis,
                                str_classifier, f_hg_19, str_email):
    """
    Request a job to occur with CRAVAT.

    * str_vcf_path : Path to VCF file of mutations to submit.
                   : string
    * str_classifier : Tissue type analyzed(see http://www.cravat.us/help.jsp).
                     : string
    * f_hg_19 : True indictes HG19 / False indicates HG38 (Current CRAVAT default).
              : boolean
    * str_email : Email of user (errors will be sent by email).
                : string
    """

    # Encode request info
    pyld_request = { "classifier": str(str_classifier),
                     "mupitinput": "on",
                     #"hg18": "off",
                     "hg19": "on", #if f_hg_19 else "off",
                     "tsvreport": "on",
                     "analyses": str_analysis,
                     "functionalannotation": "on",
                     "analysistype": "driver",
                     "email": str_email }

    # Send file over to service
    response_cravat = requests.post("http://cravat.us/CRAVAT/rest/service/submit",
                                    files={"inputfile": open(str_vcf_path)},
                                    data=pyld_request)
    print response_cravat
    json_response = response_cravat.json()
    # Get return (job id)
    return json_response.get(str_response_job_id, None)

# Ensure the extention to the output directory is zip
if not os.path.splitext(args_call.str_output_dir)[1] == ".zip":
    args_call.str_output_dir = args_call.str_output_dir + ".zip"

# Request CRAVAT service
str_job_id = func_request_cravat_service(args_call.str_input_file,
                                          args_call.str_analysis,
                                          args_call.str_classifier,
                                          args_call.f_hg_19,
                                          args_call.str_email)
if not str_job_id:
    print " ".join(["annotate_with_cravat::Error Job id not found.",
                    "Job id =" + str(str_job_id)])
    exit(100)

# Wait for result
str_download_url = func_get_cravat_response(str_job_id,
                                            args_call.i_max_attempts,
                                            args_call.i_wait)
if not str_download_url:
    print " ".join(["annotate_with_cravat::Error did no",
                    "recieve valid URL download.",
                    "Job id =" + str(str_job_id) + ".",
                    "URL =" + str(str_download_url) + "."])
    exit(101)

# Get zip file
f_success = func_download_cravat_result(str_download_url,
                                        args_call.str_output_dir)
print "annotate_with_cravat:: Success = " + str(f_success)
if not f_success:
    print " ".join(["annotate_with_cravat::Error,",
                    "could not download data in URL.",
                    "Job id =" + str(str_job_id) + ".",
                    "URL =" + str(str_download_url) + ".",
                    "Download =" + str(args_call.str_output_dir) + "."])
    exit(102)
