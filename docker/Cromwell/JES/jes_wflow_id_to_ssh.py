#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os, sys, re
import subprocess
import json
import pprint
 
def main():

    usage = "\n\n\tusage: {} jes_wflow_id\n\n".format(sys.argv[0])
    if len(sys.argv) != 2:
        print(usage, file=sys.stderr)
        sys.exit(1)

    wflow_id = sys.argv[1]

    json_info = subprocess.check_output("gcloud compute instances list --format=json", shell=True)

    json_structs_aref = json.loads(json_info)

    ggp_name = get_ggp_name(json_structs_aref, wflow_id)

    print("\n\n\tgcloud compute ssh --zone=us-central1-a {}\n\n".format(gpp_name))

    sys.exit(0)


def get_gpp_name(json_structs_aref, wflow_id):
    
    # example wflow_id:
    # EOimlIPMKhjd3M_GyfH1jXQgw7vetLsXKg9wcm9kdWN0aW9uUXVldWU
    
    for json_struct in json_structs_aref:
        
        if 'description' in json_struct and \
               re.search(wflow_id, json_struct['description']):
            return(json_struct['name'])

    raise RuntimeException("Error, no gpp name entry found for wflow_id: {}".format(wflow_id))
        



 
####################
 
if __name__ == "__main__":
    main()
