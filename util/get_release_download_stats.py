#!/usr/bin/env python3

# based on: git@github.com:mmilidoni/github-downloads-count.git

import sys, os, re
import requests
import json
from datetime import date

repos = [  "trinityrnaseq/trinityrnaseq",
           "STAR-Fusion/STAR-Fusion",
           "FusionInspector/FusionInspector",
           #"broadinstitute/infercnv", get bioconductor stats from here:  # http://bioconductor.org/packages/stats/bioc/infercnv/
           "NCIP/ctat-mutations", 
           "NCIP/CTAT-SPLICING",
           "broadinstitute/CTAT-VirusIntegrationFinder",
]


def main():

    for repo in repos:
        report_download_stats(repo)
        print() # spacer between progs.


    sys.exit(0)


def report_download_stats(repo):

    r = requests.get('https://api.github.com/repos/' + repo + '/releases')
    myobj = r.json()

    today = date.today()
    
    for p in myobj:
        #print(p)
        #print()
        release_tag = p['tag_name']
        published_time = p['published_at']
        for asset in p['assets']:
            asset_url = asset['browser_download_url']
            download_count = asset['download_count']

            print("\t".join([str(today), repo, release_tag, os.path.basename(asset_url), published_time, str(download_count)]))
            



            
if __name__=='__main__':
    main()
