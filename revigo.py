#!/usr/bin/env python3

# Python script for programtic access to Revigo. Run it with (last output file name is optional):
# python3 revigo.py example.csv result.csv

import requests
import time
import sys

# Read enrichments file
userData = open(sys.argv[1],'r').read()

# Submit job to Revigo
payload = {'cutoff':'0.7', 'valueType':'pvalue', 'speciesTaxon':'4081', 'measure':'SIMREL', 'goList':userData}
r = requests.post("http://revigo.irb.hr/StartJob.aspx", data=payload)

jobid = r.json()['jobid']

# Check job status
running = 1
while (running!=0):
    r = requests.post("http://revigo.irb.hr/QueryJob.aspx", data={'jobid':jobid,'type':'jstatus'})
    running = r.json()['running']
    print(r.json()['progress'])
    time.sleep(5)

# Fetch results
r = requests.post("http://revigo.irb.hr/QueryJob.aspx", data={'jobid':jobid, 'namespace':'1', 'type':'Table'})

# Write results to a file - if file name is not provided the default is result.csv
with open(sys.argv[2] if len(sys.argv)>=3 else 'result.tsv','w') as f:
    f.write(r.text)

