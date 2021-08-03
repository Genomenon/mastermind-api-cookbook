#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Annotate VCF and annotated VCF files with Mastermind File Annotations API.

Usage:

Before running this code, first change the API_TOKEN line below to your
Mastermind API key and save the updated file.

Then run the file directly from the command line, passing in the local file
path to the file you'd like to annotate:
    ./annotate_with_evidence.py "/path/to/file.vcf.gz"

This will do the following:
1. Create a new Mastermind API file annotation job and assign a unique job ID.
2. Upload the local gzip-compressed VCF file.
3. Poll the Mastermind API and wait for the annotation job to complete.
4. Download the gzip-compressed annotated VCF file to the local filesystem (the
   file will be saved alongside the uploaded file on the local filesystem).

This script can also continue and download a previously created job by passing
in the previously assigned Job ID as a second argument:
    ./annotate_with_evidence.py "/path/to/file.vcf.gz" "8a6830ab-2051-409c-b7d3-ad609bc70d49"
"""

import sys
import os
import re
import time
import json
import requests

URL = "https://mastermind.genomenon.com/api/v2/"
# Find your API token by logging in, visiting https://mastermind.genomenon.com/api, and clicking the link that says "Click here to fetch your API token".
API_TOKEN = "INSERT API TOKEN HERE"
ASSEMBLY = "grch37"

def api_request(endpoint, options, request_type="GET", json_request=True):
    params = options.copy()
    params.update({'api_token': API_TOKEN})

    # print("Querying API: ", endpoint, options)
    response = requests.request(request_type, url=URL+endpoint, params=params)

    if json_request:
        return json_or_print_error(response)
    else:
        return response

def json_or_print_error(response):
    if response.status_code == requests.codes.ok:
        return response.json()
    else:
        print(response.status_code)
        print(response.text)
        sys.exit(0)

def wait_for_success(job_url):
    message = "Waiting for processing success"

    response = api_request(job_url, {})
    state = response['state']

    i = 3
    sys.stdout.write('\r' + message)
    sys.stdout.flush()
    while state == "created" or state == "started":
        sys.stdout.write('\r' + message + i*'.')
        sys.stdout.flush()
        time.sleep(5)
        i += 1
        response = api_request(job_url, {})
        state = response['state']

    sys.stdout.write("\n")
    return response, state

def main(args):
    input_vcf_path = args[1]
    input_vcf_name = os.path.basename(input_vcf_path)

    if len(args) == 3:
        job_id = args[2]
        print("Continuing job")
        response = api_request("file_annotations/counts/" + job_id, {})
    else:
        job_id = None

    # 1. Create new file annotation job
    if job_id == None:
        print("Creating file annotation job for file: " + input_vcf_name)
        response = api_request('file_annotations/counts', {'assembly': ASSEMBLY, 'filename': input_vcf_name}, 'POST')

    job_id = response['job_id']
    job_url = response['job_url']
    upload_url = response['upload_url']
    state = response['state']
    print("Job ID: " + job_id)

    if state == "created":
        # 2. Upload VCF or annotated VCF file for annotation
        print("Uploading file...")
        with open(input_vcf_path, 'rb') as f:
            requests.put(upload_url, data=f, headers={'Content-Type': 'application/octet-stream'})

    if state == "created" or state == "started":
        # 3. Periodically check status until job is finished
        response, state = wait_for_success("file_annotations/counts/" + job_id)

    if state == "succeeded":
        print("Successfully annotated. " + str(response['annotated']) + " out of " + str(response['records']) + " have annotations.")
    else:
        print("There was an error processing the file for Job ID " + job_id)
        sys.exit(1)

    download_url = response['download_url']

    # 4. Download annotated file
    output_file_path = re.sub(r"\.vcf(\.gz)?$", ".annotated-" + job_id + ".vcf.gz", input_vcf_path)
    print("Downloading annotated file to " + output_file_path)
    downloaded_file = api_request("file_annotations/counts/" + job_id + "/download", {}, json_request=False)
    with open(output_file_path, 'wb') as output_file:
        output_file.write(downloaded_file.content)

if __name__ == "__main__":
    main(sys.argv)
