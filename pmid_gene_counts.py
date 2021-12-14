#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import re
import json
import requests
import urllib
from collections import Counter
import datetime
import time
import codecs

URL = "https://mastermind.genomenon.com/api/v2/"

# Find your API token by logging in, visiting https://mastermind.genomenon.com/api, and clicking the link that says "Click here to fetch your API token".
API_TOKEN = "INSERT API TOKEN HERE"

def api_get(endpoint, options, tries=0):
    params = options.copy()
    params.update({'api_token': API_TOKEN})

    # print("Querying API: ", endpoint, options)
    response = requests.get(url=URL+endpoint, params=params)

    return json_or_print_error(response, endpoint, options, tries)

def json_or_print_error(response, endpoint, options, tries):
    if response.status_code == requests.codes.ok:
        return response.json()
    else:
        sys.stdout.write('\n')
        print("ERROR ENCOUNTERED. ERROR CODE: " + str(response.status_code))
        if response.status_code != 500 and response.text:
            print("\t" + response.text)
        print("\tRESULTING FROM REQUEST: " + endpoint)
        print("\tWITH PARAMS: " + str(options))
        if response.status_code in [404]:
            print("SKIPPING DATA FOR ABOVE REQUEST")
            return
        elif response.status_code in [408, 500]:
            if tries < 1:
                print("Time out error, trying again.")
                return api_get(endpoint, options, tries+1)
            else:
                print("SKIPPING DATA FOR ABOVE REQUEST")
                return
        sys.exit(0)

def encode(str):
    if sys.version_info[0] < 3:
        return urllib.quote_plus(str)
    else:
        return urllib.parse.quote_plus(str)

def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = 'M' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()

def pmid_list():
    pmid_filename = "braf_neoplasms_multiple_sequence_alignment.pmids.txt"

    pmids = []
    with open(pmid_filename, "r") as lines:
        pmids = lines.read().splitlines()

    print("PMIDS", pmids)
    return pmids

def get_articles(pmids):
    article_info = {}

    pmid_len = len(pmids)
    current_pmid = 1
    print_progress(current_pmid, pmid_len, prefix = 'Getting ' + str(pmid_len) + ' articles:', suffix = 'Complete', bar_length = 50)

    for pmid in pmids:
        current_pmid += 1
        print_progress(current_pmid, pmid_len, prefix = 'Getting ' + str(pmid_len) + ' articles:', suffix = 'Complete', bar_length = 50)
        if pmid in article_info:
            data = article_info[pmid]
        else:
            # Get article_info for each PMID
            data = api_get("article_info", {'pmid': pmid})
            article_info[pmid] = data

    sys.stdout.write('\n')
    sys.stdout.flush()

    return article_info

def main(args):
    pmids = pmid_list()
    article_info = get_articles(pmids)

    gene_counts_by_pmid = {pmid:len(data["genes"]) for (pmid, data) in article_info.items()}
    counter = Counter(gene_counts_by_pmid)

    print(counter.most_common)

if __name__ == "__main__":
    main(sys.argv)
