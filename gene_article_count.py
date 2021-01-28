#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Using the Mastermind API, produce a CSV-formatted output with article
counts for each gene in the input file.

Usage:

Before running this code, first change the API_TOKEN line below to your
Mastermind API key and save the updated file.

Then run the file directly from the command line, replacing the first two
arguments with the path to the input genes file and the desired name of the
output CSV file.

    ./gene_article_count.py /path/to/genes.txt results.csv

As the program runs, status is printed to stdout, including any errors that happen.

Budget about ~1.5-3 seconds per gene (each gene is an API call to /v2/counts)
"""

import sys
import re
import urllib
import requests
import json

URL = "https://mastermind.genomenon.com/api/v2/"
# Find your API token by logging in, visiting https://mastermind.genomenon.com/api, and clicking the link that says "Click here to fetch your API token".
# API_TOKEN = "INSERT API TOKEN HERE"
API_TOKEN = "GSFljIDENzo8rdhcPvZonkBJbc9kyoWK7TkeVf6DvvnL92AuL8esBAQhlCMR6zi2" # dcrosby deleteme

def api_get(endpoint, options):
    params = options.copy()
    params.update({'api_token': API_TOKEN})
    response = requests.get(url=URL+endpoint, params=params)
    return json_or_print_error(response)

def json_or_print_error(response):
    if response.status_code == requests.codes.ok:
        return response.json()
    elif response.status_code == requests.codes.not_found:
        return {'article_count': 0}
    else:
        print(response.status_code)
        print(response.text)
        print(response)
        sys.exit(0)

def main(args):

    if len(args) < 3:
        sys.exit("USAGE: python gene_article_count.py <infile> <outfile>")

    infile = args[1]
    outfile = args[2]

    print("Reading gene list from input file: %s"%infile)
    print("Writing counts as CSV to: %s"%outfile)
    with open(outfile,"w") as outf:
        with open(infile, "r") as lines:
            for line in lines:
                gene = line.strip()
                print("%s: querying article count..."%gene)
                try:
                    data = api_get("counts", {'gene':gene})
                    count = data["article_count"]
                    outf.write("%s,%d\n"%(gene,count))
                except Exception as e:
                    print("ERROR getting counts for %s"%gene,e)

if __name__ == "__main__":
    main(sys.argv)
