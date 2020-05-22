#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Get evidence citation count for variant from HGVS genomic coordinate
description from Mastermind API.

Usage:

Before running this code, first change the API_TOKEN line below to your
Mastermind API key and save the updated file.

Then run the file directly from the command line:
    ./genomic_variant_citation_count.py

You will be prompted to type in the HGVS genomic coordinate description of the
desired variant and hit Enter.
"""

import sys

import base

def main():

    variant = raw_input("Type genomic variant (e.g. NC_000012.11:g.57489193T>C): ")

    data = base.api_request("suggestions", {'variant': variant})
    encoded_match = base.encode(data[0]['matched'])

    print("Matched: " + data[0]['matched'])
    print("Mastermind link: https://mastermind.genomenon.com/detail?mutation=" + encoded_match)
    print("Maps to " + str(len(data)) + " effects:")

    for match in data:
        canonical = match['canonical']
        count_data = base.api_request("counts", {'variant': canonical})
        print("\t" + canonical + " - " + str(count_data['article_count']) + " articles, direct link: " + count_data['url'])

if __name__ == "__main__":
    main()
