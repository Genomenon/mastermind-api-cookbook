#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Given a gene (or list of genes), produce a list of all the variants and their article counts.

Examples:

# Get all JAK1 variants and print CSV-style output to the console
$ python3 all_gene_variants.py --gene JAK1

# Same as above, but write results to a CSV file:
$ python3 all_gene_variants.py --gene JAK1 --outputcsv jak1_variants.csv

# Lookup all the variants of a list of genes (a text file with one-gene-per line) as input, writing to a JSON file:
$ python3 all_gene_variants.py --inputfile numerous_genes.txt --outputjson my-results.json

usage: all_gene_variants.py [-h] [-f INPUTFILE | -g GENE] [-c OUTPUTCSV] [-j OUTPUTJSON]

optional arguments:
  -h, --help            show this help message and exit
  -f INPUTFILE, --inputfile INPUTFILE
                        File containing a list of genes to query variants for
  -g GENE, --gene GENE  The gene to lookup variants for
  -c OUTPUTCSV, --outputcsv OUTPUTCSV
                        Name of the JSON output file to write
  -j OUTPUTJSON, --outputjson OUTPUTJSON
                        Name of the CSV output file to write
"""

import sys
import urllib
import requests
import json
import argparse

URL = "https://mastermind.genomenon.com/api/v2/"
# Find your API token by logging in, visiting https://mastermind.genomenon.com/api, and clicking the link that says "Click here to fetch your API token".
API_TOKEN = "GSFljIDENzo8rdhcPvZonkBJbc9kyoWK7TkeVf6DvvnL92AuL8esBAQhlCMR6zi2"

# result sets are paginated
MAX_PAGE_COUNT = 1000


def json_or_print_error(response):
    if response.status_code == requests.codes.ok:
        return response.json()
    else:
        print(response.status_code)
        print(response.text)
        sys.exit(0)


def download_all_variants(gene):
    variants = []
    page = 1
    last_page = MAX_PAGE_COUNT
    while page <= last_page:
        print(f"Gene {gene} variants - requesting page {page}")
        response = requests.get(
            url=URL + "variants",
            params={"api_token": API_TOKEN, "gene": gene, "page": page},
        )
        data = json_or_print_error(response)
        last_page = min(data["pages"], MAX_PAGE_COUNT)
        variants.extend(data["variants"])
        page += 1

    print(f"Got {len(variants)} variants for gene {gene}")
    return variants


def write_variant_csv(outf, variants):
    outf.write("Gene, Variant, # Citations, Link\n")
    for variant in variants:
        outf.write("{gene}, {key}, {article_count}, {url}\n".format(**variant))


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-f",
        "--inputfile",
        default=None,
        help="File containing a list of genes to query variants for",
    )
    group.add_argument(
        "-g", "--gene", default=None, help="The gene to lookup variants for"
    )
    parser.add_argument(
        "-c", "--outputcsv", default=None, help="Name of the JSON output file to write"
    )
    parser.add_argument(
        "-j", "--outputjson", default=None, help="Name of the CSV output file to write"
    )
    args = parser.parse_args()

    if args.inputfile:
        # Read all the genes from the input file and get variants for all of them:
        variants = []
        with open(args.inputfile) as inf:
            for line in inf.readlines():
                gene = line.strip()
                if len(gene) > 0:
                    print(f"Gene: {gene}")
                    vs = download_all_variants(gene)
                    variants.extend(vs)
    elif args.gene:
        # Lookup variants for the given gene:
        variants = download_all_variants(args.gene)
    else:
        # Prompt for gene, Lookup variants for the given gene:
        gene = input("Type gene (e.g. UCP3): ")
        variants = download_all_variants(gene)

    did_output = False
    if args.outputcsv:
        # Write CSV file:
        with open(args.outputcsv, "w") as outf:
            write_variant_csv(outf, variants)
        print(f"Wrote {args.outputcsv}")
        did_output = True

    if args.outputjson:
        # Write JSON file:
        with open(args.outputjson, "w") as outf:
            json.dump(variants, outf, indent=4)
        print(f"Wrote {args.outputjson}")
        did_output = True

    if not did_output:
        # Print CSV rows to stdout
        write_variant_csv(sys.stdout, variants)
        did_output = True


if __name__ == "__main__":
    main()
