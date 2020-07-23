#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Discover gene fusion evidence for desired disease or gene(s).

Usage:

Before running this code, first change the API_TOKEN line below to your
Mastermind API key and save the updated file.

You may also choose to increase or decrease the DEFAULT_MAX_ARTICLES value
below. Increasing the number of articles to look through for each gene-fusion
candidate will increase the sensitivity of the found associations, but will
also increase the number of API calls used and also the amount of time it takes
to run. The maximum value is 5000.

Then run the file directly from the command line:
    ./gene_fusion_evidence.py

You will be prompted to choose the desired mode for discovery, whether you're
starting from a known disease, a known gene-candidate, or a known gene-fusion
pair. Then you'll be prompted for the known inputs based on the chosen input
mode.

You will also be prompted the desired output mode, whether you want a full
report of variants, diseases, and PMIDs associated with the found gene fusion
candidate(s), or just the list of PMIDs associated with the found gene fusion
candidate(s).
"""

import sys
import re
import json
import requests
import urllib
from collections import defaultdict
import codecs

URL = "https://mastermind.genomenon.com/api/v2/"
# Find your API token by logging in, visiting https://mastermind.genomenon.com/api, and clicking the link that says "Click here to fetch your API token".
API_TOKEN = "INSERT API TOKEN HERE"
# Maxiumum value below is 5000.
DEFAULT_MAX_ARTICLES = 1000
_sensitivity = DEFAULT_MAX_ARTICLES

def api_get(endpoint, options):
    params = options.copy()
    params.update({'api_token': API_TOKEN})

    # print("Querying API: ", endpoint, options)
    response = requests.get(url=URL+endpoint, params=params)

    return json_or_print_error(response)

def json_or_print_error(response):
    if response.status_code == requests.codes.ok:
        return response.json()
    else:
        print(response.status_code)
        print(response.text)
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

def get_mode(options):
    mode = None
    options = [str(options.index(o)+1) + ". " + o for o in options]
    # Turn number of options into text, e.g.: "1, 2 or 3"
    valid = ' or '.join( ', '.join([str(i) for i in range(1, len(options)+1)]).rsplit(', ', 1))
    while mode is None:
        entered = raw_input("\t" + '\n\t'.join(options) + "\n(Type " + valid + " and hit Enter): ")
        try:
            mode = int(entered)
            if mode > len(options):
                mode = None
                raise Exception
        except:
            print("Please enter " + valid)
    return mode

def search_mode():
    print("Search by:")
    return get_mode(["Known gene pair (both genes known)", "Known gene partner (one gene known)", "Known disease"])

def output_mode():
    print("Compile and output:")
    return get_mode(["Variants, diseases, and PMIDs", "PMIDs only"])

def choose_sensitivity():
    global _sensitivity
    print("Choose maximum number of articles to fetch for each gene. The more articles, the more complete the results, but the longer the program will take to run.")
    sensitivity_entered = None
    while sensitivity_entered is None:
        sensitivity_entered = re.sub('[,\.]', '', raw_input("Enter sensitivity (between 1 and 10,000), or leave blank for default [" + str(DEFAULT_MAX_ARTICLES) + "], and press Enter: "))
        if sensitivity_entered:
            try:
                _sensitivity = int(sensitivity_entered)
                if _sensitivity > 10000 or _sensitivity < 1:
                    raise Exception
            except:
                sensitivity_entered = None
                _sensitivity = DEFAULT_MAX_ARTICLES
                print("Please enter a valid value or leave blank.")

def choose_suggestion(suggestion_type, prompt):
    value = None
    while value is None:
        options = {}
        suggestion_input = raw_input(prompt)
        options[suggestion_type] = suggestion_input
        suggestions = api_get("suggestions", options)
        if len(suggestions) == 0:
            print(suggestion_type.title() + " could not be found for \"" + options[suggestion_type] + "\"")
            print("Please make sure you've properly typed it or try searching by another nomenclature.")
        else:
            value = encode(suggestions[0]['canonical'])

    return value

def fusion_pmids(options):
    params = options.copy()
    params.update({'categories[]': ['fusion', 'breakpoint']})

    data = api_get("articles", params)
    pmids = [article['pmid'] for article in data['articles']]

    articles = min(int(data['article_count']), _sensitivity)
    pages = min(int(data['pages']), _sensitivity/5)

    print_progress(1, pages, prefix = 'Getting ' + str(articles) + ' articles for ' + str(options['gene']).upper() + ':', suffix = 'Complete', bar_length = 50)

    if pages > 1:
        for page in range(2, pages+1):
        # for page in range(2, 3):
            print_progress(page, pages, prefix = 'Getting ' + str(articles) + ' articles for ' + str(options['gene']).upper() + ':', suffix = 'Complete', bar_length = 50)

            params.update({'page': page})
            data = api_get("articles", params)
            pmids = pmids + [article['pmid'] for article in data['articles']]

    return pmids

def gene_pair_search(gene_a=None, gene_b=None, gene_a_pmids=None, gene_b_pmids=None):
    if gene_a is None:
        gene_a = choose_suggestion("gene", "Enter gene A: ")

    if gene_b is None:
        gene_b = choose_suggestion("gene", "Enter gene B: ")

    if gene_a_pmids is None:
        gene_a_pmids = fusion_pmids({'gene': gene_a})
    if gene_b_pmids is None:
        gene_b_pmids = fusion_pmids({'gene': gene_b})

    intersection = list(set(gene_a_pmids)&set(gene_b_pmids))

    pmids_by_gene_pair = {}
    pmids_by_gene_pair[gene_a + '-' + gene_b] = {'pmids': intersection}

    return gene_a, gene_b, pmids_by_gene_pair

def gene_search(gene_a=None):
    if gene_a is None:
        gene_a = choose_suggestion("gene", "Enter gene: ")

    data = api_get("genes", {'gene': gene_a})
    genes = [gene['symbol'] for gene in data['genes'] if str(gene['symbol']).upper() != str(gene_a).upper()]

    print("Found gene partner candidates for " + str(gene_a).upper() + ": " + ', '.join(genes))

    pmids_by_gene_pair = {}
    gene_a_pmids = fusion_pmids({'gene': gene_a})

    for gene_b in genes:
        _a, _b, results = gene_pair_search(gene_a, gene_b, gene_a_pmids)
        pmids_by_gene_pair.update(results)

    return gene_a, pmids_by_gene_pair

def disease_search(disease=None):
    if disease is None:
        disease = choose_suggestion("disease", "Enter disease: ")

    data = api_get("genes", {'disease': disease})
    genes = [gene['symbol'] for gene in data['genes']]
    total_genes = len(genes)

    print("Found gene partner candidates for " + str(disease).title() + ": " + ', '.join(genes))

    gene_pmids = {}
    for gene in genes:
        gene_pmids[gene] = fusion_pmids({'gene': gene})

    pmids_by_gene_pair = defaultdict(lambda: {})

    # Iterate over each unique pair
    for gene_a in genes:
        i = genes.index(gene_a)
        for l in range(i+1, total_genes):
            gene_b = genes[l]
            _a, _b, results = gene_pair_search(gene_a, gene_b, gene_pmids[gene_a], gene_pmids[gene_b])
            pmids_by_gene_pair.update(results)

    return disease, pmids_by_gene_pair

def aggregate_article_info(pmids_by_gene_pair):
    pmid_info = {}

    for gene_pair, values in pmids_by_gene_pair.items():
        diseases_by_pmids = defaultdict(lambda: [])
        variants_by_pmids = defaultdict(lambda: [])

        current = 0
        total = len(values['pmids'])

        for pmid in values['pmids']:
            current += 1
            print_progress(current, total, prefix = 'Inspecting PMID info for ' + str(gene_pair).upper() + ':', suffix = 'Complete', bar_length = 50)

            if pmid in pmid_info:
                data = pmid_info[pmid]
            else:
                # Get article_info for each PMID
                data = api_get("article_info", {'pmid': pmid})
                pmid_info[pmid] = data

            if 'diseases' in data:
                for disease in data['diseases']:
                    diseases_by_pmids[disease['key']].append(pmid)

            for gene in data['genes']:
                if 'variants' in gene:
                    for variant in gene['variants']:
                        variants_by_pmids[gene['symbol'] + ':' + variant['key']].append(pmid)

        pmids_by_gene_pair[gene_pair]['diseases'] = diseases_by_pmids
        pmids_by_gene_pair[gene_pair]['variants'] = variants_by_pmids

    return pmids_by_gene_pair

def main():
    print("Welcome to the Gene Fusion Evidence program powered by Mastermind.")
    mode = search_mode()
    output = output_mode()
    choose_sensitivity()
    if mode == 1:
        gene_a, gene_b, info = gene_pair_search()
        filename_prefix = gene_a + "-" + gene_b
    elif mode == 2:
        gene, info = gene_search()
        filename_prefix = gene
    else:
        disease, info = disease_search()
        filename_prefix = disease

    if output == 1:
        info = aggregate_article_info(info)

    output_filename = "-".join([filename_prefix, "gene-fusions", str(_sensitivity), "article-sensitivity"]) + ".txt"
    with codecs.open(output_filename, 'wb', 'utf-8') as output_file:
        for gene_pair, data in info.items():
            if len(data['pmids']) == 0:
                output_file.write("No articles found with " + str(gene_pair).upper() + "\n")
            else:
                output_file.write("Found the following associations for " + str(gene_pair).upper() + ":\n")

                output_file.write("PMIDs:\n")
                output_file.write("\t" + ', '.join(data['pmids']) + "\n")

                if output == 1:
                    output_file.write("\tVariants with supporting PMIDs:\n")
                    for values in sorted(data['variants'].items(), key=lambda item: len(item[1]), reverse=True):
                        output_file.write("\t\t" + str(values[0]) + ": " + ', '.join(values[1]) + "\n")

                    output_file.write("\tDiseases with supporting PMIDs:\n")
                    for values in sorted(data['diseases'].items(), key=lambda item: len(item[1]), reverse=True):
                        output_file.write("\t\t" + str(values[0]).title() + ": " + ', '.join(values[1]) + "\n")

if __name__ == "__main__":
    main()
