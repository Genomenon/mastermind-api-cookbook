# -*- coding: utf-8 -*-

import json, requests, urllib, sys
from collections import defaultdict

URL = "https://mastermind.genomenon.com/api/v2/"
# Find your API token by logging in, visiting https://mastermind.genomenon.com/api, and clicking the link that says "Click here to fetch your API token".
API_TOKEN = "INSERT API TOKEN HERE"
MAX_ARTICLES = 1000

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
    """
    Call in a loop to create terminal progress bar

    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = 'M' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()

def search_mode():
    mode = raw_input("Search by:\n\t1. Known gene pair (both genes known)\n\t2. Known gene partner (one gene known)\n\t3. Known disease\n(Type 1, 2, or 3 and hit Enter): ")
    return int(mode)

def choose_suggestion(options):
    suggestions = api_get("suggestions", options)
    value = encode(suggestions[0]['canonical'])

    return value

def fusion_pmids(options):
    params = options.copy()
    params.update({'categories[]': ['fusion', 'breakpoint']})

    data = api_get("articles", params)
    pmids = [article['pmid'] for article in data['articles']]
    # Only get 2000 most relevant PMIDs
    articles = min(int(data['article_count']), MAX_ARTICLES)
    pages = min(int(data['pages']), MAX_ARTICLES/5)

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
        gene_a = choose_suggestion({'gene': raw_input("Enter gene A: ")})

    if gene_b is None:
        gene_b = choose_suggestion({'gene': raw_input("Enter gene B: ")})

    if gene_a_pmids is None:
        gene_a_pmids = fusion_pmids({'gene': gene_a})
    if gene_b_pmids is None:
        gene_b_pmids = fusion_pmids({'gene': gene_b})

    intersection = list(set(gene_a_pmids)&set(gene_b_pmids))

    pmids_by_gene_pair = {}
    pmids_by_gene_pair[gene_a + '-' + gene_b] = {'pmids': intersection}

    return pmids_by_gene_pair

def gene_search(gene_a=None):
    if gene_a is None:
        gene_a = choose_suggestion({'gene': raw_input("Enter gene: ")})

    data = api_get("genes", {'gene': gene_a})
    genes = [gene['symbol'] for gene in data['genes'] if str(gene['symbol']).upper() != str(gene_a).upper()]

    print("Found gene partner candidates for " + str(gene_a).upper() + ": " + ', '.join(genes))

    pmids_by_gene_pair = {}
    gene_a_pmids = fusion_pmids({'gene': gene_a})

    for gene_b in genes:
        pmids_by_gene_pair.update(gene_pair_search(gene_a, gene_b, gene_a_pmids))

    return pmids_by_gene_pair

def disease_search(disease=None):
    if disease is None:
        disease = choose_suggestion({'disease': raw_input("Enter disease: ")})

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
            pmids_by_gene_pair.update(gene_pair_search(gene_a, gene_b, gene_pmids[gene_a], gene_pmids[gene_b]))

    return pmids_by_gene_pair

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
    if mode == 1:
        print("Gene Pair mode selected")
        pmids_by_gene_pair = gene_pair_search()
    elif mode == 2:
        print("Gene Partner mode selected")
        pmids_by_gene_pair = gene_search()
    else:
        print("Disease mode selected")
        pmids_by_gene_pair = disease_search()

    info = aggregate_article_info(pmids_by_gene_pair)

    print('-'*100)
    for gene_pair, data in info.items():
        if len(data['pmids']) == 0:
            print ("No articles found with " + str(gene_pair).upper())
        else:
            print("Found the following associations for " + str(gene_pair).upper() + ":")

            print("PMIDs:")
            print("\t" + ', '.join(data['pmids']))

            print("\tVariants with supporting PMIDs:")
            for values in sorted(data['variants'].items(), key=lambda item: len(item[1]), reverse=True):
                print("\t\t" + str(values[0]) + ": " + ', '.join(values[1]))

            print("\tDiseases with supporting PMIDs:")
            for values in sorted(data['diseases'].items(), key=lambda item: len(item[1]), reverse=True):
                print("\t\t" + str(values[0]).title() + ": " + ', '.join(values[1]))

if __name__ == "__main__":
    main()
