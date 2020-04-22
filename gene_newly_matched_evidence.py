#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Discover articles and resulting association evidence recently matched for a
given list of genes within a list of high-yield journals.

Usage:

Before running this code, first change the API_TOKEN line below to your
Mastermind API key and save the updated file.

This will find the newly matched evidence for the given list of genes, where
the list of genes are provided in a text file with one gene symbol per line.
For example, you may have a file my_gene_list.txt, which contains:

BRAF
KRAS

Then run the file directly from the command line:
    ./gene_newly_matched_evidence.py my_gene_list.txt

For my_gene_list.txt as an input, this will generate 4 output CSV files:
    my_gene_list.articles.csv
    my_gene_list.diseases.csv
    my_gene_list.phenotypes.csv
    my_gene_list.genes.csv

And one output summary file:
    my_gene_list.associations.txt

By default this looks for newly matched findings from a list of high-yield
journals which can be modified in the JOURNALS setting below.

This also defaults to articles newly matched for the genes within the past 30
days, but any date may be specified as the BEGIN setting below after Nov 26,
2019.

Another END setting may be optionally specified to filter articles matched
before a given date, which again must be after Nov 26, 2019. Note that this
will cause the script to take roughly twice as long to run, since queries must
be run twice to subtract the more recent results occurring after END from the
overall result set.

An additional specificity filter may be enabled by setting ONLY_VARIANTS =
True. This will only match articles which contain at least one genetic variant.

If testing the script or data output, the STOP_AFTER setting may be set to some
small number (typically 1-5), so that the script will exit and generate the
output data files after that many genes with non-zero results are found.
"""

import sys
import re
import json
import requests
import urllib
from collections import defaultdict
import datetime
import time
import codecs

URL = "https://mastermind.genomenon.com/api/v2/"

# Find your API token by logging in, visiting https://mastermind.genomenon.com/api, and clicking the link that says "Click here to fetch your API token".
API_TOKEN = "INSERT API TOKEN HERE"

# Whether or not to use the JOURNALS list to filter results
FILTER_JOURNALS = False

# Filter by articles from the following journals using ISO 4 abbreviated identifiers:
JOURNALS = ["Hum. Mutat.", "Am. J. Hum. Genet.", "J. Med. Genet.", "Nat. Genet.", "Hum. Genet.", "Hum. Mol. Genet.", "Eur. J. Hum. Genet.", "Clin. Genet.", "Genet. Med.", "PLOS ONE"]

# Filter by articles first matching input genes since this date.
# Defualt is 30 days before today, but can also specify a specific date:
BEGIN = datetime.date.today() - datetime.timedelta(days=14) #datetime.datetime(2019,11,26)

# Whether or not to filter out articles published before BEGIN date
FILTER_PUBLISHED_DATE = False #BEGIN - datetime.timedelta(days=60)

# Filter by articles first matching input genes before this date (set to None for no END date filtering):
END = None #datetime.datetime(2019,11,1)

# Filter by articles containing at least one variant:
ONLY_VARIANTS = False

# For testing out, you may want to stop the script after the first X genes containing data are found.
# If so, set this to the number of genes with data you'd like to stop after.
# Input genes which return no data won't count toward this number:
STOP_AFTER = False #5

def api_get(endpoint, options):
    params = options.copy()
    params.update({'api_token': API_TOKEN})

    # print("Querying API: ", endpoint, options)
    response = requests.get(url=URL+endpoint, params=params)

    return json_or_print_error(response)

def filtered_params(gene, since):
    filter_params = {'gene': gene, 'since': int(time.mktime(since.timetuple()))}
    if FILTER_JOURNALS:
        filter_params['journals[]'] = JOURNALS
    return filter_params

def skip_article(pmid_data, all_genes):
    return filter_after_published_date(pmid_data) or filter_no_variants(pmid_data, all_genes)

def filter_after_published_date(pmid_data):
    return FILTER_PUBLISHED_DATE and datetime.datetime.strptime(pmid_data['publication_date'], '%Y-%m-%d').date() < FILTER_PUBLISHED_DATE

def filter_no_variants(pmid_data, all_genes):
    if not ONLY_VARIANTS:
        return False

    for pmid_gene in pmid_data['genes']:
        if pmid_gene['symbol'].lower() in all_genes and 'variants' in pmid_gene and len(pmid_gene['variants']) > 0:
            return False

    return True

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

def get_articles(options):
    data = api_get("articles", options)
    if "articles" in data:
        pmids = [article['pmid'] for article in data['articles']]

        articles = int(data['article_count'])
        pages = int(data['pages'])

        date_string = datetime.datetime.fromtimestamp(options["since"]).strftime('%x')
        print_progress(1, pages, prefix = 'Getting ' + str(articles) + ' articles for ' + str(options['gene']).upper() + ' since ' + date_string + ':', suffix = 'Complete', bar_length = 50)

        if pages > 1:
            for page in range(2, pages+1):
            # for page in range(2, 3):
                print_progress(page, pages, prefix = 'Getting ' + str(articles) + ' articles for ' + str(options['gene']).upper() + ' since ' + date_string + ':', suffix = 'Complete', bar_length = 50)

                options.update({'page': page})
                data = api_get("articles", options)
                pmids = pmids + [article['pmid'] for article in data['articles']]
    else:
        pmids = []

    return pmids

def aggregate_article_info(gene_info, all_genes):
    article_info = {}
    disease_info = defaultdict(lambda: defaultdict(lambda: set([])))
    phenotype_info = defaultdict(lambda: defaultdict(lambda: set([])))

    for gene, values in gene_info.items():
        pmids_by_disease = defaultdict(lambda: [])
        pmids_by_phenotype = defaultdict(lambda: [])
        pmids_by_gene = defaultdict(lambda: defaultdict(lambda: []))
        pmids_by_variant = defaultdict(lambda: defaultdict(lambda: []))

        current = 0
        total = len(values['pmids'])

        gene_info[gene]['filtered_pmids'] = []

        for pmid in values['pmids']:
            current += 1
            print_progress(current, total, prefix = 'Inspecting PMID info for ' + str(gene).upper() + ':', suffix = 'Complete', bar_length = 50)

            if pmid in article_info:
                data = article_info[pmid]
                process_article = False
            else:
                # Get article_info for each PMID
                data = api_get("article_info", {'pmid': pmid})
                article_info[pmid] = data

                process_article = True
                article_info[pmid]["matched_genes"] = []
                article_info[pmid]["other_genes"] = []
                article_info[pmid]["matched_gene_variants"] = []
                article_info[pmid]["other_gene_variants"] = []


            if skip_article(data, all_genes):
                continue

            gene_info[gene]['filtered_pmids'].append(pmid)

            if 'diseases' in data:
                for disease in data['diseases']:
                    pmids_by_disease[disease['key']].append(pmid)

                    disease_info[disease['key']]['pmids'].add(pmid)

                    for other_disease in data['diseases']:
                        if other_disease['key'] != disease['key']:
                            disease_info[disease['key']]["diseases"].add(other_disease['key'])

            if 'hpo_terms' in data:
                for phenotype in data['hpo_terms']:
                    pmids_by_phenotype[phenotype['term']].append(pmid)
                    phenotype_info[phenotype['term']]['pmids'].add(pmid)
                    phenotype_info[phenotype['term']]['id'] = phenotype['key']

                    if 'diseases' in data:
                        for disease in data['diseases']:
                            disease_info[disease['key']]["phenotypes"].add(phenotype['term'])
                            phenotype_info[phenotype['term']]['diseases'].add(disease['key'])

                    for other_phenotype in data['hpo_terms']:
                        if other_phenotype['term'] != phenotype['term']:
                            phenotype_info[phenotype['term']]["phenotypes"].add(other_phenotype['term'])

            for pmid_gene in data['genes']:
                if pmid_gene['symbol'].lower() in all_genes:
                    if process_article:
                        article_info[pmid]['matched_genes'].append(pmid_gene['symbol'])

                    if 'hpo_terms' in data:
                        for phenotype in data['hpo_terms']:
                            phenotype_info[phenotype['term']]["matched_genes"].add(pmid_gene['symbol'])

                    if 'diseases' in data:
                        for disease in data['diseases']:
                            disease_info[disease['key']]['matched_genes'].add(pmid_gene['symbol'])
                else:
                    if process_article:
                        article_info[pmid]['other_genes'].append(pmid_gene['symbol'])

                    if 'hpo_terms' in data:
                        for phenotype in data['hpo_terms']:
                            phenotype_info[phenotype['term']]["other_genes"].add(pmid_gene['symbol'])

                    if 'diseases' in data:
                        for disease in data['diseases']:
                            disease_info[disease['key']]['other_genes'].add(pmid_gene['symbol'])

                if pmid_gene['symbol'].lower() != gene:
                    if pmid_gene['symbol'].lower() in all_genes:
                        pmids_by_gene['matched_genes'][pmid_gene['symbol']].append(pmid)
                    else:
                        pmids_by_gene['other_genes'][pmid_gene['symbol']].append(pmid)

                if 'variants' in pmid_gene:
                    for variant in pmid_gene['variants']:
                        variant_name = pmid_gene['symbol'] + ':' + variant['key']

                        if pmid_gene['symbol'].lower() in all_genes:
                            if process_article:
                                article_info[pmid]['matched_gene_variants'].append(variant_name)

                            if pmid_gene['symbol'].lower() == gene:
                                pmids_by_variant['gene_variants'][variant_name].append(pmid)
                            else:
                                pmids_by_variant['matched_gene_variants'][variant_name].append(pmid)

                            if 'hpo_terms' in data:
                                for phenotype in data['hpo_terms']:
                                    phenotype_info[phenotype['term']]["matched_gene_variants"].add(variant_name)

                            if 'diseases' in data:
                                for disease in data['diseases']:
                                    disease_info[disease['key']]['matched_gene_variants'].add(variant_name)

                        else:
                            if process_article:
                                article_info[pmid]['other_gene_variants'].append(variant_name)

                            pmids_by_variant['other_gene_variants'][variant_name].append(pmid)

                            if 'hpo_terms' in data:
                                for phenotype in data['hpo_terms']:
                                    phenotype_info[phenotype['term']]["other_gene_variants"].add(variant_name)

                            if 'diseases' in data:
                                for disease in data['diseases']:
                                    disease_info[disease['key']]['other_gene_variants'].add(variant_name)

        gene_info[gene]['diseases'] = pmids_by_disease
        gene_info[gene]['phenotypes'] = pmids_by_phenotype
        gene_info[gene]['genes'] = pmids_by_gene
        gene_info[gene]['variants'] = pmids_by_variant

    return gene_info, article_info, disease_info, phenotype_info

def pipe_delimited_field(values):
    return "\"" + re.sub(r"\"", "\"\"", "|".join(values)) + "\""

def main(args):
    print("Welcome to the Gene New Evidence Alerts program powered by Mastermind.")

    filename = args[1]
    genes_with_articles = 0
    gene_info = {}

    with open(filename, "r") as lines:
        # Loop through lines in gene input file
        for line in lines:
            if STOP_AFTER and genes_with_articles > STOP_AFTER:
                break
            gene_input = line.strip()
            gene_data = api_get("suggestions", {'gene': gene_input})
            if len(gene_data) > 0:
                canonical_gene = gene_data[0]['canonical']
            else:
                print "No suggestions found for " + gene_input
                continue

            begin_count = api_get("counts", filtered_params(canonical_gene, BEGIN))
            if END == None:
                end_count = {'article_count': 0}
            else:
                end_count = api_get("counts", filtered_params(canonical_gene, END))

            period_count = begin_count["article_count"] - end_count["article_count"]

            print canonical_gene + " has " + str(period_count) + " article(s)"
            if period_count > 0:
                genes_with_articles += 1

                if ONLY_VARIANTS:
                    variants = api_get("variants", filtered_params(canonical_gene, BEGIN))
                    print "Found " + str(variants["variant_count"]) + " variants"

                if (not ONLY_VARIANTS) or variants["variant_count"] > 0:
                    begin_pmids = get_articles(filtered_params(canonical_gene, BEGIN))
                    if END == None:
                        end_pmids = []
                    else:
                        end_pmids = get_articles(filtered_params(canonical_gene, END))

                    period_pmids = [item for item in begin_pmids if item not in end_pmids]

                    gene_info[canonical_gene] = {'pmids': period_pmids}

    # Aggregate all article info, from which other aggregations will be generated
    all_genes = gene_info.keys()
    gene_info, article_info, disease_info, phenotype_info = aggregate_article_info(gene_info, all_genes)

    print('-'*100)

    # Save relevant article data for each unique article across the input gene set to articles.csv file
    articles_file_path = filename + ".articles.csv"
    print("Article info in " + articles_file_path)
    with codecs.open(articles_file_path, 'wb', 'utf-8') as output_file:
        output_file.write(",".join(["PMID", "Link", "Journal", "Title", "Publication Date", "Matched Genes", "Other Genes", "Matched Gene Variants", "Other Gene Variants", "Diseases", "Phenotypes"]) + "\n")
        for pmid, article in article_info.items():
            if skip_article(article, all_genes):
                continue

            url = re.sub(r'([&\?])gene=[^&]+', r'\1gene=' + encode(article["matched_genes"][0]), article["url"])
            pmid_diseases = []
            pmid_phenotypes = []
            pmid_genes = []
            pmid_variants = []
            if "journal" in article:
                journal = article["journal"]
            else:
                journal = "[None]"
            if "title" in article:
                title = article["title"]
            else:
                title = "[None]"
            if "diseases" in article:
                pmid_diseases = [disease["key"] or "[None]" for disease in article["diseases"]]
            if "hpo_terms" in article:
                pmid_phenotypes = [phenotype["term"] or "[None]" for phenotype in article["hpo_terms"]]

            output_file.write(",".join([pmid] + [pipe_delimited_field(field) for field in [[url], [journal], [title], [article["publication_date"]], article["matched_genes"], article["other_genes"], article["matched_gene_variants"], article["other_gene_variants"], pmid_diseases, pmid_phenotypes]]) + "\n")

    print('-'*100)

    # Save relevant article data, organized by unique diseases to diseases.csv
    diseases_file_path = filename + ".diseases.csv"
    print("Disease info in " + diseases_file_path)
    with codecs.open(diseases_file_path, 'wb', 'utf-8') as output_file:
        output_file.write(",".join(["Disease", "Link", "PMIDs", "Matched Genes", "Other Genes", "Matched Gene Variants", "Other Gene Variants", "Other Diseases", "Phenotypes"]) + "\n")
        for disease, info in disease_info.items():
            url = "https://mastermind.genomenon.com/detail?gene=" + encode(next(iter(info["matched_genes"]))) + "&disease=" + encode(disease)
            output_file.write(",".join([pipe_delimited_field(field) for field in [[disease], [url], info["pmids"], info["matched_genes"], info["other_genes"], info["matched_gene_variants"], info["other_gene_variants"], info["diseases"], info["phenotypes"]]]) + "\n")

    print('-'*100)

    # Save relevant article data, organized by unique phenotypes to phenotypes.csv
    phenotypes_file_path = filename + ".phenotypes.csv"
    print("Phenotype info in " + phenotypes_file_path)
    with codecs.open(phenotypes_file_path, 'wb', 'utf-8') as output_file:
        output_file.write(",".join(["Phenotype", "Link", "PMIDs", "Matched Genes", "Other Genes", "Matched Gene Variants", "Other Gene Variants", "Diseases", "Other Phenotypes"]) + "\n")
        for phenotype, info in phenotype_info.items():
            url = "https://mastermind.genomenon.com/detail?gene=" + encode(next(iter(info["matched_genes"]))) + "&hpo=" + encode(info["id"])
            output_file.write(",".join([pipe_delimited_field(field) for field in [[phenotype], [url], info["pmids"], info["matched_genes"], info["other_genes"], info["matched_gene_variants"], info["other_gene_variants"], info["diseases"], info["phenotypes"]]]) + "\n")

    print('-'*100)

    # Save relevant article data, organized by unique genes from input gene set to genes.csv
    genes_file_path = filename + ".genes.csv"
    print("Gene info in " + genes_file_path)
    with codecs.open(genes_file_path, 'wb', 'utf-8') as output_file:
        output_file.write(",".join(["Gene", "Link", "PMIDs", "Matched Genes", "Other Genes", "Gene Variants", "Matched Gene Variants", "Other Gene Variants", "Diseases", "Phenotypes"]) + "\n")
        for gene, info in gene_info.items():
            url = "https://mastermind.genomenon.com/detail?gene=" + encode(gene)
            gene_diseases = [gene_disease for gene_disease in info["diseases"].keys()]
            gene_phenotypes = [gene_phenotype for gene_phenotype in info["phenotypes"].keys()]
            matched_genes = [matched_gene for matched_gene in info["genes"]["matched_genes"].keys()]
            other_genes = [other_gene for other_gene in info["genes"]["other_genes"].keys()]
            gene_variants = [gene_variant for gene_variant in info["variants"]["gene_variants"].keys()]
            matched_gene_variants = [gene_variant for gene_variant in info["variants"]["matched_gene_variants"].keys()]
            other_gene_variants = [gene_variant for gene_variant in info["variants"]["other_gene_variants"].keys()]
            output_file.write(",".join([pipe_delimited_field(field) for field in [[gene], [url], info["filtered_pmids"], matched_genes, other_genes, gene_variants, matched_gene_variants, other_gene_variants, gene_diseases, gene_phenotypes]]) + "\n")

    print('-'*100)

    # Save structured associations lists for input gene set to associations.txt
    associations_file_path = filename + ".associations.txt"
    print("Association info in " + associations_file_path)
    with codecs.open(associations_file_path, 'wb', 'utf-8') as output_file:
        for gene, data in gene_info.items():
            if len(data['filtered_pmids']) == 0:
                qualifier = '' if len(data['pmids']) == 0 else ' and filters'
                output_file.write("No articles found with " + str(gene).upper() + qualifier + "\n")
            else:
                qualifier = '' if len(data['pmids']) == len(data['filtered_pmids']) else ' with filters'
                output_file.write("Found the following associations for " + str(gene).upper() + ":\n")

                output_file.write("\tPMIDs:\n")
                output_file.write("\t\t" + ', '.join(data['filtered_pmids']) + "\n")

                output_file.write("\tGenes with supporting PMIDs:\n")
                if len(data['genes']) == 0:
                    output_file.write("\t\tNone\n")
                else:
                    for values in sorted(data['genes'].items(), key=lambda item: len(item[1]), reverse=True):
                        output_file.write("\t\t" + str(values[0]) + ": " + ', '.join(values[1]) + "\n")

                output_file.write("\tVariants with supporting PMIDs:\n")
                if len(data['variants']) == 0:
                    output_file.write("\t\tNone\n")
                else:
                    for values in sorted(data['variants'].items(), key=lambda item: len(item[1]), reverse=True):
                        output_file.write("\t\t" + str(values[0]) + ": " + ', '.join(values[1]) + "\n")

                output_file.write("\tDiseases with supporting PMIDs:\n")
                if len(data['diseases']) == 0:
                    output_file.write("\t\tNone\n")
                else:
                    for values in sorted(data['diseases'].items(), key=lambda item: len(item[1]), reverse=True):
                        output_file.write("\t\t" + str(values[0]).title() + ": " + ', '.join(values[1]) + "\n")

                output_file.write("\tPhenotypes with supporting PMIDs:\n")
                if len(data['phenotypes']) == 0:
                    output_file.write("\t\tNone\n")
                else:
                    for values in sorted(data['phenotypes'].items(), key=lambda item: len(item[1]), reverse=True):
                        output_file.write("\t\t" + str(values[0]).title() + ": " + ', '.join(values[1]) + "\n")

if __name__ == "__main__":
    main(sys.argv)
