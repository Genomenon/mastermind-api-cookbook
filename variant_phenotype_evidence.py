#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Find and prioritize the literature evidence for a given set of variants and
phenotypes.

Usage:

Before running this code, first change the API_TOKEN line below to your
Mastermind API key and save the updated file.

This takes an input of variants from a given case, as well as an input of
phenotypes for the same case. It then queries the genomic literature evidence
for all variants, and prioritizes variants that have been found to be
associated with each other the most with the phenotype inputs.

In addition to this prioritization, this will output the complete list of
articles matched to the variants, along with all of the other genes, variants,
phenotypes, and diseases associated with those variants from the literature
evidence.

The input lists of variants and phenotypes are supplied as two separate text
files. The variant list should contain one gene/variant combination per line,
like this:

PKHD1:c.10036T>C
npc1:A927V
acads:p.R330C
lrig2:p.Ile852Phe

The phenotype input file will contain one phenotype name per line, like this:

Ovarian neoplasm
Spatial pattern
Abnormality of the musculature
Visual impairment
Cataplexy

Then run the file directly from the command line, passing the file path to the
variants file and phenotypes file, respectively:
    ./variant_phenotype_evidence.py my_variant_list.txt my_phenotype_list.txt

For my_variant_list.txt as an input, this will generate 3 output CSV files:
    my_variant_list.articles.csv
    my_variant_list.phenotypes.csv
    my_variant_list.variants.csv

And one output Associations summary file:
    my_variant_list.associations-summary.txt
    my_variant_list.associations.txt

The top of the Associations summary file will contain any of the input variants
that were co-cited in articles, along with the PMIDs, input phenotypes, and any
diseases in those articles.

Next, the summary file will contain any of the input phenotypes that were
co-cited in articles, along with the PMIDs, input variants, and any diseases in
those articles.

The Associations file will contain the full list of genes, PMIDs, variants,
diseases, and phenotypes found to be associated with each of the input
variants.

The Articles CSV file will contain each PMID matching the input variants, some
meta-data about the publication, and all genes, variants, diseases, phenotypes
cited in each article.

The Phenotypes CSV file will contain each phenotype from the input phenotype
list found in the articles matching the input variants, and all of the
associated data from the articles from each phenotype.

The Variants CSV file will contain each variant from the input variant list
found in the literature, with all of the associated data from the articles from
each variant.

An optional parameter is provided called FILTER_ON_PHENOTYPES which when set to
True, will only match articles which contained at least one of the input
phenotypes. This will filter all other associations, including genes, variants,
and phenotypes to only those found in the phenotype-filtered set of articles.
It is recommended to leave this set to False to maximize sensitivity.

If testing the script or data output, the STOP_AFTER setting may be set to some
small number (typically 1-5), so that the script will exit and generate the
output data files after that many genes with non-zero results are found.

Another optional parmater is provided called to shrink the associations file
for large OMIT_ONE_PMID_MATCHES_FROM_ASSOCIATIONS_FILE input variant data sets
by truncated all listed associations present in only a single PMID with each of
the input variants. Set this to True for larger data sets to increase
specificity (and decrease sensitivity) of the associaitons file.
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

# Don't include articles that lack PHENOTYPES:
FILTER_ON_PHENOTYPES = False

# For testing out, you may want to stop the script after the first X genes containing data are found.
# If so, set this to the number of genes with data you'd like to stop after.
# Input genes which return no data won't count toward this number:
STOP_AFTER = False #5

# For large data sets, shrink size of associations file by ommitting associations with only one matching PMID
OMIT_ONE_PMID_MATCHES_FROM_ASSOCIATIONS_FILE = False

# If the variant inputs are known to be properly formatted, skip running them all through the Suggestions endpoint
SKIP_VARIANT_SUGGESTION_NORMALIZATION = True

# Only look up and analyze articles which have nucleotide-specific citations for the variant;
# Useful for large variant datasets when higher specificity is desired
# If setting to True, ensure SKIP_VARIANT_SUGGESTION_NORMALIZATION is set to True,
# since that step will convert nucleotide-specific variant names to their protein effects.
ONLY_NUCLEOTIDE_CITATIONS = True

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

def variant_dna_format(variant):
    return re.search(r'c\.\d+', variant) or re.search(r'g\.\d+', variant) or re.search(r'rs\d+', variant) or re.search(r'IVS\d', variant)

def specificity_match(variant_dna_specificity, article):
    return (not ONLY_NUCLEOTIDE_CITATIONS or not variant_dna_specificity) or ('matched_dna' in article and article['matched_dna'] == True)

def get_articles(options):
    data = api_get("articles", options)
    variant_dna_specificity = 'variant' in options and variant_dna_format(options['variant'])

    if data and "articles" in data:
        pmids = [article['pmid'] for article in data['articles'] if specificity_match(variant_dna_specificity, article)]

        articles = int(data['article_count'])
        pages = int(data['pages'])

        print_progress(1, pages, prefix = 'Getting ' + str(articles) + ' articles for ' + str(options['variant']) + ':', suffix = 'Complete', bar_length = 50)

        if pages > 1:
            if specificity_match(variant_dna_specificity, data['articles'][-1]):
                for page in range(2, pages+1):
                    print_progress(page, pages, prefix = 'Getting ' + str(articles) + ' articles for ' + str(options['variant']) + ':', suffix = 'Complete', bar_length = 50)

                    options.update({'page': page})
                    data = api_get("articles", options)
                    pmids = pmids + [article['pmid'] for article in data['articles'] if specificity_match(variant_dna_specificity, article)]
                    if not specificity_match(variant_dna_specificity, data['articles'][-1]):
                        sys.stdout.write('\n')
                        sys.stdout.flush()
                        break;
            else:
                sys.stdout.write('\n')
                sys.stdout.flush()
    else:
        articles = 0
        pmids = []

    return articles, pmids

def has_phenotypes(phenotypes, article_hpo_terms):
    article_hpo_ids = [article_hpo_term['key'] for article_hpo_term in article_hpo_terms]
    intersection = [hpo_id for hpo_id in phenotypes.values() if hpo_id in article_hpo_ids]
    return len(intersection) > 0

def variant_key(pmid_gene, pmid_variant):
    return pmid_gene['symbol'] + ':' + pmid_variant['key']

def aggregate_article_info(variant_info, phenotypes):
    article_info = {}
    filtered_article_info = {}
    phenotype_info = defaultdict(lambda: defaultdict(lambda: set([])))
    comentioned_variants = defaultdict(lambda: defaultdict(lambda: set([])))
    comentioned_phenotypes = defaultdict(lambda: defaultdict(lambda: set([])))

    for variant, values in variant_info.items():
        pmids_by_disease = defaultdict(lambda: [])
        pmids_by_phenotype = defaultdict(lambda: [])
        pmids_by_variant = defaultdict(lambda: [])
        pmids_by_gene = defaultdict(lambda: [])
        matching_phenotypes = set([])
        pmids_matching_phenotypes = []
        diseases_matching_phenotypes = set([])

        current = 0
        total = len(values['pmids'])

        for pmid in values['pmids']:
            current += 1
            pmid_matches_phenotypes = False
            print_progress(current, total, prefix = 'Inspecting PMID info for ' + str(variant) + ':', suffix = 'Complete', bar_length = 50)

            if pmid in article_info:
                data = article_info[pmid]
            else:
                # Get article_info for each PMID
                data = api_get("article_info", {'pmid': pmid})
                article_info[pmid] = data

            if data == None:
                print("COULDN'T GET INFO FOR ARTICLE. SKIPPING ANALYSIS FOR PMID " + pmid)
                continue

            if FILTER_ON_PHENOTYPES and ('hpo_terms' not in data or not has_phenotypes(phenotypes, data['hpo_terms'])):
                # print("Aritcle PMID " + pmid + " does not contain phenotypes. Skipping.")
                continue
            else:
                filtered_article_info[pmid] = data

            pmid_diseases = []
            pmid_genes = []
            pmid_variants = []
            pmid_phenotypes = []

            if 'diseases' in data:
                for disease in data['diseases']:
                    pmids_by_disease[disease['key']].append(pmid)
                    pmid_diseases.append(disease['key'])

            for pmid_gene in data['genes']:
                pmids_by_gene[pmid_gene['symbol']].append(pmid)
                pmid_genes.append(pmid_gene['symbol'])

                if 'variants' in pmid_gene:
                    for pmid_variant in pmid_gene['variants']:
                        pmids_by_variant[variant_key(pmid_gene, pmid_variant)].append(pmid)
                        pmid_variants.append(variant_key(pmid_gene, pmid_variant))

            if 'hpo_terms' in data:
                for phenotype in data['hpo_terms']:
                    pmids_by_phenotype[phenotype['term']].append(pmid)
                    pmid_phenotypes.append(phenotype['term'])

                    if phenotype['term'] in phenotypes.keys():
                        matching_phenotypes.add(phenotype['term'])
                        pmid_matches_phenotypes = True
                        phenotype_info[phenotype['term']]['pmids'].add(pmid)
                        if 'diseases' in data:
                            phenotype_info[phenotype['term']]['diseases'].update(pmid_diseases)
                        # All PMIDS will have at least one gene and variant, since we returned all PMIDs by input variants
                        phenotype_info[phenotype['term']]['genes'].update(pmid_genes)
                        phenotype_info[phenotype['term']]['variants'].update(pmid_variants)

            pmid_comentioned_variants = [input_variant for input_variant in variant_info.keys() if input_variant.upper() in [pmid_variant.upper() for pmid_variant in pmid_variants]]
            pmid_comentioned_variants.sort()
            pmid_comentioned_phenotypes = [input_phenotype_term for input_phenotype_term in phenotypes.keys() if input_phenotype_term in pmid_phenotypes]
            pmid_comentioned_phenotypes.sort()

            if len(pmid_comentioned_variants) > 1:
                comentioned_key = "; ".join(pmid_comentioned_variants)
                comentioned_variants[comentioned_key]['pmids'].add(pmid)
                comentioned_variants[comentioned_key]['phenotypes'].update(pmid_comentioned_phenotypes)
                comentioned_variants[comentioned_key]['diseases'].update(pmid_diseases)

            if len(pmid_comentioned_phenotypes) > 1:
                comentioned_key = "; ".join(pmid_comentioned_phenotypes)
                comentioned_phenotypes[comentioned_key]['pmids'].add(pmid)
                comentioned_phenotypes[comentioned_key]['variants'].update(pmid_comentioned_variants)
                comentioned_phenotypes[comentioned_key]['diseases'].update(pmid_diseases)

            if pmid_matches_phenotypes:
                matching_phenotypes.update(pmid_comentioned_phenotypes)
                pmids_matching_phenotypes.append(pmid)
                diseases_matching_phenotypes.update(pmid_diseases)

        # Else leave pmids the same
        if FILTER_ON_PHENOTYPES:
            variant_info[variant]['pmids'] = pmids_matching_phenotypes
        variant_info[variant]['diseases'] = pmids_by_disease
        variant_info[variant]['phenotypes'] = pmids_by_phenotype
        variant_info[variant]['genes'] = pmids_by_gene
        variant_info[variant]['variants'] = pmids_by_variant
        variant_info[variant]['matching_phenotypes'] = matching_phenotypes
        variant_info[variant]['pmids_matching_phenotypes'] = pmids_matching_phenotypes
        variant_info[variant]['diseases_matching_phenotypes'] = diseases_matching_phenotypes

    return variant_info, filtered_article_info, phenotype_info, comentioned_variants, comentioned_phenotypes

def pipe_delimited_field(values):
    return "\"" + re.sub(r"\"", "\"\"", "|".join(values)) + "\""

def main(args):
    print("Welcome to the Variant Phenotype Evidence program powered by Mastermind.")

    if ONLY_NUCLEOTIDE_CITATIONS and not SKIP_VARIANT_SUGGESTION_NORMALIZATION:
        print("")
        print("Misconfiguration: ONLY_NUCLEOTIDE_CITATIONS is set to True, so SKIP_VARIANT_SUGGESTION_NORMALIZATION must also be set to True.")
        print("Exiting.")
        sys.exit(0)

    variants_filename = args[1]
    phenotypes_filename = args[2]
    variants_with_articles = 0
    variant_info = {}
    phenotypes = {}
    phenotype_inputs = []

    with open(phenotypes_filename, "r") as lines:
        phenotype_inputs = lines.read().splitlines()

    phenotypes_parsed = 1
    total_phenotypes = len(phenotype_inputs)
    for phenotype in phenotype_inputs:
        print_progress(phenotypes_parsed, total_phenotypes, prefix = 'Parsing phenotypes:', suffix = 'Complete', bar_length = 50)
        phenotype_data = api_get("suggestions", {'hpo': phenotype})
        if len(phenotype_data) > 0:
            phenotypes[phenotype_data[0]['name']] = phenotype_data[0]['canonical']
        else:
            print("\nPhenotype " + phenotype + " could not be found. Skipping.")
        phenotypes_parsed += 1

    with open(variants_filename, "r") as lines:
        # Loop through lines in variant input file
        for line in lines:
            if STOP_AFTER and variants_with_articles > STOP_AFTER:
                break
            variant_input = line.strip()
            if SKIP_VARIANT_SUGGESTION_NORMALIZATION:
                canonical_variant = variant_input
            else:
                variant_data = api_get("suggestions", {'variant': variant_input})
                if len(variant_data) > 0:
                    canonical_variant = variant_data[0]['canonical']
                else:
                    # Try again with gene suggestion
                    gene_input = variant_input.split(":", 1)
                    gene_data = api_get("suggestions", {'gene': gene_input[0]})
                    if len(gene_data) > 0:
                        canonical_gene = gene_data[0]['canonical']
                        variant_data = api_get("suggestions", {'variant': canonical_gene + ':' + gene_input[1]})
                        if len(variant_data) > 0:
                            canonical_variant = variant_data[0]['canonical']
                        else:
                            print variant_input + " could not be matched to a valid transcript."
                            continue
                    else:
                        print gene_input[0] + " could not be matched to a valid transcript."
                        continue

            if canonical_variant in variant_info:
                print("Articles already fetched for " + canonical_variant)
            else:
                count, articles = get_articles({'variant': canonical_variant})
                variant_info[canonical_variant] = {'pmids': articles}

                if count > 0:
                    variants_with_articles += 1
                else:
                    print("No articles found for " + canonical_variant + " (" + variant_input + ").")

    if variants_with_articles == 0:
        print("Mastermind searched 30 million abstracts and 7 million genomic full-text articles, and no articles cite these variants."
        return

    # Aggregate all article info, from which other aggregations will be generated
    variant_info, article_info, phenotype_info, comentioned_variants, comentioned_phenotypes = aggregate_article_info(variant_info, phenotypes)

    print('-'*100)

    # Save relevant article data for each unique article across the input gene set to articles.csv file
    articles_file_path = variants_filename + ".articles.csv"
    print("Article info in " + articles_file_path)
    with codecs.open(articles_file_path, 'wb', 'utf-8') as output_file:
        output_file.write(",".join(["PMID", "Journal", "Title", "Publication Date", "Genes", "Variants", "Diseases", "Phenotypes"]) + "\n")
        for pmid, article in article_info.items():
            pmid_diseases = []
            pmid_phenotypes = []
            pmid_genes = []
            pmid_variants = []
            if "diseases" in article:
                pmid_diseases = [disease["key"] or "[None]" for disease in article["diseases"]]
            if "hpo_terms" in article:
                pmid_phenotypes = [phenotype["term"] or "[None]" for phenotype in article["hpo_terms"]]
            for gene in article["genes"]:
                pmid_genes.append(gene["symbol"])
                if "variants" in gene:
                    pmid_variants.extend([gene["symbol"] + ":" + variant["key"] for variant in gene["variants"]])

            output_file.write(",".join([pipe_delimited_field(field) for field in [[pmid], [article["journal"]], [article["title"]], [article["publication_date"]], pmid_genes, pmid_variants, pmid_diseases, pmid_phenotypes]]) + "\n")

    print('-'*100)

    # Save relevant article data, organized by unique phenotypes to phenotypes.csv
    phenotypes_file_path = variants_filename + ".phenotypes.csv"
    print("Phenotype info in " + phenotypes_file_path)
    with codecs.open(phenotypes_file_path, 'wb', 'utf-8') as output_file:
        output_file.write(",".join(["Phenotype", "PMIDs", "Genes", "Variants", "Diseases", "Other Phenotypes"]) + "\n")
        for phenotype, info in phenotype_info.items():
            output_file.write(",".join([pipe_delimited_field(field) for field in [[phenotype], info["pmids"], info["genes"], info["variants"], info["diseases"], info["phenotypes"]]]) + "\n")

    print('-'*100)

    # Save relevant article data, organized by unique variants from input variant set to variants.csv
    variants_file_path = variants_filename + ".variants.csv"
    print("Variant info in " + variants_file_path)
    with codecs.open(variants_file_path, 'wb', 'utf-8') as output_file:
        output_file.write(",".join(["Variant", "PMIDs", "Genes", "Other Variants", "Diseases", "Phenotypes"]) + "\n")
        for variant, info in variant_info.items():
            variant_diseases = [variant_disease for variant_disease in info["diseases"].keys()]
            variant_phenotypes = [variant_phenotype for variant_phenotype in info["phenotypes"].keys()]
            variant_genes = [variant_gene for variant_gene in info["genes"].keys()]
            variant_variants = [variant_variant for variant_variant in info["variants"].keys()]
            output_file.write(",".join([pipe_delimited_field(field) for field in [[variant], info["pmids"], variant_genes, variant_variants, variant_diseases, variant_phenotypes]]) + "\n")

    print('-'*100)

    # Save structured associations lists for input gene set to associations-summary.txt
    associations_summary_file_path = variants_filename + ".associations-summary.txt"
    print("Association summary info in " + associations_summary_file_path)
    with codecs.open(associations_summary_file_path, 'wb', 'utf-8') as output_file:

        if len(comentioned_variants) > 0:
            output_file.write("Found co-cited variants in the literature:\n")
            for variants, data in sorted(comentioned_variants.items(), key=lambda item: 1 if len(item[1]['phenotypes']) == 0 else 1/len(item[1]['phenotypes'])):
                output_file.write("\t" + variants + ":\n")
                output_file.write("\t\tPMIDS\n:")
                output_file.write("\t\t\t" + ", ".join(data["pmids"]) + "\n")
                output_file.write("\t\tMatching input phenotypes:\n")
                output_file.write("\t\t\t" + "; ".join(data["phenotypes"]) + "\n")
                output_file.write("\t\tMatching diseases:\n")
                output_file.write("\t\t\t" + "; ".join(data["diseases"]) + "\n")
        else:
            output_file.write("Found no co-cited variants in the literature.\n")

        if len(comentioned_phenotypes) > 0:
            output_file.write("Found co-cited phenotypes in the literature:\n")
            for phenotypes, data in sorted(comentioned_phenotypes.items(), key=lambda item: 1/len(item[1]['variants'])):
                output_file.write("\t" + phenotypes + ":\n")
                output_file.write("\t\tPMIDS:\n")
                output_file.write("\t\t\t" + ", ".join(data["pmids"]) + "\n")
                output_file.write("\t\tMatching input variants:\n")
                output_file.write("\t\t\t" + "; ".join(data["variants"]) + "\n")
                output_file.write("\t\tMatching diseases:\n")
                output_file.write("\t\t\t" + "; ".join(data["diseases"]) + "\n")
        else:
            output_file.write("Found no co-cited phenotypes in the literature.\n")

    # Save structured associations lists for input gene set to associations.txt
    associations_file_path = variants_filename + ".associations.txt"
    print("Association info in " + associations_file_path)
    with codecs.open(associations_file_path, 'wb', 'utf-8') as output_file:

        for variant, data in variant_info.items():
            if len(data['pmids']) == 0:
                output_file.write("No articles found with " + str(variant) + "\n")
            else:
                output_file.write("Found the following associations for " + str(variant) + ":\n")

                output_file.write("\tPMIDs:\n")
                output_file.write("\t\t" + ', '.join(data['pmids']) + "\n")

                output_file.write("\tPhenotypes matched:\n")
                output_file.write("\t\t" + ", ".join(data["matching_phenotypes"]) + "\n")

                output_file.write("\tPMIDs with matching phenotypes:\n")
                output_file.write("\t\t" + ", ".join(data["pmids_matching_phenotypes"]) + "\n")

                output_file.write("\tDiseases from PMIDs with matching phenotypes:\n")
                output_file.write("\t\t" + ", ".join(data["diseases_matching_phenotypes"]) + "\n")

                output_file.write("\tGenes with supporting PMIDs:\n")
                if len(data['genes']) == 0:
                    output_file.write("\t\tNone\n")
                else:
                    for values in sorted(data['genes'].items(), key=lambda item: len(item[1]), reverse=True):
                        if OMIT_ONE_PMID_MATCHES_FROM_ASSOCIATIONS_FILE and len(values[1]) <= 1:
                            output_file.write("\t\t(Truncated associations with only one PMID)\n")
                            break
                        output_file.write("\t\t" + str(values[0]) + ": " + ', '.join(values[1]) + "\n")

                output_file.write("\tVariants with supporting PMIDs:\n")
                if len(data['variants']) == 0:
                    output_file.write("\t\tNone\n")
                else:
                    for values in sorted(data['variants'].items(), key=lambda item: len(item[1]), reverse=True):
                        if OMIT_ONE_PMID_MATCHES_FROM_ASSOCIATIONS_FILE and len(values[1]) <= 1:
                            output_file.write("\t\t(Truncated associations with only one PMID)\n")
                            break
                        output_file.write("\t\t" + str(values[0]) + ": " + ', '.join(values[1]) + "\n")

                output_file.write("\tDiseases with supporting PMIDs:\n")
                if len(data['diseases']) == 0:
                    output_file.write("\t\tNone\n")
                else:
                    for values in sorted(data['diseases'].items(), key=lambda item: len(item[1]), reverse=True):
                        if OMIT_ONE_PMID_MATCHES_FROM_ASSOCIATIONS_FILE and len(values[1]) <= 1:
                            output_file.write("\t\t(Truncated associations with only one PMID)\n")
                            break
                        output_file.write("\t\t" + str(values[0]).title() + ": " + ', '.join(values[1]) + "\n")

                output_file.write("\tPhenotypes with supporting PMIDs:\n")
                if len(data['phenotypes']) == 0:
                    output_file.write("\t\tNone\n")
                else:
                    for values in sorted(data['phenotypes'].items(), key=lambda item: len(item[1]), reverse=True):
                        if OMIT_ONE_PMID_MATCHES_FROM_ASSOCIATIONS_FILE and len(values[1]) <= 1:
                            output_file.write("\t\t(Truncated associations with only one PMID)\n")
                            break
                        output_file.write("\t\t" + str(values[0]).title() + ": " + ', '.join(values[1]) + "\n")

if __name__ == "__main__":
    main(sys.argv)
