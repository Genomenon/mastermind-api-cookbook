#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Using the Mastermind API, produce a CSV-formatted output with evidence
citation counts for variants from input VCF file with focus on specified
disease, with additional evidence citation counts for each variant's gene and
other most commonly associated diseases from the medical literature.

Usage:

Before running this code, first change the API_TOKEN line below to your
Mastermind API key and save the updated file.

Then run the file directly from the command line, replacing the first two
arguments with the path to the input VCF file and the desired disease MeSH
term, respectively:
    ./gene_variant_disease_counts_csv.py "/path/to/input_file.vcf" "melanoma"

If you would like to save the results, just direct the output to the desired
output filename (be careful, this command will overwrite the specified file if
it already exists):
    ./gene_variant_disease_counts_csv.py "/path/to/input_file.vcf" "melanoma" > output_file.csv
"""

import sys
import re

import base

GRCH37_ACCESSION_NUMBERS = {
        '1': 'NC_000001.10',
        '2': 'NC_000002.11',
        '3': 'NC_000003.11',
        '4': 'NC_000004.11',
        '5': 'NC_000005.9',
        '6': 'NC_000006.11',
        '7': 'NC_000007.13',
        '8': 'NC_000008.10',
        '9': 'NC_000009.11',
        '10': 'NC_000010.10',
        '11': 'NC_000011.9',
        '12': 'NC_000012.11',
        '13': 'NC_000013.10',
        '14': 'NC_000014.8',
        '15': 'NC_000015.9',
        '16': 'NC_000016.9',
        '17': 'NC_000017.10',
        '18': 'NC_000018.9',
        '19': 'NC_000019.9',
        '20': 'NC_000020.10',
        '21': 'NC_000021.8',
        '22': 'NC_000022.10',
        'X': 'NC_000023.10',
        'Y': 'NC_000024.9'
        }

def generate_lines(variant, canonical_disease):
    lines = []

    variant_data = base.api_request("suggestions", {'variant': variant})

    if len(variant_data) == 0:
        sys.stderr.write("Variant " + variant + " not found. Skipping.\n")
        sys.stderr.flush()
        return []

    encoded_variant = base.encode(variant_data[0]['matched'])

    for match in variant_data:
        new_line = []
        canonical_variant = match['canonical']
        canonical_gene, just_variant = canonical_variant.split(':', 1)

        new_line.append(canonical_gene)
        new_line.append(just_variant)
        new_line.append(canonical_variant)
        new_line.append(match['url'])

        variant_disease_count_data = base.api_request("counts", {'variant': canonical_variant, 'disease': canonical_disease})
        new_line.append(str(variant_disease_count_data['article_count']))

        variant_count_data = base.api_request("counts", {'variant': canonical_variant})
        new_line.append(str(variant_count_data['article_count']))

        variant_diseases_with_counts = []
        variant_disease_data = base.api_request("diseases", {'variant': canonical_variant})

        if 'diseases' in variant_disease_data:
            for disease in variant_disease_data['diseases']:
                variant_diseases_with_counts.append(str(disease['key']) + "(" + str(disease['article_count']) + ")")

        new_line.append('"' + '|'.join(variant_diseases_with_counts) + '"')

        gene_count_data = base.api_request("counts", {'gene': canonical_gene})
        new_line.append(gene_count_data['url'])

        gene_disease_count_data = base.api_request("counts", {'gene': canonical_gene, 'disease': canonical_disease})
        new_line.append(str(gene_disease_count_data['article_count']))

        new_line.append(str(gene_count_data['article_count']))

        gene_diseases_with_counts = []
        gene_disease_data = base.api_request("diseases", {'gene': canonical_gene})

        if 'diseases' in gene_disease_data:
            for disease in gene_disease_data['diseases']:
                gene_diseases_with_counts.append(str(disease['key']) + "(" + str(disease['article_count']) + ")")

        new_line.append('"' + '|'.join(gene_diseases_with_counts) + '"')

        lines.append(new_line)

    return lines

def main(args):
    filename = args[1]
    disease = args[2]
    output = []

    disease_data = base.api_request("suggestions", {'disease': disease})
    canonical_disease = disease_data[0]['canonical']
    encoded_disease = base.encode(canonical_disease)

    output.append(["SYMBOL", "Variant", "MM Code", "MM Variant Link", "Variant " + "\"" + re.sub(r"\"", "\"\"", canonical_disease) + "\"" + " Articles in MM (Article Count", "Variant Article Count in MM", "Variant Diseases in MM (Article Count)", "MM Gene Link", "Gene " + canonical_disease + " Articles in MM (Article Count)", "Gene Article Count in MM", "Gene Diseases in MM (Article Count)"])

    with open(filename, "r") as lines:
        readlines = False
        for line in lines:
            if readlines:
                out = re.split(r'\t', line)
                chrom, pos, ref, alt, data = out[0], out[1], out[3], out[4], out[7]
                if ref == '.':
                    ref = ''
                if alt == '.':
                    alt = ''

                if chrom in GRCH37_ACCESSION_NUMBERS:
                    accession = GRCH37_ACCESSION_NUMBERS[chrom]
                else:
                    accession = chrom
                variant = accession + ':g.' + pos
                if len(ref) == 1 and len(alt) == 1:
                    variant += ref + '>'
                elif len(ref) == 0:
                    variant += 'ins'
                else:
                    if len(ref) > 1:
                        variant += '_' + str(int(pos) + len(ref)-1)

                    if len(alt) == 0:
                        variant += 'del'
                    else:
                        variant += 'delins'
                variant += alt

                lines = generate_lines(variant, canonical_disease)
                for line in lines:
                    output.append(line)

            else:
                if line.find("#CHROM") == 0:
                    readlines = True

    for line in output:
        print ','.join(line)

if __name__ == "__main__":
    main(sys.argv)
