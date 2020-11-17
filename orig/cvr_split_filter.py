#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Split the Mastermind Cited Variants Reference (CVR) VCF file into separate
VCF files per variant category, and filter out more complex Multiple Nucleotide
Variants (MNVs) with only protein-level citations in the literature.

Usage:

Run the file with the path to the input CVR VCF file on the local filesystem:
    ./cvr_split_filter.py "/path/to/mastermind_cited_variants_reference-grch37.vcf"

This will produce several separate gzip-compressed VCF output files.
"""

import sys
import re
import gzip
import fnmatch


def open_file(filename):
    if fnmatch.fnmatch(filename, "*.gz"):
        with gzip.open(filename, "r") as lines:
            parse(lines, filename)
    else:
        with open(filename, "r") as lines:
            parse(lines, filename)


def parse(lines, filename):
    readlines = False
    split = {
        "non-coding": [],
        "coding": [],
        "substitution-snvs-cdna": [],
        "substitution-snvs-protein": [],
        "substitution-mnvs-cdna": [],
        "substitution-mnvs-protein": [],
        "other": [],
    }
    split_substitution_mnv_possibilities = {}
    header_lines = []

    i = 0
    for line in lines:
        if readlines:
            i += 1
            # Only process first 1000 lines of VCF file for faster debugging
            # if i > 10000:
            # break
            out = re.split(r"\t", line)
            chrom, pos, ref, alt, data = out[0], out[1], out[3], out[4], out[7]
            if ref == ".":
                ref = ""
            if alt == ".":
                alt = ""
            data = re.split(";", data)
            mmcnt1 = int(re.split("=", data[2])[1])
            mmcnt2 = int(re.split("=", data[3])[1])
            mmcnt3 = int(re.split("=", data[4])[1])
            mmid3 = re.split("=", data[-2])[1].rstrip()
            mmid3_list = re.split(",", mmid3)
            hgvs = re.split("=", data[0])[1].rstrip()
            gene = mmid3.split(":", 1)[0]

            if re.search(r":.*fs(,|$)", mmid3, flags=re.IGNORECASE):
                split["coding"].append(line)
            elif re.search(r":.*dup(,|$)", mmid3, flags=re.IGNORECASE):
                split["coding"].append(line)
            elif re.search(r":.*delins(,|$)", mmid3, flags=re.IGNORECASE):
                split["coding"].append(line)
            elif re.search(r":.*del(,|$)", mmid3, flags=re.IGNORECASE):
                split["coding"].append(line)
            elif re.search(r":.*ins(,|$)", mmid3, flags=re.IGNORECASE):
                split["coding"].append(line)
            elif re.search(r":.*inv(,|$)", mmid3, flags=re.IGNORECASE):
                split["coding"].append(line)
            elif re.search(r":.*ext(,|$)", mmid3, flags=re.IGNORECASE):
                split["coding"].append(line)
            elif re.search(
                r":[ARNDCQEGHILKMFPSTWYVXMU]\d+[ARNDCQEGHILKMFPSTWYVXMU](,|$)",
                mmid3,
                flags=re.IGNORECASE,
            ):
                if len(ref) == 1:
                    if mmcnt1 > 0:
                        split["substitution-snvs-cdna"].append(line)
                    else:
                        split["substitution-snvs-protein"].append(line)
                else:
                    if mmid3 in split_substitution_mnv_possibilities:
                        split_substitution_mnv_possibilities[mmid3].append(
                            [ref, mmcnt1, line]
                        )
                    else:
                        split_substitution_mnv_possibilities[mmid3] = [
                            [ref, mmcnt1, line]
                        ]
            elif re.search(r":.*sa(,|$)", mmid3, flags=re.IGNORECASE):
                split["non-coding"].append(line)
            elif re.search(r":.*sd(,|$)", mmid3, flags=re.IGNORECASE):
                split["non-coding"].append(line)
            elif re.search(r":.*3utr(,|$)", mmid3, flags=re.IGNORECASE):
                split["non-coding"].append(line)
            elif re.search(r":.*5utr(,|$)", mmid3, flags=re.IGNORECASE):
                split["non-coding"].append(line)
            elif re.search(r":.*int(,|$)", mmid3, flags=re.IGNORECASE):
                split["non-coding"].append(line)
            elif re.search(r":.*ugv(,|$)", mmid3, flags=re.IGNORECASE):
                split["non-coding"].append(line)
            elif re.search(r":.*multi-intron", mmid3, flags=re.IGNORECASE):
                split["non-coding"].append(line)
            else:
                split["other"].append(line)

        else:
            header_lines.append(line)

            if line.find("#CHROM") == 0:
                readlines = True

    for mmid3, values in split_substitution_mnv_possibilities.iteritems():
        min_ref = min(map(lambda x: len(x[0]), values))
        for ref, mmcnt1, line in values:
            if mmcnt1 > 0:
                split["substitution-mnvs-cdna"].append(line)
            elif len(ref) == min_ref:
                split["substitution-mnvs-protein"].append(line)

    for cat, lines in split.iteritems():
        if len(lines) > 0:
            print cat + ": " + str(len(lines))
            output_file_path = re.sub(r"\.vcf(\.gz)?$", "." + cat + ".vcf.gz", filename)

            print "\tSaving annotations to " + output_file_path
            with gzip.open(output_file_path, "wb") as output_file:
                output_file.write("".join(header_lines))
                output_file.write("".join(lines))


def main(args):
    open_file(args[-1])


if __name__ == "__main__":
    main(sys.argv)
