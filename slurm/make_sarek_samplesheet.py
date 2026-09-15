#!/usr/bin/env python3
"""
Scans a directory of paired FASTQ files and writes a sarek-compatible
samplesheet.csv (patient,sample,lane,fastq_1,fastq_2).

Handles three common naming patterns:
  <sample>_S<N>_L<LLL>_R1_001.fastq.gz / _R2_001.fastq.gz   (bcl2fastq/Illumina)
  <sample>_R1.fastq.gz / _R2.fastq.gz
  <sample>_1.fastq.gz  / _2.fastq.gz

Each distinct sample becomes its own patient (one patient == one sample):
fine for unrelated germline exomes; if you have trios/families, edit the
`patient` column afterwards to group them (same patient ID, one row per
family member as a distinct `sample`).

Usage:
  python3 make_sarek_samplesheet.py /path/to/fastq/dir > samplesheet.csv
  python3 make_sarek_samplesheet.py /path/to/fastq/dir -o samplesheet.csv
"""
import argparse
import csv
import re
import sys
from pathlib import Path

PATTERNS = [
    # sample, lane group, read group
    re.compile(r'^(?P<sample>.+?)_S\d+_L(?P<lane>\d+)_R(?P<read>[12])_001\.fastq\.gz$'),
    re.compile(r'^(?P<sample>.+?)_R(?P<read>[12])\.fastq\.gz$'),
    re.compile(r'^(?P<sample>.+?)_(?P<read>[12])\.fastq\.gz$'),
]


def match_fastq(name):
    for pat in PATTERNS:
        m = pat.match(name)
        if m:
            gd = m.groupdict()
            lane = gd.get('lane') or '1'
            return gd['sample'], lane.lstrip('0') or '0', gd['read']
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_dir', type=Path)
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()

    pairs = {}  # (sample, lane) -> {'1': path, '2': path}
    unmatched = []

    for f in sorted(args.fastq_dir.glob('*.fastq.gz')):
        m = match_fastq(f.name)
        if not m:
            unmatched.append(f.name)
            continue
        sample, lane, read = m
        pairs.setdefault((sample, lane), {})[read] = str(f.resolve())

    writer = csv.writer(args.output)
    writer.writerow(['patient', 'sample', 'lane', 'fastq_1', 'fastq_2'])

    incomplete = []
    for (sample, lane), reads in sorted(pairs.items()):
        if '1' not in reads or '2' not in reads:
            incomplete.append((sample, lane, reads))
            continue
        writer.writerow([sample, sample, lane, reads['1'], reads['2']])

    if unmatched:
        print(f"WARNING: {len(unmatched)} file(s) didn't match any known naming "
              f"pattern, skipped: {unmatched}", file=sys.stderr)
    if incomplete:
        print(f"WARNING: {len(incomplete)} sample/lane group(s) missing an R1 or R2 "
              f"mate, skipped: {incomplete}", file=sys.stderr)
    if not pairs:
        print("ERROR: no FASTQ pairs recognized in this directory. "
              "Check filenames against the patterns in this script's docstring.",
              file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
