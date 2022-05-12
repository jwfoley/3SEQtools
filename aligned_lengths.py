#! /usr/bin/env python3

import pysam
from collections import Counter

length_counts = Counter()

align_file = pysam.AlignmentFile('-')
for alignment in align_file: 
	if not (
		alignment.is_unmapped or
		alignment.is_supplementary or
		alignment.is_secondary or
		alignment.is_qcfail
	): length_counts[alignment.query_alignment_end] += 1

print('length\treads')
for (length, count) in sorted(length_counts.items()):
	print('%i\t%i' % (length, count))

