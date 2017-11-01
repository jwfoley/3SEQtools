#! /usr/bin/env python3

# given a GTF annotation file, return pairs of gene_id, gene_name

import sys, re

line_number = 0
for line in sys.stdin:
	line_number += 1
	if line.startswith('#'): continue
	line_split = line.rstrip().split('\t')
	if line_split[2] != 'gene': continue
	gene_id = gene_name = None
	for match in re.finditer('(\w+) "(.+?)"', line_split[8]):
		if match.group(1) == 'gene_id':
			gene_id = match.group(2)
		elif match.group(1) == 'gene_name':
			gene_name = match.group(2)
	print(gene_id, gene_name)

