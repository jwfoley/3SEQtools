#! /usr/bin/env python3

import re, collections, sys, argparse, pysam


class GtfParser:
	'''
	generator that yields dictionaries containing genes' feature data plus lists of exons
	'''

	def __init__ (self, gtf_file, ref_order):
		self.ref_index = dict((ref_order[i], i) for i in range(len(ref_order)))
		self.file = gtf_file
		self.next_feature = None
		self.readline()
		
	def readline (self):
		line = next(self.file).rstrip() # StopIteration is passed all the way up through __next__
		while line.startswith('#'): line = next(self.file).rstrip()
		fields = line.split('\t')
		new_feature = dict([
			('reference_id',  self.ref_index[fields[0]]),
			('feature_type',  fields[2]),
			('left',          int(fields[3])),
			('right',         int(fields[4])),
			('strand',        fields[6])
		])
			
		# quick and dirty parsing of specific other features instead of parsing all of them
		for attribute in ['gene_id', 'gene_type']: new_feature[attribute] = re.search('%s "(.*?)"' % attribute, fields[8]).groups()[0]
		
		assert new_feature['left'] <= new_feature['right'], 'invalid coordinates in GTF: %s' % new_feature['gene_id']
		assert self.next_feature is None or (
			new_feature['reference_id'] > self.next_feature['reference_id'] or
			(new_feature['reference_id'] == self.next_feature['reference_id'] and new_feature['left'] >= self.next_feature['left'])
		), 'features out of order: %s' % new_feature['gene_id']
		
		self.next_feature = new_feature
	
	def __next__ (self):
		if self.next_feature is None: raise StopIteration
		assert self.next_feature['feature_type'] == 'gene'
		gene = self.next_feature
		gene['exons'] = []
		while True:
			try:
				self.readline()
			except StopIteration:
				self.next_feature = None
				break
			if self.next_feature['feature_type'] == 'gene':
				break
			elif self.next_feature['feature_type'] == 'exon':
				for attribute in ['reference_id', 'strand', 'gene_id', 'gene_type']: assert gene[attribute] == self.next_feature[attribute]
				gene['exons'] += [(self.next_feature['left'], self.next_feature['right'])]
		if gene['strand'] == '-': gene['exons'] = gene['exons'][::-1] # reverse-strand genes have exons from TSS to TTS, not left to right
#		assert gene['exons'][0][0] == gene['left'] and gene['exons'][-1][1] == gene['right']
		
		return gene
	
	def __iter__ (self):
		return self


parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-d', '--debug', action = 'store_true')
parser.add_argument('gtf_file', action = 'store', type = argparse.FileType('r'))
parser.add_argument('bam_file', action = 'store', nargs = '?', type = str, default = '-')
args = parser.parse_args()


sam = pysam.Samfile(args.bam_file)
gtf = GtfParser(args.gtf_file, sam.references)


genes = collections.deque()
counts = collections.OrderedDict((category, 0) for category in ['total alignments', 'no annotated gene', 'in exon'])

for alignment in sam:
	unannotated = True
	in_exon = False

	left, right = alignment.reference_start + 1, alignment.reference_end + 1 # pysam is 0-based but GTF is 1-based, so let's agree on 1-based
	
	if args.debug: print('%s\t%s:%i-%i' % (alignment.query_name, alignment.reference_name, left, right), file = sys.stderr)
	if args.debug: print('\tbuffer: %i' % len(genes), file = sys.stderr)
	
	# remove genes that have already been passed
	while genes and	(
		genes[0]['reference_id'] < alignment.reference_id or
		(genes[0]['reference_id'] == alignment.reference_id and genes[0]['right'] < left)
	):
		if args.debug: print('\tpopleft:\t%s\t%s:%i-%i' % (genes[0]['gene_id'], sam.references[genes[0]['reference_id']], genes[0]['left'], genes[0]['right']), file = sys.stderr)
		genes.popleft() # discard genes we have safely passed
	
	# add more genes until they're safely past this alignment
	while not genes or genes[-1]['reference_id'] == alignment.reference_id and genes[-1]['left'] <= right:
		try:
			new_gene = next(gtf)
			assert not genes or (
				next_gene['reference_id'] > genes[-1]['reference_id'] or
				(next_gene['reference_id'] == genes[-1]['reference_id'] and next_gene['left'] >= genes[-1]['left'])
			), 'genes out of order in GTF: %s' % new_gene['gene_id']
			
			if (
				new_gene['reference_id'] > alignment.reference_id or
				(new_gene['reference_id'] == alignment.reference_id and new_gene['left'] > right)
			):
				genes.append(new_gene)
				if args.debug: print('\tappend:\t%s\t%s:%i-%i' % (new_gene['gene_id'], sam.references[new_gene['reference_id']], new_gene['left'], new_gene['right']), file = sys.stderr)
			else: # gene has already been passed
				if args.debug: print('\tskip:\t%s\t%s:%i-%i' % (new_gene['gene_id'], sam.references[new_gene['reference_id']], new_gene['left'], new_gene['right']), file = sys.stderr)
		except StopIteration:
			break

	for gene in genes:
		if (
			gene['reference_id'] > alignment.reference_id or
			(gene['reference_id'] == alignment.reference_id and gene['left'] > right)
		):
			break # stop looking when the next gene is past this alignment
			
		# identify overlap (all possible cases)	
		if (
			(left >= gene['left'] and left <= gene['right']) or
			(right >= gene['left'] and right <= gene['right']) or
			(left < gene['left'] and right > gene['right']) # gene is entirely inside the alignment!			
		):
			unannotated = False
			if args.debug: print('\thit:\t%s\t%s:%i-%i' % (gene['gene_id'], sam.references[gene['reference_id']], gene['left'], gene['right']), file = sys.stderr)
			
			for exon in gene['exons']:
				if left >= exon[0] and right <= exon[1]: # completely inside this exon
					in_exon = True
					break
				elif right < exon[0]: # completely to the left of this exon, so stop looking
					break
	
	counts['total alignments'] += 1
	counts['no annotated gene'] += unannotated
	counts['in exon'] += in_exon


for category, count in counts.items(): print('%s\t%i' % (category, count))

