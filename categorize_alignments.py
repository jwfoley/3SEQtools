#! /usr/bin/env python3

import re, collections, sys, argparse, pysam


GenomeFeature = collections.namedtuple('GenomeFeature', ['reference_id', 'feature_type', 'left', 'right', 'is_reverse', 'gene_id', 'gene_type', 'segments']) # this might be easier with a Python 3.7 data class

def feature_starts_before (feature1, feature2):
	'''
	test whether the first GenomeFeature starts before the second GenomeFeature starts, in sorting order by reference_id and left position
	'''
	return (
		feature1.reference_id < feature2.reference_id or
		(feature1.reference_id == feature2.reference_id and feature1.left < feature2.left)
	)

def feature_completely_before (feature1, feature2):
	'''
	test whether the first GenomeFeature ends before the second GenomeFeature starts
	'''
	return (
		feature1.reference_id < feature2.reference_id or
		(feature1.reference_id == feature2.reference_id and feature1.right < feature2.left)
	)


class GtfParser:
	'''
	generator that yields GenomeFeature instances containing genes' feature data plus lists of exons
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
		
		# parse the basic fields
		reference_id =  self.ref_index[fields[0]]
		feature_type =  fields[2]
		left =          int(fields[3])
		right =         int(fields[4])
		if fields[6] == '+':
			is_reverse = False
		elif fields[6] == '-':
			is_reverse = True
		else:
			raise RuntimeError('bad strand')
			
		# quick and dirty parsing of specific other features instead of parsing all of them
		gene_id = re.search('gene_id "(.*?)"', fields[8]).groups()[0]
		gene_type = re.search('gene_type "(.*?)"', fields[8]).groups()[0]
		
		assert left <= right, 'invalid coordinates in GTF: %s' % gene_id
		assert self.next_feature is None or reference_id >= self.next_feature.reference_id, 'features out of order: %s' % gene_id # don't check left coordinate because it depends whether this is a child or parent feature and that's hard
		
		self.next_feature = GenomeFeature(reference_id = reference_id, feature_type = feature_type, left = left, right = right, is_reverse = is_reverse, gene_id = gene_id, gene_type = gene_type, segments = [])
	
	def __next__ (self):
		if self.next_feature is None: raise StopIteration
		assert self.next_feature.feature_type == 'gene'
		gene = self.next_feature
		exons = []
		while True:
			try:
				self.readline()
			except StopIteration:
				self.next_feature = None
				break
			if self.next_feature.feature_type == 'gene':
				break
			elif self.next_feature.feature_type == 'exon':
				for attribute in ['reference_id', 'is_reverse', 'gene_id', 'gene_type']: assert getattr(gene, attribute) == getattr(self.next_feature, attribute), 'non-matching attribute %s for %s' % (attribute, self.next_feature.gene_id)
				exons += [(self.next_feature.left, self.next_feature.right)]
		if gene.is_reverse: exons = exons[::-1] # reverse-strand genes have exons from TSS to TTS, not left to right
#		assert gene.segments[0][0] == gene.left and gene.segments[-1][1] == gene.right
		
		return GenomeFeature(reference_id = gene.reference_id, feature_type = gene.feature_type, left = gene.left, right = gene.right, is_reverse = gene.is_reverse, gene_id = gene.gene_id, gene_type = gene.gene_type, segments = exons)
	
	def __iter__ (self):
		return self


parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-d', '--ignoredup', action = 'store_true', help = 'ignore reads marked as duplicate')
parser.add_argument('--debug', action = 'store_true')
parser.add_argument('gtf_file', action = 'store', type = argparse.FileType('r'))
parser.add_argument('bam_file', action = 'store', nargs = '?', type = str, default = '-')
args = parser.parse_args()


sam = pysam.Samfile(args.bam_file)
gtf = GtfParser(args.gtf_file, sam.references)


genes = collections.deque()
counts = collections.OrderedDict((category, 0) for category in ['total alignments', 'no annotated gene', 'in exon'])
gene_hit_counter = collections.Counter()


previous_alignment = None # saved only to verify sorting
for raw_alignment in sam:
	if (
		raw_alignment.is_qcfail or
		raw_alignment.is_secondary or
		raw_alignment.is_supplementary or
		raw_alignment.is_unmapped or
		(args.ignoredup and raw_alignment.is_duplicate)
	): continue

	hit_gene = hit_sense = hit_exon = hit_intron = False
	
	alignment = GenomeFeature(reference_id = raw_alignment.reference_id, feature_type = 'alignment', left = raw_alignment.reference_start + 1, right = raw_alignment.reference_end + 1, is_reverse = raw_alignment.is_reverse, gene_id = None, gene_type = None, segments = []) # left, right: pysam is 0-based but GTF is 1-based, so let's agree on 1-based
	
	assert previous_alignment is None or (not feature_starts_before(alignment, previous_alignment)), 'alignments out of order: %s' % raw_alignment.query_name
	previous_alignment = alignment
	
	if args.debug: print('%s\t\t\t%s\t%i\t%i' % (raw_alignment.query_name, raw_alignment.reference_name, alignment.left, alignment.right), file = sys.stderr)
	if args.debug: print('\tbuffer: %i' % len(genes), file = sys.stderr)
	
	# remove genes that have already been passed
	while genes and feature_completely_before(genes[0], alignment):
		if args.debug: print('\tpop:\t%s\t%s\t%i\t%i' % (genes[0].gene_id, sam.references[genes[0].reference_id], genes[0].left, genes[0].right), file = sys.stderr)
		genes.popleft() # discard genes we have safely passed
	
	# add more genes until they're safely past this alignment
	while (not genes) or (not feature_completely_before(alignment, genes[-1])):
		try:
			new_gene = next(gtf)
			assert (not genes) or (not feature_starts_before(new_gene, genes[-1])), 'genes out of order in GTF: %s' % new_gene.gene_id
			if not feature_completely_before(new_gene, alignment):
				genes.append(new_gene)
				if args.debug: print('\tappend:\t%s\t%s\t%i\t%i' % (new_gene.gene_id, sam.references[new_gene.reference_id], new_gene.left, new_gene.right), file = sys.stderr)
			else: # gene has already been passed
				if args.debug: print('\tskip:\t%s\t%s\t%i\t%i' % (new_gene.gene_id, sam.references[new_gene.reference_id], new_gene.left, new_gene.right), file = sys.stderr)
		except StopIteration:
			break
	
	# search for gene hits
	for gene in genes:
		if feature_completely_before(alignment, gene): break # stop looking when the next gene is past this alignment
			
		# identify overlap (all possible cases)	
		if (
			(alignment.left >= gene.left and alignment.left <= gene.right) or
			(alignment.right >= gene.left and alignment.right <= gene.right) or
			(alignment.left < gene.left and alignment.right > gene.right) # gene is entirely inside the alignment!			
		):
			hit_gene = True
			if args.debug: print('\thit:\t%s\t%s\t%i\t%i' % (gene.gene_id, sam.references[gene.reference_id], gene.left, gene.right), file = sys.stderr)
			
			for exon in gene.segments:
				if alignment.left >= exon[0] and alignment.right <= exon[1]: # completely inside this exon
					hit_exon = True
					if args.debug: print('\t\texon:\t\t%i\t%i' % exon, file = sys.stderr)
					break
				elif alignment.right < exon[0]: # completely to the left of this exon, so stop looking
					break
	
	counts['total alignments'] += 1
	counts['no annotated gene'] += not hit_gene
	counts['in exon'] += hit_exon


for category, count in counts.items(): print('%s\t%i' % (category, count))

