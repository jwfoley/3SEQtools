#! /usr/bin/env python3

import re, collections, sys, argparse, pysam

DEFAULT_END_DISTANCE = 500 # maximum distance from gene end that a read can start to be counted as 3' end


GenomeFeature = collections.namedtuple('GenomeFeature', ['reference_id', 'feature_type', 'left', 'right', 'is_reverse', 'gene_id', 'transcript_id', 'gene_type', 'children']) # this might be easier with a Python 3.7 data class

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

def get_introns (exons):
	'''
	return a list of introns, given exons, both as GenomeFeature instances
	'''
	introns = []
	for i in range(len(exons) - 1):
		for attribute in ['reference_id', 'is_reverse', 'gene_type', 'gene_id', 'transcript_id']: assert getattr(exons[i], attribute) == getattr(exons[i + 1], attribute), 'non-matching attribute %s for %s' % (attribute, exons[i].transcript_id)
		assert exons[i + 1].left > exons[i].right + 1 # intron must have length of at least 1
		introns += [GenomeFeature(
			reference_id =   exons[i].reference_id,
			feature_type =   'intron',
			left =           exons[i].right + 1,
			right =          exons[i + 1].left - 1,
			is_reverse =     exons[i].is_reverse,
			gene_type =      exons[i].gene_type,
			gene_id =        exons[i].gene_id,
			transcript_id =  exons[i].transcript_id,
			children =       []
		)]
	return introns

def make_transcript (transcript, exons):
	'''
	return a GenomeFeature for a transcript with introns as children
	assumes exons are in right-to-left order for reverse-strand transcripts, and reverses that
	'''
	if transcript.is_reverse: exons = exons[::-1]
	assert exons[0].left == transcript.left and exons[-1].right == transcript.right, 'incomplete exons in GTF: %s' % transcript.transcript_id
	return GenomeFeature(reference_id = transcript.reference_id, feature_type = transcript.feature_type, left = transcript.left, right = transcript.right, is_reverse = transcript.is_reverse, gene_type = transcript.gene_type, gene_id = transcript.gene_id, transcript_id = transcript.transcript_id, children = get_introns(exons))

def assert_attributes(attributes, feature1, feature2):
	'''
	assert two features match in all the given attributes
	'''
	for attribute in attributes: assert getattr(feature1, attribute) == getattr(feature2, attribute), 'non-matching %s in GTF: %s' % (attribute, feature1.gene_id)
		

class GtfParser:
	'''
	generator that yields GenomeFeature instances containing each gene's feature data plus a list of its transcripts, each of which is itself a GenomeFeature instance containing a list of its introns (not exons as you might expect)
	gtf_file must be a file object (or stream of lines), not a file name
	'''

	def __init__ (self, gtf_file, ref_order):
		self.ref_index = dict((ref_order[i], i) for i in range(len(ref_order)))
		self.file = gtf_file
		self.next_feature = None
		self.read_feature()
		
	def read_feature (self): # read another feature from the file and put it into self.next_feature
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
		gene_type = re.search('gene_type "(.*?)"', fields[8]).groups()[0]
		gene_id = re.search('gene_id "(.*?)"', fields[8]).groups()[0]
		try:
			transcript_id = re.search('transcript "(.*?)"', fields[8]).groups()[0]
		except AttributeError:
			transcript_id = None
		
		assert left <= right, 'invalid coordinates in GTF: %s' % gene_id
		assert self.next_feature is None or reference_id >= self.next_feature.reference_id, 'features out of order in GTF: %s' % gene_id # don't check left coordinate because it depends whether this is a child or parent feature and that's hard
		
		self.next_feature = GenomeFeature(reference_id = reference_id, feature_type = feature_type, left = left, right = right, is_reverse = is_reverse, gene_type = gene_type, gene_id = gene_id, transcript_id = transcript_id, children = [])
	
	def __next__ (self):
		if self.next_feature is None: raise StopIteration
		assert self.next_feature.feature_type == 'gene'
		gene = self.next_feature
		transcript = None
		transcripts = []
		exons = []
		while True:
			try:
				self.read_feature()
			except StopIteration:
				self.next_feature = None
				break
			if self.next_feature.feature_type == 'gene':
				assert not feature_starts_before(self.next_feature, gene), 'genes out of order in GTF: %s' % self.next_feature.gene_id
				if transcript is not None: transcripts += [make_transcript(transcript, exons)]
				break
						
			elif self.next_feature.feature_type == 'transcript':
				if transcript is not None: transcripts += [make_transcript(transcript, exons)]
				assert_attributes(['reference_id', 'is_reverse', 'gene_type', 'gene_id'], gene, self.next_feature)
				transcript = self.next_feature
				exons = []
			
			elif self.next_feature.feature_type == 'exon':
				assert transcript is not None, 'missing transcript definition for %s' % gene.gene_id
				assert_attributes(['reference_id', 'is_reverse', 'gene_type', 'gene_id', 'transcript_id'], transcript, self.next_feature)
				exons += [self.next_feature]
		
		return GenomeFeature(reference_id = gene.reference_id, feature_type = gene.feature_type, left = gene.left, right = gene.right, is_reverse = gene.is_reverse, gene_type = gene.gene_type, transcript_id = None, gene_id = gene.gene_id, children = transcripts)
	
	def __iter__ (self):
		return self


parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-d', '--ignoredup', action = 'store_true', help = 'ignore reads marked as duplicate')
parser.add_argument('-e', '--end_distance', action = 'store', type = int, default = DEFAULT_END_DISTANCE, help = 'maximum distance from gene end to be considered 3\' ends, default %i' % DEFAULT_END_DISTANCE)
parser.add_argument('--debug', action = 'store_true')
parser.add_argument('gtf_file', action = 'store', type = argparse.FileType('r'))
parser.add_argument('bam_file', action = 'store', nargs = '?', type = str, default = '-')
args = parser.parse_args()


sam = pysam.Samfile(args.bam_file)
gtf = GtfParser(args.gtf_file, sam.references)


genes = collections.deque()
counts = collections.OrderedDict((category, 0) for category in ['total alignments', 'no annotated gene', 'ribosomal', 'wrong strand', 'intron', '3\' end'])
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

	n_hit_gene = n_sense = n_ribosomal = n_hit_intron = n_hit_end = 0
	
	alignment = GenomeFeature(reference_id = raw_alignment.reference_id, feature_type = 'alignment', left = raw_alignment.reference_start + 1, right = raw_alignment.reference_end + 1, is_reverse = raw_alignment.is_reverse, gene_type = None, gene_id = None, transcript_id = None, children = []) # left, right: pysam is 0-based but GTF is 1-based, so let's agree on 1-based
	
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
			
		# identify overlap
		if raw_alignment.get_overlap(gene.left + 1, gene.right + 1) > 0:
			n_hit_gene += 1
			hit_sense = alignment.is_reverse == gene.is_reverse	
			n_sense += hit_sense
			n_hit_end += hit_sense and (
				(not gene.is_reverse and alignment.left >= gene.right - args.end_distance + 2) or
				(gene.is_reverse and alignment.right <= gene.left + args.end_distance)
			)
			if gene.gene_type == 'rRNA':
				n_ribosomal += 1
				break # not interested in intron hits with rRNA genes
			
			if args.debug: print('\thit (%s):\t%s\t%s\t%i\t%i' % (('sense' if hit_sense else 'antisense'), gene.gene_id, sam.references[gene.reference_id], gene.left, gene.right), file = sys.stderr)
			
			# find intron hit (only counts one per transcript)
			for transcript in gene.children:
				if feature_completely_before(transcript, alignment): continue
				if feature_completely_before(alignment, transcript): break
				for intron in transcript.children:
					if raw_alignment.get_overlap(intron.left + 1, intron.right + 1) > 0:
						n_hit_intron += 1
						if args.debug: print('\t\tintron:\t\t%i\t%i' % intron, file = sys.stderr)
						break
					elif alignment.right < intron[0]: # completely to the left of this intron, so stop looking
						break
	
	# update tallies (fix logic)
	counts['total alignments'] +=   1
	counts['no annotated gene'] +=  n_hit_gene == 0
	counts['wrong strand'] +=       n_hit_gene > 0 and n_sense == 0 # only if there were no hits on correct strand
	counts['ribosomal'] +=          n_ribosomal == n_hit_gene > 0
	counts['intron'] +=             n_hit_intron == n_sense > 0 # only if there were no transcripts that were hit without hitting introns (not sure if foolproof)
	counts['3\' end'] +=						n_hit_end > 0

#assert counts['total alignments'] == sum(count for category, count in counts.items() if not category == 'total alignments')
for category, count in counts.items(): print('%s\t%i' % (category, count))

