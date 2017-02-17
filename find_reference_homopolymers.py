#!/usr/bin/env python3

# simplified first pass:
# don't look for complement on opposite strand
# don't allow IUPAC base codes outside the usual ACGT

# rules:
# homopolymer region must be at least args.length bases long
# homopolymer region may not contain two consecutive mismatches


import re, argparse, sys

COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

CHR_HEADER = re.compile('> *(.+)')


parser = argparse.ArgumentParser(description = 'given reference sequence in FASTA format, produce a list in BED format of homopolymers of a given base')
parser.add_argument('-l', '--length', action = 'store', type = int, default = 10, help = 'minimum length of homopolymer')
parser.add_argument('-b', '--base', action = 'store', type = str, default = 'A', help = 'target base')
parser.add_argument('in_fastq', action = 'store', nargs = '?', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('out_bed', action = 'store', nargs = '?', type = argparse.FileType('w'), default = sys.stdout)
args = parser.parse_args()

assert len(args.base) == 1
target_base = args.base.upper()
#try:
#	complement = COMPLEMENTS[target_base]
#except KeyError:
#	raise NotImplementedError(('unknown base %s; options: ' % target_base) + ' '.join(COMPLEMENTS.keys()))


def read_base (fastq): # simple generator that spits out (chr, pos, base) one at a time
	for fastq_line in fastq:
		fastq_line = fastq_line.rstrip()
		if len(fastq_line) == 0: continue
		
		match = re.match(CHR_HEADER, fastq_line)
		
		if match: # new chromosome header
			current_pos = 0
			chrom = match.group(1)
		
		else: # sequence data
			for base in fastq_line:
				current_pos += 1
				yield (chrom, current_pos, base.upper())

def output_region (chrom, start, end): # write a region in BED format
	args.out_bed.write('%s\t%i\t%i\n' % (chrom, start - 1, end))

def process_region (chrom, start, region): # find the end of a region and output if it qualifies; assume the region starts on a target base
	assert(region[0])
	region_length = len(region) - 1 + region[-1] # subtract 1 if the final base is a mismatch (if it is, we can assume the previous base is a match)
	if region_length >= args.length:
		output_region(chrom, start, start + region_length - 1)
		for i in range(region_length): sys.stderr.write(str(int(region[i]))) # test
		sys.stderr.write('\n') # test


current_region = []
current_chrom = None
region_start = None

for chrom, pos, base in read_base(args.in_fastq):
	if chrom != current_chrom: # new chromosome
		if len(current_region) != 0: process_region(current_chrom, region_start, current_region)		
		current_chrom = chrom
		region_start = None
		current_region.clear()
	
	base_bool = (base == target_base)
	
	if len(current_region) == 0:
		if base_bool: # beginning of region
			current_region.append(base_bool)
			region_start = pos
	
	else:
		if not base_bool and not current_region[-1]: # end of region
			process_region(chrom, region_start, current_region[:-1])
			current_region.clear()
			region_start = None
		else: # region continues
			current_region.append(base_bool)

if len(current_region) != 0: process_region(current_chrom, region_start, current_region)

