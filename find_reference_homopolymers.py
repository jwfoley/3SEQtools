#!/usr/bin/env python3

# simplified first pass:
# don't look for complement on opposite strand
# don't allow mismatches (!!!)
# don't allow IUPAC base codes outside the usual ACGT

# rules:
# homopolymer region must be at least args.length bases long
# at least args.proportion of bases in homopolymer region must be the target base
# the first base of the homopolymer region must be the target base
# the last base of the homopolymer region must be the target base

from __future__ import division
import collections, re, argparse, sys

COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

CHR_HEADER = re.compile('> *(.+)')


parser = argparse.ArgumentParser(description = 'given reference sequence in FASTA format, produce a list in BED format of homopolymers of a given base')
parser.add_argument('-l', '--length', action = 'store', type = int, default = 10, help = 'minimum length of homopolymer')
parser.add_argument('-p', '--proportion', action = 'store', type = float, default = 0.8, help = 'minimum proportion of target base in homopolymer stretch')
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

def process_region (chrom, start, window): # find the end of a region and output if it qualifies; assume the given window starts at the start of the region and ends at least one base past the end of it; returns the length of the found region
	for region_length in range(len(window) - 1, args.length - 1, -1):
		if window[region_length -1]:
			output_region(chrom, start, start + region_length - 1)
			for i in range(region_length): sys.stderr.write(str(int(window[i]))) # test
			sys.stderr.write('\n') # test
			return region_length
	return 0

window = collections.deque()
current_chrom = None
target_start = None
n_target = 0

for chrom, pos, base in read_base(args.in_fastq):
	if chrom != current_chrom: # new chromosome
		if target_start is not None: process_region(current_chrom, target_start, window)
		current_chrom = chrom
		target_start = None
		n_target = 0
		window.clear()
	
	base_bool = (base == target_base)
	window.append(base_bool)
	n_target += base_bool
	
	if target_start is None:
		if len(window) >= args.length:
			# determine whether this is the beginning of a target region
			if n_target / len(window) >= args.proportion and window[0]:
				target_start = pos - len(window) + 1
			else:
				n_target -= window.popleft()
		
	else:
		# determine whether this is the end of a target region
		if n_target / len(window) < args.proportion:
			for i in range(process_region(current_chrom, target_start, window)): n_target -= window.popleft() # clear the current region from the window but leave the remainder
			target_start = None
if target_start is not None: process_region(current_chrom, target_start, window)

