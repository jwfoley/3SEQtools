#!/usr/bin/env python3

# simplified first pass:
# don't allow IUPAC base codes outside the usual ACGT

# rules:
# homopolymer region must be at least args.length bases long
# homopolymer region may not contain two consecutive mismatches

# warning! output may be slightly out of order because regions on different strands can overlap


import re, argparse, sys

COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

CHR_HEADER = re.compile('> *(.+)')


parser = argparse.ArgumentParser(description = 'given reference sequence in FASTA format, produce a list in BED format of homopolymers of a given base', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-l', '--length', action = 'store', type = int, default = 10, help = 'minimum length of homopolymer')
parser.add_argument('-b', '--base', action = 'store', type = str, default = 'A', choices = sorted(COMPLEMENTS.keys()), help = 'target base')
parser.add_argument('-u', '--upstream', action = 'store', type = int, default = 0, help = 'extend BED features this many bases upstream (left for + strand, right for - strand)')
parser.add_argument('-d', '--downstream', action = 'store', type = int, default = 0, help = 'extend BED features this many bases downstream (right for + strand, left for - strand)')
parser.add_argument('in_fasta', action = 'store', nargs = '?', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('out_bed', action = 'store', nargs = '?', type = argparse.FileType('w'), default = sys.stdout)
args = parser.parse_args()

target_base = args.base.upper()
complement = COMPLEMENTS[target_base]


def read_base (fasta): # simple generator that spits out (chr, pos, base) one at a time
	for fasta_line in fasta:
		fasta_line = fasta_line.rstrip()
		if len(fasta_line) == 0: continue
		
		match = re.match(CHR_HEADER, fasta_line)
		
		if match: # new chromosome header
			current_pos = 0
			chrom = match.group(1)
		
		else: # sequence data
			for base in fasta_line:
				current_pos += 1
				yield (chrom, current_pos, base.upper())

def output_region (chrom, strand, start, end): # write a region in BED format
	if strand == 0:
		left_extend, right_extend = args.upstream, args.downstream
	elif strand == 1:
		left_extend, right_extend = args.downstream, args.upstream
	else:
		raise NotImplementedError
	args.out_bed.write('%s\t%i\t%i\t.\t.\t%s\n' % (
		chrom,
		max(1, start - left_extend) - 1,
		end + right_extend,
		'+-'[strand]
	))

def process_region (chrom, strand, start, length, last_base_match): # find the length of a region and output if it qualifies; assume the region starts on a target base
	length = length - 1 + last_base_match # shorten if the last base is a mismatch (which implies the previous base is not)
	if length >= args.length: output_region(chrom, strand, start, start + length - 1)


# strands processed separately, as region_start[strand] etc. where 0 = forward, 1 = reverse
current_chrom = None
region_start = [None, None]
region_length = [0, 0]
last_base_match = [False, False]

for chrom, pos, base in read_base(args.in_fasta):
	if chrom != current_chrom: # new chromosome
		for strand in (0, 1):
			process_region(current_chrom, strand, region_start[strand], region_length[strand], last_base_match[strand])
			region_start[strand] = None
			region_length[strand] = 0
			last_base_match[strand] = False
		current_chrom = chrom
	
	match = (base == target_base, base == complement)
	
	for strand in (0, 1):		
		if region_length[strand] == 0:
			if match[strand]: # beginning of region
				region_start[strand] = pos
				region_length[strand] = 1
				last_base_match[strand] = True
		
		else:
			if not match[strand] and not last_base_match[strand]: # end of region
				process_region(current_chrom, strand, region_start[strand], region_length[strand], last_base_match[strand])
				region_start[strand] = None
				region_length[strand] = 0
				last_base_match[strand] = False
			else: # region continues
				region_length[strand] += 1
				last_base_match[strand] = match[strand]

for strand in (0, 1):
	process_region(current_chrom, strand, region_start[strand], region_length[strand], last_base_match[strand])

