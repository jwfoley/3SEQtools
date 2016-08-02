#!/usr/bin/env python3

import argparse, sys
from Bio import SeqIO
from collections import Counter

# parse arguments
parser = argparse.ArgumentParser(description = 'for 3SEQ-like reads in FASTQ format, trim 3\' homopolymers (e.g. polyA tails) and whatever follows them')
parser.add_argument('-b', '--homopolymer_base', action = 'store', type = str, default = 'A', help = 'base to look for in 3\' homopolymer')
parser.add_argument('-p', '--homopolymer_length', action = 'store', type = int, default = 8, help = 'minimum homopolymer length for detecting a stretch to trim')
parser.add_argument('-m', '--mismatches', action = 'store', type = int, default = 1, help = 'number of allowed mismatches in homopolymer')
parser.add_argument('-l', '--min_length', action = 'store', type = int, default = 1, help = 'minimum length of sequences to keep')
parser.add_argument('-L', '--max_length', action = 'store', type = int, default = 0, help = 'maximum length of sequences to keep (truncate if necessary), 0 for no limit')
parser.add_argument('infile', action = 'store', nargs = '?', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('outfile', action = 'store', nargs = '?', type = argparse.FileType('w'), default = sys.stdout)
args = parser.parse_args()
assert(len(args.homopolymer_base) == 1)

read_counter = 0
read_length_counter = Counter()
for read in SeqIO.parse(args.infile, 'fastq'):
	
	# find A stretch
	polya_pos = len(read.seq)
	if not args.homopolymer_length == 0:
		for pos in range(len(read.seq)):
			if read.seq[pos] == args.homopolymer_base and read.seq[pos:(pos + args.homopolymer_length)].count(args.homopolymer_base) >= args.homopolymer_length - args.mismatches:
				polya_pos = pos
				break
	read_length_counter[polya_pos] += 1
	
	# write trimmed read
	if polya_pos >= args.min_length: # don't write if too short
		qualities = read.letter_annotations['phred_quality']
		read.letter_annotations = {} # letter annotations must be emptied before changing sequence
		trim_pos = min(polya_pos, args.max_length) if args.max_length else polya_pos
		read.seq = read.seq[:trim_pos]
		read.letter_annotations['phred_quality'] = qualities[:trim_pos]
		args.outfile.write(read.format('fastq'))
	
	read_counter += 1
args.infile.close()
args.outfile.close()

# print read count
sys.stderr.write('%i reads processed\n\n' % read_counter)

# print insert length frequency
if not args.homopolymer_length == 0:
	sys.stderr.write('trimmed read length counts\n')
	for (length, count) in sorted(read_length_counter.items()):
		sys.stderr.write('%i\t%i\n' % (length, count))

