#!/usr/bin/env python3

import argparse, sys
from collections import Counter

# from Readfq by Heng Li, https://github.com/lh3/readfq
# yields reads as tuples of name, sequence, qualities (all strings)
def readfq(fp): # this is a generator function
	last = None # this is a buffer keeping the last unprocessed line
	while True: # mimic closure; is it a bad idea?
		if not last: # the first record or a record following a fastq
			for l in fp: # search for the start of the next record
				if l[0] in '>@': # fasta/q header line
					last = l[:-1] # save this line
					break
		if not last: break
		name, seqs, last = last[1:].partition(" ")[0], [], None
		for l in fp: # read the sequence
			if l[0] in '@+>':
				last = l[:-1]
				break
			seqs.append(l[:-1])
		if not last or last[0] != '+': # this is a fasta record
			yield name, ''.join(seqs), None # yield a fasta record
			if not last: break
		else: # this is a fastq record
			seq, leng, seqs = ''.join(seqs), 0, []
			for l in fp: # read the quality
				seqs.append(l[:-1])
				leng += len(l) - 1
				if leng >= len(seq): # have read enough quality
					last = None
					yield name, seq, ''.join(seqs); # yield a fastq record
					break
			if last: # reach EOF before reading enough quality
				yield name, seq, None # yield a fasta record instead
				break
# end of Readfq

def writefq (name, seq, qualities):
	return '@%s\n%s\n+\n%s\n' % (name, seq, qualities)


# parse arguments
parser = argparse.ArgumentParser(description = 'for 3SEQ-like reads in FASTQ format, trim 3\' homopolymers (e.g. polyA tails) and whatever follows them')
parser.add_argument('-b', '--homopolymer_base', action = 'store', type = str, default = 'A', help = 'base to look for in 3\' homopolymer')
parser.add_argument('-p', '--homopolymer_length', action = 'store', type = int, default = 8, help = 'minimum homopolymer length for detecting a stretch to trim')
parser.add_argument('-m', '--mismatches', action = 'store', type = int, default = 1, help = 'number of allowed mismatches in homopolymer')
parser.add_argument('-l', '--min_length', action = 'store', type = int, default = 1, help = 'minimum length of sequences to keep')
parser.add_argument('-L', '--truncate_length', action = 'store', type = int, default = 0, help = 'length to truncate sequences (before trimming) , 0 for no truncation')
parser.add_argument('infile', action = 'store', nargs = '?', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('outfile', action = 'store', nargs = '?', type = argparse.FileType('w'), default = sys.stdout)
args = parser.parse_args()
assert(len(args.homopolymer_base) == 1)

read_counter = 0
read_length_counter = Counter()
for name, seq, qualities in readfq(args.infile):

	# truncate read
	if args.truncate_length != 0:
		seq = seq[:args.truncate_length]
		qualities = qualities[:args.truncate_length]
	
	# find A stretch
	polya_pos = len(seq)
	if not args.homopolymer_length == 0:
		for pos in range(len(seq)):
			if seq[pos] == args.homopolymer_base and seq[pos:(pos + args.homopolymer_length)].count(args.homopolymer_base) >= args.homopolymer_length - args.mismatches:
				polya_pos = pos
				break
	read_length_counter[polya_pos] += 1
	
	# write trimmed read
	if polya_pos >= args.min_length: # don't write if too short
		seq = seq[:polya_pos]
		qualities = qualities[:polya_pos]
		args.outfile.write(writefq(name, seq, qualities))
	
	read_counter += 1

# print read count
sys.stderr.write('%i reads processed\n\n' % read_counter)

# print insert length frequency
if not args.homopolymer_length == 0:
	sys.stderr.write('trimmed read length counts\n')
	for (length, count) in sorted(read_length_counter.items()):
		sys.stderr.write('%i\t%i\n' % (length, count))

