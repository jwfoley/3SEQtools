#!/usr/bin/env pypy3

from __future__ import print_function
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
parser = argparse.ArgumentParser(description = 'for Smart-3SEQ reads in FASTQ format starting with a UMI and then a G-overhang (e.g. 5\'-NNNNNGGG-3\'), cut off both, adding the UMI to the read name and discarding the G\'s', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-u', '--umi', action = 'store', type = int, default = 5, help = 'length of UMI')
parser.add_argument('-g', '--g_overhang', action = 'store', type = int, default = 3, help = 'length of G-overhang')
parser.add_argument('-l', '--min_length', action = 'store', type = int, default = 1, help = 'minimum length of sequences to keep (after trimming)')
parser.add_argument('-L', '--truncate_length', action = 'store', type = int, default = 0, help = 'length to truncate sequences (before trimming), 0 for no truncation')
parser.add_argument('--bcl2fastq_trim_length', action = 'store', type = int, default = 35, help = 'length of the all-N reads generated when bcl2fastq trims below the threshold (--minimum-trimmed-read-length)')
parser.add_argument('infile', action = 'store', nargs = '?', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('outfile', action = 'store', nargs = '?', type = argparse.FileType('w'), default = sys.stdout)
args = parser.parse_args()

header_length = args.umi + args.g_overhang

read_counter = 0
overtrimmed_counter = 0
umi_counter = [Counter() for i in range(args.umi)]
read_length_counter = Counter()
for name, seq, qualities in readfq(args.infile):
	read_counter += 1

	# check for overtrimmed read
	if len(seq) == args.bcl2fastq_trim_length == seq.count('N'):
		overtrimmed_counter += 1
		continue
	
	# truncate read
	if args.truncate_length != 0:
		seq = seq[:args.truncate_length]
		qualities = qualities[:args.truncate_length]
	
	read_length_counter[len(seq)- header_length] += 1
	
	if len(seq) - header_length >= args.min_length:
		# get UMI
		umi_seq = seq[0:args.umi]
		for i in range(len(umi_seq)): umi_counter[i][umi_seq[i]] += 1
		
		# write trimmed read
		seq = seq[header_length:]
		name += ':' + umi_seq
		args.outfile.write(writefq(name, seq, qualities[header_length:]))

# print read count
print('%i reads processed\n%i trimmed below threshold by bcl2fastq\n' % (read_counter, overtrimmed_counter), file = sys.stderr)

# print UMI base frequency
alphabet = sorted(set.union(*(set(i) for i in umi_counter)))
print('UMI base frequency by position', file = sys.stderr)
print('\t'.join(alphabet), file = sys.stderr)
for pos in umi_counter:
	print('\t'.join(str(pos[base]) for base in alphabet), file = sys.stderr)
print('', file = sys.stderr)

# print insert length frequency
print('trimmed read length counts', file = sys.stderr)
for (length, count) in sorted(read_length_counter.items()):
	print('%i\t%i' % (length, count), file = sys.stderr)

