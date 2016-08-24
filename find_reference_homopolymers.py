#!/usr/bin/env python3

# simplified first pass:
# don't look for complement on opposite strand
# don't allow mismatches (!!!)
# don't allow IUPAC base codes outside the usual ACGT

import re, argparse, sys

COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

CHR_HEADER = re.compile('> *(.+)')


parser = argparse.ArgumentParser(description = 'given reference sequence in FASTA format, produce a list in BED format of homopolymers of a given base')
parser.add_argument('-l', '--length', action = 'store', type = int, default = 10, help = 'minimum length of homopolymer')
parser.add_argument('-p', '--proportion', action = 'store', type = float, default = 0.8, help = 'minimum proportion of target base in homopolymer stretch')
parser.add_argument('base', action = 'store', type = str, help = 'target base')
parser.add_argument('in_fastq', action = 'store', nargs = '?', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('out_bed', action = 'store', nargs = '?', type = argparse.FileType('w'), default = sys.stdout)
args = parser.parse_args()

target_base = args.base.upper()
try:
	complement = COMPLEMENTS[target_base]
except KeyError:
	raise NotImplementedError(('unknown base %s; options: ' % target_base) + ' '.join(COMPLEMENTS.keys()))

current_pos = 0
buffer_start_pos = None
chrom = None
seq_buffer = ''
target_count = non_target_count = 0
homopolymers_found = 0
total_homopolymer_length = 0
total_reference_length = 0


def flush_buffer ():
	global seq_buffer, args, chrom, buffer_start_pos, seq_buffer, homopolymers_found, total_homopolymer_length, target_count, non_target_count # yes I know this is ugly
	if len(seq_buffer) >= args.length:
		print('\t'.join((chrom, str(buffer_start_pos - 1), str(buffer_start_pos + len(seq_buffer))))) # output in BED format; note the 0-based start, 1-based end
		homopolymers_found += 1
		total_homopolymer_length += len(seq_buffer)
	buffer_start_pos = None
	seq_buffer = ''
	target_count = non_target_count = 0


for fastq_line in args.in_fastq:
	fastq_line.rstrip()
	if len(fastq_line) == 0: continue # empty line
	
	match = re.match(CHR_HEADER, fastq_line)
	if match: # chromosome header
		flush_buffer()
		current_pos = 0
		chrom = match.groups()[0]
	
	else: # sequence data
		for base in fastq_line:
			current_pos += 1
			total_reference_length += 1
			if base.upper() == target_base:
				if len(seq_buffer) == 0: buffer_start_pos = current_pos # starting new run
				seq_buffer += base
				target_count += 1
			else:
				flush_buffer()

flush_buffer()
sys.stderr.write('found %i homopolymers comprising %i bases (%.6f%% of total reference length)\n' % (homopolymers_found, total_homopolymer_length, 100 * total_homopolymer_length / total_reference_length))

