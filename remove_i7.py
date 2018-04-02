#! /usr/bin/env python3

import argparse, sys

DATA_SECTION_HEADER = '[Data]'
INDEX_ID_HEADER = 'I7_Index_ID'
INDEX_SEQ_HEADER = 'index'
DEFAULT_DUMMY_SEQ = 'NN'

parser = argparse.ArgumentParser('remove i7 index sequences from an Illumina Experiment Manager sample sheet')
parser.add_argument('-s', '--replacement_sequence', action = 'store', type = str, default = DEFAULT_DUMMY_SEQ, help = 'dummy sequence to replace i7 sequences with, default %s' % DEFAULT_DUMMY_SEQ)
parser.add_argument('original_sheet', action = 'store', nargs = '?', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('new_sheet', action = 'store', nargs = '?', type = argparse.FileType('w'), default = sys.stdout)
args = parser.parse_args()

for line in args.original_sheet:
	args.new_sheet.write(line)
	if line.startswith(DATA_SECTION_HEADER): break

header_line = next(args.original_sheet)
data_headers = header_line.rstrip().split(',')
args.new_sheet.write(header_line)

n_columns = len(data_headers)
which_index_id = [i for i in range(len(data_headers)) if data_headers[i] == INDEX_ID_HEADER]
assert(len(which_index_id) == 1)
which_index_seq = [i for i in range(len(data_headers)) if data_headers[i] == INDEX_SEQ_HEADER]
assert(len(which_index_seq) == 1)

for line in args.original_sheet:
	fields = line.rstrip().split(',')
	assert(len(fields) == n_columns)
	fields[which_index_id[0]] = 'none'
	fields[which_index_seq[0]] = args.replacement_sequence
	args.new_sheet.write(','.join(fields) + '\n')

