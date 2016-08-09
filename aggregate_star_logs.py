#! /usr/bin/env python3

import sys, re, argparse

parser = argparse.ArgumentParser(description = 'given a list of STAR log files ("Log.final.out"), aggregate the numbers into one table\nexpects they all have the same fields in the same order')
parser.add_argument('log_files', nargs = '*')
args = parser.parse_args()

fields = []
values = []
first_file = True

parse_expr = re.compile('^ +(.+) \|\t(.+)$')

# parse input
for log_file in args.log_files:
	i = 0
	for line in open(log_file):
		expr_match = parse_expr.match(line)
		if expr_match is not None:
			field, value = expr_match.groups()
			if first_file:
				assert len(fields) == len(values) == i
				fields += [field]
				values += [[value]]
			else:
				assert fields[i] == field
				values[i] += [value]
			i += 1
	first_file = False

# write output
print('\t'.join(['library'] + fields))
for i in range(len(args.log_files)):
	print('\t'.join([args.log_files[i]] + [value_list[i] for value_list in values]))

