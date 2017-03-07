#! /usr/bin/env python3

import collections, argparse, warnings, sys, pysam

# parse arguments
parser = argparse.ArgumentParser(description = 'given a BED file of genome regions expected to produce unwanted reads and a BAM file containing read alignments, mark any reads inside the unwanted genome regions', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-u', '--upstream', action = 'store', type = int, default = 500, help = 'number of bases upstream to extend unwanted regions')
parser.add_argument('-d', '--downstream', action = 'store', type = int, default = 0, help = 'number of bases downstream to extend unwanted regions')
parser.add_argument('-r', '--remove', action = 'store_true', help = 'remove blacklisted alignments instead of marking them')
parser.add_argument('-c', '--count', action = 'store_true', help = 'display count of hits per region')
parser.add_argument('region_bed', action = 'store', type = argparse.FileType('r'), help = 'BED file containing target regions')
parser.add_argument('in_file', action = 'store', nargs = '?', default = '-', help = 'input BAM')
parser.add_argument('out_file', action = 'store', nargs = '?', default = '-', help = 'output BAM')
args = parser.parse_args()

# parse region list into fancy dictionary
# top level: tuple by strand (forward, reverse)
# next level: dictionary by chromosome
# next level: list of regions
# each region is a tuple of (extended start position, extended end position, ID) where the ID is the chromosome name and original coordinates
regions = [collections.defaultdict(list), collections.defaultdict(list)]
region_names = collections.OrderedDict()
region_hits = collections.Counter()
for line in args.region_bed:
	if line.startswith('#'): continue # comment line
	fields = line.rstrip().split()
	
	# coordinates
	if len(fields) < 3: continue # no coordinates
	chrom, left, right = (fields[0], int(fields[1]) + 1, int(fields[2])) # adjust left by 1 because BED is 0-based start, 1-based end
	
	# name
	if len(fields) >= 4 and fields[3] and fields[3] != '.':
		name = fields[3]
	else:
		name = None
	
	# strand
	is_reverse = None
	if len(fields) >= 6:
		if fields[5] == '+':
			is_reverse = False
		elif fields[5] == '-':
			is_reverse = True
	
	# ID
	if is_reverse:
		region_id = '%s:%i-%i' % (chrom, right, left) # reverse coordinate order for reverse strand
	else:
		region_id = '%s:%i-%i' % (chrom, left, right)
	
	# add to region name lookup
	if region_id in region_names:
		warnings.warn('warning: %s defined twice; keeping %s and ignoring %s' % (region_id, ('first occurrence' if region_names[region_id] is None else region_names[region_id]), ('second occurrence' if name is None else name)))
	else:
		region_names[region_id] = name
	
	# add to region coordinate structure; non-directional regions go on both strands
	if is_reverse is None or is_reverse is False:
		regions[0][chrom].append((left - args.upstream, right + args.downstream, region_id))
	if is_reverse is None or is_reverse is True:
		regions[1][chrom].append((left - args.downstream, right + args.upstream, region_id))

print('read %i regions' % len(region_names), file = sys.stderr)


# read and mark alignments
in_bam = pysam.Samfile(args.in_file, 'rb')
out_bam = pysam.Samfile(args.out_file, 'wb', template = in_bam)

n_alignment = n_blacklisted = 0
for alignment in in_bam:
	n_alignment += 1
	blacklisted = False
	for region_left, region_right, region_id in regions[alignment.is_reverse][alignment.reference_name]:
		if alignment.reference_start >= region_left and alignment.reference_end <= region_right:
			blacklisted = True
			region_hits[region_id] += 1
	n_blacklisted += blacklisted
	if not (args.remove and blacklisted): out_bam.write(alignment)

print('read %i alignments, blacklisted %i' % (n_alignment, n_blacklisted), file = sys.stderr)

if args.count and n_blacklisted > 0:
	print('\ncoordinates', 'hits', file = sys.stderr)
	for region_id, region_name in region_names.items():
		hits = region_hits[region_id]
		if hits > 0: print('\t'.join((region_id, str(hits))), file = sys.stderr)

