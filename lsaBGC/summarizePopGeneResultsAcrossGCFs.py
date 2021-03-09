import os
import sys
import argparse

def popgen_summarize(parent_dir, output):
	try:
		assert(os.path.isdir(parent_dir))
	except:
		sys.stderr.write('Error validating %s with BGReport results exists. Exiting now ...\n')

	out_handle = open(output, 'w')
	first = True
	for d in os.listdir(parent_dir):
		bgreport_result = parent_dir + d + '/Ortholog_Group_Information.txt'
		if os.path.isfile(bgreport_result):
			with open(bgreport_result) as obr:
				for i, line in enumerate(obr):
					line = line.strip()
					ls = line.split('\t')
					if i == 0 and first:
						out_handle.write('\t'.join(['gcf_id'] + ls) + '\n')
					elif i > 0:
						out_handle.write('\t'.join([d, d + '_|_' + ls[0]] + ls[1:]) + '\n')
			first = False
	out_handle.close()

if __name__ == '__main__':
	#Parse arguments.

	parser = argparse.ArgumentParser(description="""
	Program to summarize and aggregate results from BGReport run for each GCF.
	""")

	parser.add_argument('-p', '--parent_dir', help='Path to the parent directory where to search for results of BGReport runs.', required=True)
	parser.add_argument('-o', '--output', help='Path to output aggregate file.', required=False, default='bgreport_summary.txt')
	args = parser.parse_args()

	popgen_summarize(args.parent_dir, args.output)