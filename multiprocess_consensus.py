#! /usr/bin/python

import os;
import sys;
import multiprocessing;

from consensus_from_mpileup import *;

if __name__ == "__main__":
	if (len(sys.argv) < 5):
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\t%s <reference_file_path> coverage_threshold <collective_output_file> <{sb}am_file_1> [<{sb}am_file_2> <{sb}am_file_3> ...]\n' % sys.argv[0]);
		sys.stderr.write('\t(If <collective_output_file> is equal to "-", no files will be written to disk.)\n');
		exit(1);
	
	if (len(sys.argv) >= (5 + multiprocessing.cpu_count())):
		sys.stderr.write('Number of tasks to run in parallel is larger than number of cores.\n');
		sys.stderr.write('This case is currently not handled by the script.\n');
		sys.stderr.write('Please reduce the list of SAM/BAM files to contain at most %d files.\n' % multiprocessing.cpu_count());
		exit(1);
	
	reference_file = sys.argv[1];
	coverage_threshold = int(sys.argv[2]);
	collective_output_file = sys.argv[3];
	sam_files = sys.argv[4:];
	
	processes = [];
	
	i = 0;
	for sam_file in sam_files:
		output_prefix = os.path.splitext(sam_file)[0] if (collective_output_file != '-') else '';
		p = multiprocessing.Process(target=main, args=(sam_file, reference_file, coverage_threshold, output_prefix, i,))
		processes.append(p);
		p.start();
		i += 1;
	
	for p in processes:
		p.join();

	if (collective_output_file != '-'):
		CollectSummaries(sam_files, collective_output_file);
