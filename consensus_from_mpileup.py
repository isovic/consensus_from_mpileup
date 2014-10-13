#! /usr/bin/python

import os;
import sys;
import operator;
import subprocess;

def increase_in_dict(dict_counter, value):
	try:
		dict_counter[value] += 1;
	except:
		dict_counter[value] = 1;

def process_mpileup_line(line, line_number, ret_variant_list, ret_snp_count, ret_insertion_count, ret_deletion_count, ret_num_undercovered_bases, ret_num_called_bases, ret_num_correct_bases, ret_coverage_sum, coverage_threshold):
	# Split the line, and perform a sanity check.
	split_line = line.strip().split('\t');
	if (len(split_line) != 6):
		return;
	
	ref_name = split_line[0];
	position = split_line[1];
	ref_base = split_line[2];
	coverage = split_line[3];
	original_bases = split_line[4];
	qualities = split_line[5];
	
	ret_coverage_sum[0] += int(coverage);
	
	bases = '';
	
	# Replace the '.' and ',' signs with the actual reference base.
	i = 0;
	while (i < len(original_bases)):
		if (original_bases[i] == '.' or original_bases[i] == ','):
			bases += ref_base;
		else:
			bases += original_bases[i];
		i += 1;

	base_counts = {};
	insertion_count = 0;
	current_base_deletion_count = 0;
	deletion_count = 0;
	insertion_event_counts = {};
	deletion_event_counts = {};
	end_counts = 0;
	
	i = 0;
	while (i < len(bases)):
		base = bases[i];
		
		if (base == r'^'):
			# This is the starting position of a read. It encodes two
			# symbols: '^' marking the read start and a char marking the
			# mapping quality of the read.
			#increase_in_dict(base_counts, bases[i + 1].upper());
			i += 1;			# Increase only by 1, because we have i += 1 down there.
		elif (base == r'$'):
			# This marks the end of a read.
			end_counts += 1;
		elif (base == r'*'):
			# This is a deletion, just count it.
			current_base_deletion_count += 1;
		elif (base == r'-'):
			# This marks the occurance of deletions. It is a composite object
			# consisting of: the special character '-', the number of the deleted bases
			# and the actual bases that are deleted (these bases follow the current position).
			# In our approach, we ignore this case, because we count deletions one by one
			# through the '*' character.
			
			# Get the number of bases that need to be skipped in the string.
			j = (i + 1);
			while (bases[j] in '0123456789'):
				j += 1;
			num_bases = int(bases[(i + 1):j]);
			skip_bases = (j - i) + num_bases - 1;
			deletion_count += 1;
			deletion = bases[j : (j + num_bases)].upper();
			increase_in_dict(deletion_event_counts, deletion);
			# Skip the length of the numeric entry plus the actual number of bases
			# that need to be skipped.
			i += skip_bases;
		elif (base == r'+'):
			# This marks the occurance of an insertion. It is a composite object
			# consisting of: the special character '+', the number of the inserted bases
			# and the actual bases that are inserted (these bases follow the current position).
			# Similar to the deletion marking, but here we actually care about the bases,
			# and we need to make an allele aware count.

			# Get the number of bases that are inserted;
			j = (i + 1);
			while (bases[j] in '0123456789'):
				j += 1;
			num_bases = int(bases[(i + 1):j]);
			skip_bases = (j - i) + num_bases - 1;
			insertion_count += 1;
			insertion = bases[j : (j + num_bases)].upper();
			increase_in_dict(insertion_event_counts, insertion);
			i += skip_bases;
		else:
			increase_in_dict(base_counts, bases[i].upper());
		i += 1;
	
	if (int(coverage) < coverage_threshold):
		ret_num_undercovered_bases[0] += 1;
	else:
		ret_num_called_bases[0] += 1;
		
		### Handling base consensus.
		sorted_base_counts = sorted(base_counts.items(), key=operator.itemgetter(1));
		most_common_base_count = sorted_base_counts[-1][1];
		# Allow for the case where there are multiple equally good choices.
		# In this case, we prefer the choice which is equal to the reference.
		is_good = False;
		for base_count in sorted_base_counts:
			if (base_count[1] == most_common_base_count):
				if (base_count[0] == ref_base):
					is_good = True;
					break;
		if (is_good == False):
			ret_snp_count[0] += 1;
			ret_variant_list.append(line_number);
		else:
			ret_num_correct_bases[0] += 1;

		### Handling indel consensus.
		non_indel_coverage = int(coverage) - end_counts - deletion_count - insertion_count;

		# Sanity check, just to see if there actually were any insertions (to avoid index out of bounds error).
		# If there are insertions, get the most common one.
		if (len(insertion_event_counts.keys()) > 0):
			sorted_insertion_counts = sorted(insertion_event_counts.items(), key=operator.itemgetter(1));
			most_common_insertion_count = sorted_insertion_counts[-1][1];
			most_common_insertion_length = len(sorted_insertion_counts[-1][0]);
			insertion_unique = True if (sum([int(insertion_count[1] == most_common_insertion_count) for insertion_count in sorted_insertion_counts]) == 1) else False;
		else:
			most_common_insertion_count = 0;
			most_common_insertion_length = 0;
			insertion_unique = False;
		
		# Sanity check, just to see if there actually were any deletions (to avoid index out of bounds error).
		# If there are deletions, get the most common one.
		if (len(deletion_event_counts.keys()) > 0):
			sorted_deletion_counts = sorted(deletion_event_counts.items(), key=operator.itemgetter(1));
			most_common_deletion_count = sorted_deletion_counts[-1][1];
			most_common_deletion_length = len(sorted_deletion_counts[-1][0]);
			deletion_unique = True if (sum([int(deletion_count[1] == most_common_deletion_count) for deletion_count in sorted_deletion_counts]) == 1) else False;
		else:
			most_common_deletion_count = 0;
			most_common_deletion_length = 0;
			deletion_unique = False;
		
		if (most_common_insertion_count > most_common_deletion_count and most_common_insertion_count > non_indel_coverage):
			# In this case, insertions are a clear winner.
			if (insertion_unique == True):
				ret_insertion_count[0] += most_common_insertion_length;
				ret_num_called_bases[0] += most_common_insertion_length;
		elif (most_common_deletion_count > most_common_insertion_count and most_common_deletion_count > non_indel_coverage):
			# In this case, deletions are a clear winner.
			if (deletion_unique == True):
				ret_deletion_count[0] += most_common_deletion_length;
				return most_common_deletion_length;
		else:
			# In this case, either the base count consensus wins, or the
			# insertion/deletion count is ambiguous.
			pass;

	return 0;

def process_mpileup(alignments_path, mpileup_path, coverage_threshold, output_prefix, thread_id=0):
	fp = None;
	try:
		fp = open(mpileup_path, 'r');
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % mpileup_path;
		return;
	
	ret_variant_list = [];
	ret_snp_count = [0];
	ret_insertion_count = [0];
	ret_deletion_count = [0];
	ret_num_undercovered_bases = [0];
	ret_num_called_bases = [0];
	ret_num_correct_bases = [0];
	ret_coverage_sum = [0];
	
	lines = fp.readlines();
	fp.close();
	
	i = 0;
	j = 0;
	while (i < len(lines)):
		line = lines[i];
		
		if (thread_id == 0):
			if ((j % 100) == 0):
				sys.stdout.write('\r[%d] snps = %d, insertions = %d, deletions = %d, undercovered = %d, coverage = %.2f' % (i, ret_snp_count[0], ret_insertion_count[0], ret_deletion_count[0], ret_num_undercovered_bases[0], (float(ret_coverage_sum[0])/float((i + 1)))));
				sys.stdout.flush();
		
		num_bases_to_skip = process_mpileup_line(line, i, ret_variant_list, ret_snp_count, ret_insertion_count, ret_deletion_count, ret_num_undercovered_bases, ret_num_called_bases, ret_num_correct_bases, ret_coverage_sum, coverage_threshold);
		
		i += num_bases_to_skip;
		i += 1;
		j += 1;

	summary_lines = '';
	summary_lines += 'alignments_file: %s\n' % alignments_path;
	summary_lines += 'mpileup_file: %s\n' % mpileup_path;
	summary_lines += 'coverage_threshold: %d\n' % coverage_threshold;
	summary_lines += 'snp_count: %d\n' % ret_snp_count[0];
	summary_lines += 'insertion_count: %d\n' % ret_insertion_count[0];
	summary_lines += 'deletion_count: %d\n' % ret_deletion_count[0];
	summary_lines += 'num_undercovered_bases: %d\n' % ret_num_undercovered_bases[0];
	summary_lines += 'num_called_bases: %d\n' % ret_num_called_bases[0];
	summary_lines += 'num_correct_bases: %d\n' % ret_num_correct_bases[0];
	summary_lines += 'average_coverage: %.2f\n' % ((float(ret_coverage_sum[0])/float((i + 1))));
	
	sys.stdout.write(summary_lines);
	sys.stdout.write('\n');
	
	if (output_prefix != ''):
		summary_file = output_prefix + '.conssum';
		try:
			fp_sum = open(summary_file, 'w');
			fp_sum.write(summary_lines);
			fp_sum.close();
		except IOError:
			sys.stderr.write('ERROR: Could not open file "%s" for writing!\n' % (summary_file));

def main(alignments_path, reference_path, coverage_threshold, output_prefix, thread_id=0):
	# Sanity checking the existence of the file, and the correctness of its extension.
	# Also, if input file is a SAM file, then convert it to a sorted BAM.
	alignments_path_bam = alignments_path;
	if (os.path.exists(alignments_path) == False):
		print 'ERROR: File "%s" does not exist!' % alignments_path;
		return;
	if (alignments_path.endswith('sam')):
		# Determine the path where the new BAM file will be generated.
		alignments_path_bam = os.path.dirname(alignments_path) + '/' + os.path.splitext(os.path.basename(alignments_path))[0] + '.bam'
		# Check if a BAM file with the given name already exists.
		if (os.path.exists(alignments_path_bam) == False):
			# Convert the SAM file to a sorted BAM file.
			command = 'samtools view -bS %s | samtools sort - %s' % (alignments_path, os.path.splitext(alignments_path_bam)[0]);
			subprocess.call(command, shell='True');
			# Create the BAM index file.
			command = 'samtools index %s %s.bai' % (alignments_path_bam, alignments_path_bam);
			subprocess.call(command, shell='True');		
	elif (alignments_path.endswith('bam') == False):
		print 'ERROR: File extension needs to be either .sam or .bam! Input file path: "%s".' % alignments_path;
		return;
	
	# Convert the sorted BAM file to a mpileup file if it doesn't exist yet.
	if (os.path.exists(('%s.mpileup' % alignments_path_bam)) == False):
		command = 'samtools mpileup -B -d 1000000 -Q 0 -A -f %s %s > %s.mpileup' % (reference_path, alignments_path_bam, alignments_path_bam);
		subprocess.call(command, shell='True');

	print 'Processing file "%s"...' % alignments_path;
	print 'Coverage threshold: %d' % coverage_threshold;
	process_mpileup(alignments_path, ('%s.mpileup' % alignments_path_bam), coverage_threshold, output_prefix, thread_id);

def CollectSummaries(sam_files, collective_output_file):
	fp_collect = None;
	
	try:
		fp_collect = open(collective_output_file, 'w');
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for writing!\n' % collective_output_file);
		return;
	
	for sam_file in sam_files:
		summary_file = os.path.splitext(sam_file)[0] + '.conssum';
		
		try:
			fp_sum = open(summary_file, 'r');
			lines = fp_sum.readlines();
			fp_sum.close();
		except IOError:
			sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % summary_file);
			continue;
		
		fp_collect.write(''.join(lines) + '\n');
	
	fp_collect.close();

if __name__ == "__main__":
	if (len(sys.argv) < 5):
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\t%s <reference_file_path> coverage_threshold <collective_output_file> <{sb}am_file_1> [<{sb}am_file_2> <{sb}am_file_3> ...]\n' % sys.argv[0]);
		sys.stderr.write('\t(If <collective_output_file> is equal to "-", no files will be written to disk.)\n');
		exit(1);
	
	reference_file = sys.argv[1];
	coverage_threshold = int(sys.argv[2]);
	collective_output_file = sys.argv[3];
	sam_files = sys.argv[4:];
	
	processes = [];
	
	for sam_file in sam_files:
		output_prefix = os.path.splitext(sam_file)[0] if (collective_output_file != '-') else '';
		main(sam_file, reference_file, coverage_threshold, output_prefix, 0);

	if (collective_output_file != '-'):
		CollectSummaries(sam_files, collective_output_file);
