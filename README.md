Script that computes a simple consensus from given sequence alignments.
Alignments can be provided either as SAM or as sorted BAM files, after which they are converted to the mpileup format using SAMtools.
SAMtools needs to be installed.

Variants are called by a majority vote, and include no fancy statistics, just simple counting.
The resulting report contains counts for SNPs, insertions and deletions, as well as the number of called bases, number of bases under the coverage threshold, and the average coverage of the reference sequence.
Output is provided in two files: a summary '.txt' file with these measurements, '*.csv' file containing details on the detected variants/uncalled bases. The CSV file is not in the VCF format.

Usage:
     ./consensus_from_mpileup.py <reference_file_path> coverage_threshold <collective_output_file> <{sb}am_file_1> [<{sb}am_file_2> <{sb}am_file_3> ...]

Multiple SAM/BAM files can be provided, and will be ran sequentially one after another.
If instead of the concrete path provided by <collective_output_file> a character "-" is given, no output files will be written to disk, and only the summary will be output to stdout.

Example usages:
     ./consensus_from_mpileup.py escherichia_coli.fa 10 collective_summary.txt alignments.sam
     ./consensus_from_mpileup.py escherichia_coli.fa 10 collective_summary.txt alignments1.sam alignments2.bam alignments3.sam
     ./consensus_from_mpileup.py escherichia_coli.fa 10 - alignments1.sam alignments2.sam alignments3.sam

