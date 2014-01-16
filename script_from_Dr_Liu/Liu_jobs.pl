#!/usr/bin/perl -w
use strict;
# process bisulfite sequencing data
# step 1: check quality
#`perl check_quality.pl`;
#`perl move_quality_folder.pl quality_h20`;

# try -h 17
#perl check_quality.pl 17
#perl move_quality_folder.pl quality_h17

# step 2: trim by quality
#`perl trim_quality.pl`;

# step 3: convert fastq files to fasta files
#`perl change_fastq_to_fasta.pl`;

# step 4: remove adaptors & calculate base composition
#`perl batch_remove_adaptor.pl`;
#`perl calc_base_composition.pl`;

# step 5: make brat input files
#`perl batch_get_pairs.pl`;

# step 6: brat map
`perl batch_brat_map.pl`;
`perl batch_brat_single.pl`;

#brat_single did not finish for Lane1 due to power outage. re-run as follows:
#/home/renyi/archives/brat/brat-1.1.18/brat -r /mnt/disk4/renyi/bis_seq/col_ros2_11_09_10/references_names.txt -s Lane1_JKZ1_Col0/s_1_1_single.txt -m 2 -bs -o Lane1_JKZ1_Col0/s_1_1_single_brat.txt >> brat_single.log 2>&1
#/home/renyi/archives/brat/brat-1.1.18/brat -r /mnt/disk4/renyi/bis_seq/col_ros2_11_09_10/references_names.txt -A -s Lane1_JKZ1_Col0/s_1_2_single.txt -m 2 -bs -o Lane1_JKZ1_Col0/s_1_2_single_brat.txt >> brat_single.log 2>&1

# step 7: calculate coverage
`perl batch_calc_coverage.pl`;

#step 8: count on each position
`perl batch_count.pl`;

# step 9: calculate conversion error
`perl batch_error.pl`;

#step 10: prepare wig files for methylation positions
`perl batch_wig.pl`;

