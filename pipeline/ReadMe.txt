#######################
#
#	step 1: from compressed pair-end *fq.gz files to remove-dupl acgt-count results
#		1_brat_bw_one_step_one_dir_nodupl_Col0.pl
########################################################################################
1_brat_bw_one_step_one_dir_nodupl_Col0.pl <indir> <output_pre>

The raw data for each sample downloaded from BGI is named
1.fq.gz
2.fa.gz

So the script I wrote is specifically for these two file names. You can change it a little
for other sample(file) name.

Each sample has its own directory and should only contains the two files.

Finally, we can get 
<output_pre>_forw.txt
<output_pre>_rev.txt
The actg-count results. 

NOTE:
1, we need change the path of
$brat_bw
$trim_tool
$ref
$fas_ref
$remove_dupl
$acgt_count
to use the script.

for trim, the trim of brat_bw-2.0.1 itself does NOT work for pair-ended reads that mate1 and mate2 have different
read length. For that case, we need trim for latest version of BRAT.

2, we can run several samples at one time, seems one sample takes ~600M memory.

 


#######################
#
#	step 1.5: manually move all *_forw.txt *_rev.txt files to a directory (assume we name it nodupl_acgt_count)
#		
########################################################################################



#######################
#
#	step 2: calculate cytosine coverage 
#		2_batch_cal_cytosines_coverage.pl	
########################################################################################

2_batch_cal_cytosines_coverage.pl takes one input_dir and one output_dir
input_dir is nodupl_acgt_count from step2, 
assume output_dir is cytosine_coverage_nodupl

2_batch_cal_cytosines_coverage.pl <nodupl_acgt_count> <cytosine_coverage_nodupl>

each sample will generate two files in cytosine_coverage_nodupl:
one is *_cytosines_coverage_info.txt
the other is *_depth_distribution.txt

You can see real output files in sample_file/cytosine_coverage_nodupl directory.

From *_depth_distribution.txt we can know max depth of the sample. As I see, after remove-dupl, usually max_depth <= 200 for one single sample. 

NOTE:
2_batch_cal_cytosines_coverage.pl will use the script 
cal_cytosines_coverage_for_acgt_count_output_v0.1.pl


#######################
#
#	step 3: summary cytosine coverage 
#		3_summary_report_Kaust.pl	
########################################################################################

3_summary_report_Kaust.pl <cytosine_coverage_nodupl> <output_file_name>

3_summary_report_Kaust.pl use cytosine_coverage_nodupl directory in step3 as input directory to summarize coverage, methylation level and 
error rate(methylation level for chic).
The sample of output file is sample_file/SH_C24_summary_nodupl.txt


#######################
#
#	step 4: batch call methylation
#		4_batch_call_methylation.pl	
########################################################################################
4_batch_call_methylation.pl <input_summary_file> <indir> <outdir> 

input_summary_file is the output_file_name from step4, <indir> is nodupl_acgt_count, 
<outdir> is directory for output. We can name it isMeth_nodupl.

three files for each sample will be generated:

1, *_chrC_error_cutoff_separately_called.txt
2, *_depth_mC.txt
3, *_isMeth_chrC_error_separately_called.txt

The sample files are in sample_file/isMeth_nodupl


NOTE:
1, call_methylation_CG_CHG_CHH_separated_v0.2_chrC.pl will be used

2, we may need modify $max_dep (default is 200) in 4_batch_call_methylation.pl to make sure $max_dep is really >= real_max_depth, 
especially for pooled samples(combine two replicates toghter).

################################################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################
#
#			PART II : DMR
#
################################################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################


#######################
#
#	step 5: generate DMC
#		5_generate_p_value_Fisher_for_brat_bw_acgt_count_output_v1.1.pl	
########################################################################################

5_generate_p_value_Fisher_for_brat_bw_acgt_count_output_v1.1.pl <WT_isMeth> <WT_pre> <mut_isMeth> <mut_pre>  <outdir>

use isMeth files from step4

For example
5_generate_p_value_Fisher_for_brat_bw_acgt_count_output_v1.1.pl /path/to/colA_nodupl2_isMeth_chrC_error_separately_called.txt   colA_nodupl   /path/to/nrpd1-3A_nodupl_isMeth_chrC_error_separately_called.txt   nrpd1-3A   .

And the output is named as colA_nodupl_vs_nrpd1-3A_diff_bases_P0.05_OnlyBothMinDep4.txt defaultly.

Note:
1,In this step, we can run for several samples.

2, In this script, 
default
$p_cutoff = 0.05;
$min_depth = 4;

means that DMC is defined as p-value of Fisher's exact test <= 0.05
$min_depth = 4 means that only those cytosines with depth in WT and mutant both >= 4 will be considered for Fisher's exact test.


#######################
#
#	step 5.5:  manually move all  *_diff_bases_P0.05_OnlyBothMinDep4.txt into one directory (assume named DMC_db). DMC_db should
#			contains *_diff_bases_P0.05_OnlyBothMinDep4.txt only, no other files.			
########################################################################################



#######################
#
#	step 6: batch generate raw list
#		6_batch_select_DMR_IDM1_paper_method_v1.0_sliding_win.pl	
########################################################################################

6_batch_select_DMR_IDM1_paper_method_v1.0_sliding_win.pl <indir> <outdir> <p_value> <dep> 

<indir> is DMC_db,
<outdir> is just directory for output list

<p_value> should be 0.05 which is consistent with step5
<dep> should be 4 which is also consistent with step5


NOTE:
select_DMR_output_list_v1.0_sliding_window.pl
will be called.

#######################
#
#	step 7: calculation methylation level in raw list regions
#			
########################################################################################

7_batch_cal_meth_level_in_Meth_file_for_bed_list_modify.pl <bed_indir> <outdir> 



<bed_indir> is the <outdir> of step6

NOTE:
1, We need most modification for this script to complete the task
(1), line 33, $WT_isMeth_file need be modified
(2), line 39, the hash %pres need to be modified for each sample
(3), line 75, the hash %posts need to be modified to be consist with raw_list file name.

2, For hyper list and hypo list, the output format is different that high methylation level is before the low methylation level.

For example:
For hyper list, mut_methylation_level is prior to WT_methylation_level.
For hypo list, WT_methylation_level is prior to mut_methylation_level.


3, count_meth_level_in_Meth_file_for_bed_list_v0.1.pl
will be called.


#######################
#
#	step 8: filter
#			8_batch_filter_based_on_len100_CG30_CHG15_CHH10.pl
########################################################################################
8_batch_filter_based_on_len100_CG30_CHG15_CHH10.pl

batch_filter_based_on_len100_CG30_CHG15_CHH10.pl  <indir> <outdir>

the <indir> is outdir of Step7.

The filter rules are;
(1) length >= 100bp

AND
(2)
diff_CG_meth_level >= 0.3 
or 
diff_CHG_meth_level >= 0.15 
or 
diff_CHH_meth_level >= 0.1 


NOTE:
filter_len100_CG30_CHG15_CHH10.pl
will be called

