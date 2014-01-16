###########
#	1 brat bw
#####################




###########
#	2 cal cytosine coverage
#####################
time perl ~/Kai_BS/myown/cytosines_coverage/batch_cal_cytosines_coverage.pl acgt_count/ cytosine_coverage/


###########
#	3 summary cytosine coverage
#####################
perl ~/Kai_BS/myown/cytosines_coverage/summary_report_Kaust.pl cytosine_coverage/ low_ros1_4_pool_meth_level.txt


###########
#	4 batch call isMeth
#####################
/Users/tang58/Kai_BS/myown/call_methylation_CG_CHG_CHH/batch_call_methylation.pl 
 <input_summary_file> <indir> <outdir>
