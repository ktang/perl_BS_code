Gene:
(]
All_Gene
0_1kb
1_2kb
2_3kb
3_4kb
4kb_up

~/Kai_BS/for_publish/9_feature_meth_level_db/Gene_13:34:57_N=568$ less all_protein_coding_gene_coordinate_bed.txt | perl -lane ' $l = $F[2] - $F[1] + 1; print  if ($F[0] eq "chr" or ( $l>0 and $l <= 1000 ))' > Gene_0_1kb_coordinate_bed.txt
~/Kai_BS/for_publish/9_feature_meth_level_db/Gene_13:35:23_N=569$ less all_protein_coding_gene_coordinate_bed.txt | perl -lane ' $l = $F[2] - $F[1] + 1; print  if ($F[0] eq "chr" or ( $l>1000 and $l <= 2000 ))' > Gene_1_2kb_coordinate_bed.txt
~/Kai_BS/for_publish/9_feature_meth_level_db/Gene_13:35:44_N=570$ less all_protein_coding_gene_coordinate_bed.txt | perl -lane ' $l = $F[2] - $F[1] + 1; print  if ($F[0] eq "chr" or ( $l>2000 and $l <= 3000 ))' > Gene_2_3kb_coordinate_bed.txt
~/Kai_BS/for_publish/9_feature_meth_level_db/Gene_13:35:57_N=571$ less all_protein_coding_gene_coordinate_bed.txt | perl -lane ' $l = $F[2] - $F[1] + 1; print  if ($F[0] eq "chr" or ( $l>3000 and $l <= 4000 ))' > Gene_3_4kb_coordinate_bed.txt
~/Kai_BS/for_publish/9_feature_meth_level_db/Gene_13:36:13_N=572$ less all_protein_coding_gene_coordinate_bed.txt | perl -lane ' $l = $F[2] - $F[1] + 1; print  if ($F[0] eq "chr" or ( $l>4000  ))' > Gene_4kb_up_coordinate_bed.txt
~/Kai_BS/for_publish/9_feature_meth_level_db/Gene_13:36:32_N=573$ mv all_protein_coding_gene_coordinate_bed.txt Gene_all_coordinate_bed.txt



TE:
(]
All_TE
0_0.5kb
0.5_1kb
1_2kb
2_4kb
4kb_up
~/Kai_BS/for_publish/9_feature_meth_level_db/TE_13:56:10_N=580$ mv all_TE_coordinate_bed.txt TE_all_coordinate_bed.txt
~/Kai_BS/for_publish/9_feature_meth_level_db/TE_13:58:07_N=587$ less TE_all_coordinate_bed.txt |  perl -lane ' $l = $F[2] - $F[1] + 1; print  if ($F[0] eq "chr" or ( $l>0 and $l <= 500 )) ' > TE_0_0.5kb_coordinate_bed.txt
~/Kai_BS/for_publish/9_feature_meth_level_db/TE_14:00:12_N=588$ less TE_all_coordinate_bed.txt |  perl -lane ' $l = $F[2] - $F[1] + 1; print  if ($F[0] eq "chr" or ( $l>500 and $l <= 1000 )) ' > TE_0.5_1kb_coordinate_bed.txt
~/Kai_BS/for_publish/9_feature_meth_level_db/TE_14:00:41_N=589$ less TE_all_coordinate_bed.txt |  perl -lane ' $l = $F[2] - $F[1] + 1; print  if ($F[0] eq "chr" or ( $l>1000 and $l <= 2000 )) ' > TE_1_2kb_coordinate_bed.txt
~/Kai_BS/for_publish/9_feature_meth_level_db/TE_14:00:54_N=590$ less TE_all_coordinate_bed.txt |  perl -lane ' $l = $F[2] - $F[1] + 1; print  if ($F[0] eq "chr" or ( $l>2000 and $l <= 4000 )) ' > TE_2_4kb_coordinate_bed.txt
~/Kai_BS/for_publish/9_feature_meth_level_db/TE_14:01:09_N=591$ less TE_all_coordinate_bed.txt |  perl -lane ' $l = $F[2] - $F[1] + 1; print  if ($F[0] eq "chr" or ( $l>4000  )) ' > TE_4kb_up_coordinate_bed.txt



~/Kai_BS/for_publish/9_feature_meth_level_db_21:24:32_N=505$ time perl ../oneUse_db_src/count_feature_number_density_using_sliding_window.pl Gene/Gene_all_coordinate_bed.txt . TAIR10_all_Gene

real	0m0.238s
user	0m0.228s
sys	0m0.005s
~/Kai_BS/for_publish/9_feature_meth_level_db_21:24:46_N=506$ less TAIR10_all_Gene_number_WinSize100000_sliding100000.txt 
~/Kai_BS/for_publish/9_feature_meth_level_db_21:25:13_N=507$ time perl ../oneUse_db_src/count_feature_number_density_using_sliding_window.pl TE . TAIR10_all_TE
TE/                         TE_classify_by_SuperFamily/ 
~/Kai_BS/for_publish/9_feature_meth_level_db_21:25:13_N=507$ time perl ../oneUse_db_src/count_feature_number_density_using_sliding_window.pl TE/TE_all_coordinate_bed.txt . TAIR10_all_TE

real	0m0.269s
user	0m0.258s
sys	0m0.005s

