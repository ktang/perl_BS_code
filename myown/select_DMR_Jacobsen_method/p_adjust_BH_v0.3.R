#setwd("/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/figure")
#my $cmd = "R --slave --vanilla --args $indir $input $outdir $label < R_cmd "
args = commandArgs()
p_file      = args[5]
p_adj_file  = args[6]

if(!file.exists(p_file) ){ print ("file\n"); q("no")}
if(file.exists(p_adj_file)) { print("adj_file exits\n\n");q("no")}

raw_p = read.table(p_file, head=T, sep = "\t");
# 1       2          3    4       5            6           7            8       9
#"chr", "start", "end", "p_CG","CG_hyper",  "p_CHG", "CHG_hyper", "p_CHH", "CHH_hyper"
p_adj_CG = p.adjust(raw_p[,4 ], method = "BH");
p_adj_CHG = p.adjust(raw_p[,6 ], method = "BH");
p_adj_CHH = p.adjust(raw_p[,8 ], method = "BH");

p_adj = cbind(raw_p[,1:3], p_adj_CG, raw_p[,5] ,  p_adj_CHG, raw_p[,7] , p_adj_CHH, raw_p[,9] );

colnames(p_adj) = c("chr", "start", "end", "p_adj_CG","CG_hyper",  "p_adj_CHG", "CHG_hyper", "p_adj_CHH", "CHH_hyper" )

print ("write file\n\n")
write.table(p_adj, file =p_adj_file, quote = F , row.names =F, col.names=T, sep="\t")

print ("done\n\n");


q("no");
