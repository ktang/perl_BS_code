#	my $R_cmd = "R --slave --vanilla --args $raw_file $fdr_file < $R_script_sub ";

args = commandArgs(trailingOnly = TRUE)
p_file      = args[1]
p_adj_file  = args[2]

if(!file.exists(p_file) ){ print ("file\n"); q("no")}
if(file.exists(p_adj_file)) { print("adj_file exits\n\n");q("no")}

raw_p = read.table(p_file, head=T, sep = "\t");

p_adj = p.adjust(raw_p[,11 ], method = "BH");
p_adj = cbind(raw_p, p_adj );

colnames(p_adj) = c( colnames(raw_p), "FDR");

print ("write file\n\n")
write.table(p_adj, file =p_adj_file, quote = F , row.names =F, col.names=T, sep="\t")

print ("done\n\n");


q("no");
