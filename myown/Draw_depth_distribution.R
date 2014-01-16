# barplot for nucleitide distribution
#this is version 1.1
# for position distribution, the ylim should be a little above the max.

#         1    2        3       4       5      6         7
#my $cmd = "R --slave --vanilla --args  $input  $output $label ". "< /Users/tang58/Kai_BS/myown/Draw_depth_distribution.R";

args = commandArgs()
#base_dist = read.table("base_dist.txt" ,header =T, sep ="\t")
#pos_dist = read.table("tail.txt" ,header =T, sep ="\t")
#len_dist = read.table("length_dist.txt" ,header =T, sep ="\t")

#pre = args[5]
#jpeg("len_dist.jpeg",width = 1300)

input  = args[5]
output = args[6]
label  = args[7]
main_cov = paste(label ,"cytosine depth distribution",sep=" ")
#output_per = paste(label,"_methylated_percentage.jpg",sep="")
main_per = paste("methylated percentage at each coverage level",label,sep="\n")


y = read.table(input,header=T , sep="\t")

#
#30-100
#101-200
#>=200

jpeg(filename = output,width = 800)

y31 = sum(y[ y[,1] >= 31 & y[,1] <= 100 ,2]);
y101 =sum(y[ y[,1] >= 101 & y[,1] <= 200 ,2]);
y201 = sum(y[ y[,1] >= 201  ,2]);

barplot(c(y[y[,1] <= 30,2] , y31, y101, y201),space = 0.2,col = "blue"
, names.arg = c(y[y[,1] <= 30,1], "31-100", "101-200", ">200") , main = main_cov,
xlab = "depth", ylab = "count", las = 3);


dev.off()




if(FALSE){
jpeg(filename = output_per,width = 800)
met30 = sum(y[ y[,1] >= 30,3])	
max = max(y[ y[,1] < 30,3]/y[y[,1] < 30,2] * 100, met30/y30 *100)
barplot(c(y[ y[,1] < 30,3]/y[y[,1] < 30,2] * 100, met30/y30 *100), space = 1,width = 10,col = "blue", names.arg = c(y[y[,1] < 30,1], ">=30") , main = main_per, ylim=c(0,max + 5) , ylab = "%"
,xlab = "coverage level")
dev.off()

# y50 = sum(y[ y[,1] >= 50,2])
# met50 = sum(y[ y[,1] >= 50,3])
#barplot(c(y[1:50,3]/y[1:50,2] * 100, met50/y50 *100), space = 1,width = 10,col = "blue", names.arg = c(0:49, ">=50") , main = "distribution_of_coverage", ylim=c(0,100) , ylab = "%")

#setwd("/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/downstream")
#len_dist = read.table("JKZ13_phd_coverage_info_perl.txt",head=T)
#y= len_dist
#y30 = sum(y[ y[,1] >= 30,2])
#barplot(c(y[y[,1] < 30,2] , y30),space = 2,width = 1,col = "blue", names.arg = c(y[y[,1] < 30,1], ">=30"))

###############################
if(FALSE){
jpeg("task2_try_2.jpg",width = 800)
y30 = sum(y[ y[,1] >= 30,2])
barplot(c(y[y[,1] < 30,2] , y30),space = 0.2,col = "blue", names.arg = c(y[y[,1] < 30,1], ">=30") , main = "histgram of coverage XXX", xlab = "coverage", ylab = "frequency")
dev.off()
	
met30 = sum(y[ y[,1] >= 30,3])	
max = max(y[ y[,1] < 30,3]/y[y[,1] < 30,2] * 100, met30/y30 *100)
barplot(c(y[ y[,1] < 30,3]/y[y[,1] < 30,2] * 100, met30/y30 *100), space = 1,width = 10,col = "blue", names.arg = c(y[y[,1] < 30,1], ">=30") , main = "methylated percentage at each coverage level", ylim=c(0,max + 5) , ylab = "%")
	
}