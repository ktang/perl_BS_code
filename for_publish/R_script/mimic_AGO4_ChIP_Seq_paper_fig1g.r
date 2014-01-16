# us letter 8.5 X 11 inch
#1   2        3         4       5          6 

# R --slave --vanilla --args   DMR_file png_pre chr < R_script.txt

args = commandArgs()  

file = args[5]
png_pre = args[6]
chr_lable = args[7];

x_label = paste( chr_lable, " (Mb)", sep="")


width_val = 4# 3.8; 
#height_val = 4;
height_val = 3
units_val = "in"; res_val = 100; pointsize_val =8;

file_name = paste(png_pre, "_", chr_lable , "_w", width_val, "_h", height_val, "_res", res_val, "_pt", pointsize_val, ".png", sep="")

gene_distribution_file = "/Users/tang58/Kai_BS/for_publish/9_feature_meth_level_db/TAIR10_all_Gene_number_WinSize100000_sliding100000.txt";
TE_distribution_file = "/Users/tang58/Kai_BS/for_publish/9_feature_meth_level_db/TAIR10_all_TE_number_WinSize100000_sliding100000.txt";

png(file_name, width = width_val, height = height_val, units = units_val , res = res_val, pointsize = pointsize_val );


par(oma = c(3, 3, 0, 3) ) #down left up rigth

data_raw = read.delim(file, head=T, sep ="\t")
data_chr = data_raw[ data_raw[,1] ==  chr_lable  ,]

data_x = data_chr[, 1:3];
data_x[,4] = data_x[,3] - data_x[,2] + 1;

gene_raw =  read.table( gene_distribution_file , head=T, sep ="\t")
TE_raw	 =  read.table( TE_distribution_file , head=T, sep ="\t")

gene_chr = gene_raw [ gene_raw[,1] == chr_lable , ]
TE_chr	 = TE_raw [ TE_raw[,1] == chr_lable , ]

max_x =  trunc( max(data_x[,3] ) / 1000000 ) # ceiling floor trunc round signif

layout(matrix(c(1, 2)), heights=c(2, 1)) # split screen to 2 parts

#par(xpd=F)
par(xpd=NA) #A logical value or NA. If FALSE, 
#all plotting is clipped to the plot region, if TRUE, all plotting is clipped to the figure region, and if NA, all plotting is clipped to the device region

par(mar = c(0, 2, 3, 2) )
plot(  data_x[,2], data_x[,4], type = "p", pch= ".",  xaxt = "n", las = 2, cex = 1.6, ylab = "DMRs length (nts)", xlab = "", font.lab = 2, cex.lab = 1.2);

par(mar = c(2, 2, 0, 2) )

plot( TE_chr[,2], TE_chr[,4] , col = "blue", lwd= 0.8 , las =2, type ="l",col.axis = "blue", xaxt = "n", xlab = x_label, ylab="", font.lab = 2, cex.lab = 2)

#axis(1, col.axis='black',  las = 1 )

seq_array = seq(0, max_x, by = 5)

axis(1, col.axis='black',  las = 1, labels =  seq_array , at =  seq_array * 1000000 , font.axis = 2  ) # bottom

par(new=T)

plot( gene_chr[,2], gene_chr[,4] , yaxt = "n", xaxt = "n", type ="l", col = "red",  xlab = "", ylab="" )
axis(4, col.axis='red',  las =2 )

par(xpd=NA)
mtext("TE/100kb", col="blue", side = 2, line = 2)
mtext("Gene/100kb", col="red", side = 4, line = 2)

#par(mar = mar_ori)

dev.off()
