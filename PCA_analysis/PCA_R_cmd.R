# us letter 8.5 X 11 inch

# R --slave --vanilla --args  file_pattern png_pre  type_TE_Gene dir < /Users/tang58/Kai_BS/for_publish/10_gene_TE_bin_files_R_cmd_for_all_command_v0.1.r
# TE_all; Gene_all
args = commandArgs( trailingOnly =  T)  

lwd_line = 1;#0.8 # 2


file_pattern = args[1]
png_pre		 = args[2]
type_TE_Gene = args[3]
indir        = args[4]

setwd(indir)

width_val = 8# 3.8; 
#height_val = 2.5;
height_val = 2
units_val = "in"; res_val = 500; pointsize_val =8;

files = dir(pattern= file_pattern)
#setwd("/Volumes/Macintosh_HD_2/PCA_analysis/5kb/output/CG")
i = 1
tmp = read.delim (files[i])
matrix_data = matrix(data = tmp[,2],nrow= dim(tmp)[1],ncol=1)
col_lable = sub("(.+)_mC.+_level", "\\1", colnames(tmp)[2], perl =T)

colnames(matrix_data)[i] = col_lable

for (i in 2:length(files)){
	tmp = read.delim (files[i])
	matrix_data = cbind(matrix_data, tmp[,2])
	
	col_lable = sub("(.+)_mC.+_level", "\\1", colnames(tmp)[2], perl =T)
	colnames(matrix_data)[i] = col_lable

}

m_na = na.omit(matrix_data)

mt = t(m_na)
mt.pca = prcomp(mt)

pdf("CG_debug.pdf"); biplot(mt.pca,var.axes=F, cex = c(0.4, 0.00000000000000000001));dev.off()


x_plot = mt.pca$x[,1:2]


library(calibrate)

ext = 200
plot(x_plot[,1],x_plot[,2],pch=".", cex = 5, col="red", xlim = c(min(x_plot[,1]) - ext, max(x_plot[,1])+ ext))
textxy(x_plot[,1], x_plot[,2], rownames(x_plot), cx = 0.6)


plot(x_plot[,1],x_plot[,2],pch=".", cex = 5, col="red", xlim = c(min(x_plot[,1]) - ext, max(x_plot[,1])+ ext))
textxy(x_plot[,1], x_plot[,2], rownames(x_plot), cx = 0.6, m=c( mean(x_plot[,1]) , mean(x_plot[,2]) )  )