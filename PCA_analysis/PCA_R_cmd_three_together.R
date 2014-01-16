# us letter 8.5 X 11 inch

# R --slave --vanilla --args  dir pre  < /Users/tang58/Kai_BS/PCA_analysis/PCA_R_cmd_three_together.R
args = commandArgs( trailingOnly =  T)  

lwd_line = 1;#0.8 # 2

ext = 200
cx_val = 0.6

indir = args[1]
png_pre = args[2]


setwd(indir)
library(calibrate)

if( file.exists("CG") &  file.exists("CHG") &  file.exists("CHH") ){
	
}else{
	print("wrong")
	q("no")
}



width_val = 5# 3.8; 
#height_val = 2.5;
height_val = 5
units_val = "in"; res_val = 500; pointsize_val =8;

types = c("CG", "CHG", "CHH");

for (type in types){
#	print (type)
	setwd(type)
	files = dir(pattern= "dep4.txt")
#	print (files);
	
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
	
	print (dim(matrix_data))
	
	m_na = na.omit(matrix_data)
	print (dim(m_na))

	mt = t(m_na)
	mt.pca = prcomp(mt)
	x_plot = mt.pca$x[,1:2]

	
	setwd("..")
	file_name = paste(png_pre , "_", type,  "_w", width_val, "_h", height_val, "_res", res_val, "_pt", pointsize_val, "_ext",ext , "_cx_val",cx_val ,  ".png", sep="")
	png(file_name, width = width_val, height = height_val, units = units_val , res = res_val, pointsize = pointsize_val );
#	print(file_name)
	
	plot(x_plot[,1],x_plot[,2],pch=".", cex = 5, col="red", xlim = c(min(x_plot[,1]) - ext, max(x_plot[,1])+ ext), main=type, xlab= "PC1",ylab="PC2")
	textxy(x_plot[,1], x_plot[,2], rownames(x_plot), cx =cx_val, m=c( mean(x_plot[,1]) , mean(x_plot[,2]) )  )

	dev.off()
}



if(FALSE){
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

#pdf("CG_debug.pdf"); biplot(mt.pca,var.axes=F, cex = c(0.4, 0.00000000000000000001));dev.off()


x_plot = mt.pca$x[,1:2]


library(calibrate)

ext = 200
plot(x_plot[,1],x_plot[,2],pch=".", cex = 5, col="red", xlim = c(min(x_plot[,1]) - ext, max(x_plot[,1])+ ext))
textxy(x_plot[,1], x_plot[,2], rownames(x_plot), cx = 0.6)


plot(x_plot[,1],x_plot[,2],pch=".", cex = 5, col="red", xlim = c(min(x_plot[,1]) - ext, max(x_plot[,1])+ ext))
textxy(x_plot[,1], x_plot[,2], rownames(x_plot), cx = 0.6, m=c( mean(x_plot[,1]) , mean(x_plot[,2]) )  )

}