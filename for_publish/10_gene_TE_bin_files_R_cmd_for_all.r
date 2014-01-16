#source("~/Kai_BS/for_publish/10_gene_TE_bin_files_R_cmd_for_all.r")
# us letter 8.5 X 11 inch

width_val = 8# 3.8; 
height_val = 2.5;
units_val = "in"; res_val = 100; pointsize_val =8;


cex_lab_val = 2.5 #cex.lab =  1.5, 
cex_axis_val = 2; #cex.axis = 1.5, 
cex.main = 2

las_val = 1; # axis labels horizontal or vertical

print("modify file name\n")
file_name = paste("ibm2Paper_all_Gene", "_w", width_val, "_h", height_val, "_res", res_val, "_pt", pointsize_val, ".png", sep="")


png(file_name, width = width_val, height = height_val, units = units_val , res = res_val, pointsize = pointsize_val );

#par(mar=c(1,3,1,1)+0.1)

par(oma = c(5, 5, 1, 0) ) #down left up rigth
par(mar= rep(2, 4))



mar_ori = par("mar")

#par(mar = c(5, 6 , 4, 2) + 0.1)

print("modify label TE/GENE")

files = dir(pattern=".*in.+txt")

#print (files)

length = length(files)

x_data = list()

for (i in 1:length(files)){ x_data[[i]] = read.table(files[i], head=T, sep = "\t") }

#x_lables = sub(".*TE_.+_in_(.+)_binNum20_flaking2000bp_binNum20.txt", "\\1", files, perl =T)
#x_lables = sub(".*idm2-1_targets_in_(.+)_binNum20_flaking2000bp_binNum20.txt", "\\1", files, perl =T)
x_lables = sub(".+_in_(.+)_binNum20_flaking2000bp_binNum20.+.txt", "\\1", files, perl =T)


lwd_line = 2

label = c("mCG", "mCHG", "mCHH", "mC")

#head
if(FALSE){
	cols_db =  c( 
				 "blue",
#rgb( 0, 255, 255, maxColorValue = 255),   #blue
				 rgb( 255,0 , 255, maxColorValue = 255),  #pink
				 
#rgb( 181, 204, 82,  maxColorValue = 255),  #qing zi se
#rgb( 242, 221, 51,  maxColorValue = 255),  #yellow
				 "darkgreen",
				 
				 rgb(153, 153, 153, maxColorValue = 255), #grey
				 
				 
#rgb( 102, 255, 204, maxColorValue = 255), #light blue
				 rgb( 181, 204, 82,  maxColorValue = 255),  #qing zi se
				 
				 rgb( 255, 102,0, maxColorValue = 255 ) #orange
#rgb( 255,0 , 255, maxColorValue = 255),  #pink
#rgb( 0, 255, 255, maxColorValue = 255)   #blue
				 
				 )
	
}
#cut

cols_db =  c( 
rgb(153, 153, 153, maxColorValue = 255), #grey
rgb ( 222,147,142 ,maxColorValue = 255 ), #edm2 red

#rgb (200,226, 126  ,maxColorValue = 255 ), # ibm1 green
rgb (177,199,112  ,maxColorValue = 255 ), # ibm1 green for ibm2paper


#rgb ( 160, 132, 200  ,maxColorValue = 255 ) # rdd purple

rgb ( 142,120,176  ,maxColorValue = 255 ), # ibm2paper ibm2 purple

#rgb ( 160, 132, 200  ,maxColorValue = 255 ) # rdd purple for edm2Paper
rgb ( 127,184,207 ,maxColorValue = 255 ) # rdd blue for edm2Paper


#rgb( 0, 255, 255, maxColorValue = 255)   #blue

)
cols = cols_db[1:length]
#cols = c("black", "blue", "darkgreen","blueviolet", "red", "pink", "lawngreen")
layoutmat =  matrix(1:4,ncol = 4) ;layout(layoutmat)

par(xpd=F)
for (i in 1:4){
# max_y = max( pro_WT[, 3*i+1 ],  pro_JKZ3[, 3*i+1],  pro_JKZ4[, 3*i+1],  pro_ros1[, 3*i+1], pro_newJKZ4[, 3*i+1] )
	max_y = max (unlist(lapply(x_data, function(x) max(x[,3*i+1 ]) ) ) )  #
	
	y_lab = "";
	if(i == 1){
		y_lab = "% Methylation ";
		par(xpd=NA)

	}else{
		y_lab = "";
	}
	
	plot(x_data[[1]][,1], x_data[[1]][,3*i+1], type="l", 
		 col = cols[1], 
		 
		 ylim = c(0, 1.1*max_y ) ,
		 ylab = y_lab, 

		 main = label[i], 
		 xaxt="n", 
		 xlab="", 

    	 cex.lab =  cex_lab_val ,
		 cex.main = 2, 
		 cex.axis = cex_axis_val,

		 lwd = lwd_line,
		 las = las_val
		 )
	par(xpd=F)
	for(j in 2:length){
			lines(  x_data[[j]][,1], x_data[[j]][,3*i+1],
					col = cols[j],
				    pch = j, 
				    type = "l", 
				    cex = 0.5 , 
				    lwd = lwd_line,
				    las = las_val, 
				    cex.axis =  cex_axis_val
				  )
	}
	
	abline(v = 21, lty = 2); 	abline(v = 40, lty = 2)
	axis(1, labels = c("-2kb", "Gene", "+2kb"), at = c(11,31,51), tick =F, cex.axis = 1.4, font =2)
#	axis(1, labels = c("-2kb", "TE", "+2kb"), at = c(11,31,51), tick =F, cex.axis = 1.4, font =2)
}

par(xpd=NA)
legend(-250, -0.8,legend= x_lables , col = cols, lwd = 3, bty = "n",  cex = 2 ,   ncol = length)  #Gene
#legend(-250, -3.5,legend= x_lables , col = cols, lwd = 3, bty = "n",  cex = 2 ,   ncol = length)  #Gene

par(mar = mar_ori)

dev.off()
