# us letter 8.5 X 11 inch

# R --slave --vanilla --args  file_pattern(Gene_all TE_all) png_pre  type_TE_Gene dir < /Users/tang58/Kai_BS/useful_R_script/10_gene_TE_bin_files_R_cmd_for_all_command_v0.1.r 
#  /Users/tang58/Kai_BS/useful_R_script/10_gene_TE_bin_files_R_cmd_for_all_command_v0.1.r 
#/Users/tang58/Kai_BS/for_publish/10_gene_TE_bin_files_R_cmd_for_all_command_v0.1.r
# TE_all; Gene_all
args = commandArgs( trailingOnly =  T)  

lwd_line = 2.5;#0.8 # 2


file_pattern = args[1]
png_pre		 = args[2]
type_TE_Gene = args[3]
indir        = args[4]

setwd(indir)

width_val = 8# 3.8; 
#height_val = 2.5;
height_val = 2
units_val = "in"; res_val = 500; pointsize_val =8;


cex_lab_val = 2 #cex.lab =  1.5, 
cex_axis_val = 1.5; #cex.axis = 1.5, 
cex.main = 2

las_val = 1; # axis labels horizontal or vertical

#print("modify file name\n")
# file_name = paste("ibm2Paper_all_Gene", "_w", width_val, "_h", height_val, "_res", res_val, "_pt", pointsize_val, ".png", sep="")
file_name = paste(png_pre , "_w", width_val, "_h", height_val, "_res", res_val, "_pt", pointsize_val, ".png", sep="")


files = dir(pattern= file_pattern)

#print (files)

length = length(files)

x_data = list()

for (i in 1:length(files)){ x_data[[i]] = read.table(files[i], head=T, sep = "\t") }

#x_lables = sub(".*TE_.+_in_(.+)_binNum20_flaking2000bp_binNum20.txt", "\\1", files, perl =T)
#x_lables = sub(".*idm2-1_targets_in_(.+)_binNum20_flaking2000bp_binNum20.txt", "\\1", files, perl =T)
x_lables = sub(".+_in_(.+)_binNum20_flaking2000bp_binNum20.+.txt", "\\1", files, perl =T)

setwd("..")
png(file_name, width = width_val, height = height_val, units = units_val , res = res_val, pointsize = pointsize_val );
#par(mar=c(1,3,1,1)+0.1)
par(oma = c(5, 5, 1, 0) ) #down left up rigth
par(mar= rep(2, 4))


label = c("mCG", "mCHG", "mCHH", "mC")

#cols_db =  c( 
#rgb(115,145,199, maxColorValue = 255), # WT blue
#rgb ( 184,105, 98 ,maxColorValue = 255 ), # ros1-4 red
#rgb ( 177,199, 112 ,maxColorValue = 255 ) # rdd green
#)

#cols_db = c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
#cols_db = c("black", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
cols_db = c(
"gray50",
rgb( 187, 118, 187  , maxColorValue = 255),#asi1, zise
rgb( 140, 255, 119  , maxColorValue = 255),#ibm1, green
"blue",
"black", "darkgreen","blueviolet", "red", "pink", "lawngreen")
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
		 
# ylim = c(0, 1.1*max_y ) ,
 ylim = c(0, 1.02*max_y ) ,
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
#	axis(1, labels = c("-2kb", "TE", "+2kb"), at = c(11,31,51), tick =F, cex.axis = 1.4, font =2) #type_TE_Gene
 	axis(1, labels = c("-2kb", type_TE_Gene, "+2kb"), at = c(11,31,51), tick =F, cex.axis = 1.4, font =2) #type_TE_Gene
	
}

par(xpd=NA)
#legend(-250, -3.5,legend= x_lables , col = cols, lwd = 3, bty = "n",  cex = 2 ,   ncol = length)  #Gene
if( type_TE_Gene == "TE" ){
	legend(-180, -8.5,legend= x_lables , col = cols, lwd = 3, bty = "n",  cex = 0.6 ,   ncol = length)  #TE
}
if( type_TE_Gene == "Gene" ){
	legend(-250, -1.5,legend= x_lables , col = cols, lwd = 3, bty = "n",  cex = 0.6 ,   ncol = length)  #Gene
}
#par(mar = mar_ori)

dev.off()
