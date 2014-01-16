# setwd("/Users/tang58/misc/Mingguang_Lei/paper/edm2/paper_figures/TE_Gene_meth_level_dist/TE")

# source("~/Kai_BS/for_publish/10_gene_TE_bin_files_R_cmd_TE.r")

# us letter 8.5 X 11 inch
width_val = 8# 3.8; 
height_val = 5;
units_val = "in"; res_val = 100; pointsize_val =8;

cex_lab_val = 1.8 #cex.lab =  1.5, 
# cex_axis_val = 2; #cex.axis = 1.5, 
cex_main = 1.8

las_val = 1; # axis labels horizontal or vertical

print("modify file name\n")
file_name = paste("ibm2Paper_TE_length", "_w", width_val, "_h", height_val, "_res", res_val, "_pt", pointsize_val, ".png", sep="")

png(file_name, width = width_val, height = height_val, units = units_val , res = res_val, pointsize = pointsize_val );

#par(mar=c(1,3,1,1)+0.1)

#TE:

#All_TE
#0_0.5kb
#0.5_1kb
#1_2kb
#2_4kb
#4kb_up

length_labels = c( "0_0.5kb", "0.5_1kb","1_2kb", "2_4kb", "4kb_up" )
length_labels_in_fig = c( "<0.5kb", "0.5-1kb","1-2kb", "2-4kb", ">4kb" )

#length_labels = c( "0_1kb", "1_2kb","2_3kb", "3_4kb", "4kb_up" )
#length_labels_in_fig = c( "<1kb", "1-2kb","2-3kb", "3-4kb", ">4kb" )

y_labs = c( "mCG level (%)", "mCHG level (%)", "mCHH level (%)" );


files_list = list();

for (i in 1:length( length_labels ) ){
	files_list[[i]] = dir(pattern = length_labels[i] )
}




#print (files)

lwd_line = 2

#label = c("CG", "CHG", "CHH", "C")
label = c("mCG", "mCHG", "mCHH")

#ago4 double ago6 wt

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

max_CG = vector();
max_CHG = vector();
max_CHH = vector();

layoutmat =  matrix(1: (3 * length( length_labels )), nrow = 3) ;layout(layoutmat)
#> matrix(1:15,nrow = 3)
#[,1] [,2] [,3] [,4] [,5]
#[1,]    1    4    7   10   13
#[2,]    2    5    8   11   14
#[3,]    3    6    9   12   15

par(oma = c(6, 6, 6, 0.5) ) #down left up rigth
#par(mar= rep(0.8, 4))
mar_ori = par("mar")


for (i_file_list in 1:length(files_list ) ){
	files = files_list[[ i_file_list ]];
#	length = length(files)
	x_data = list()
	for (i in 1:length(files)){ x_data[[i]] = read.table(files[i], head=T, sep = "\t") }
	
	max_y_CG = max (unlist(lapply(x_data, function(x) max(x[,3* 1 +1 ]) ) ) ) 
	max_y_CHG = max (unlist(lapply(x_data, function(x) max(x[,3* 2+1 ]) ) ) ) 
	max_y_CHH = max (unlist(lapply(x_data, function(x) max(x[,3* 3+1 ]) ) ) ) 
	
	max_CG = c(max_CG, max_y_CG );
	max_CHG = c(max_CHG, max_y_CHG );
	max_CHH = c(max_CHH, max_y_CHH );
}

max_y_vector = c( max(max_CG), max(max_CHG), max(max_CHH) )

for (i_file_list in 1:length(files_list ) ){
	files = files_list[[ i_file_list ]];
	
	length = length(files)
	x_data = list()
	for (i in 1:length(files)){ x_data[[i]] = read.table(files[i], head=T, sep = "\t") }
	
	x_lables = sub("TE.+_in_(.+)_binNum20_flaking2000bp_binNum20.+.txt", "\\1", files, perl =T);
	
	cols = cols_db[1:length]
	
		par(xpd=F)
		for (i in 1:3){
		
#	max_y = max (unlist(lapply(x_data, function(x) max(x[,3*i+1 ]) ) ) )  #
			max_y = max_y_vector[i];
			y_lab = "";
			if(i_file_list == 1){
#y_lab = "% Methylation ";
				y_lab = y_labs[i]
				par(xpd=NA)
			}else{
				y_lab = "";
			}
			
			if( i == 1){
				main_label = length_labels_in_fig[i_file_list]
				par(mar= c(0.6, 0.6, 2, 0.6) )  #down left up rigth

			}else{
				main_label = ""
				par(mar= rep(0.6, 4))

			}
			
			if( i_file_list == 1 ){cex_axis_val = 1.5}
			else{ cex_axis_val =  0.0000001}
				
	
#plot
			plot(x_data[[1]][,1], x_data[[1]][,3*i+1], type="l", 
			
				 col = cols[1], 
			 
				 ylim = c(0, 1.1*max_y ) ,
				 ylab = y_lab, 
			 
#	 main = label[i], 
				 main = main_label,
				 xaxt="n", 
				 xlab="", 
				 cex.main = cex_main, 
			 
				 cex.lab =  cex_lab_val ,
				 font.lab = 2 ,# bold
				 
				 cex.axis = cex_axis_val,
			 
				 lwd = lwd_line,
				 las = las_val
				 
			)
#			plot
			
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
			if(i == 3){
				axis(1, labels = c("-2kb", "TE", "+2kb"), at = c(11,31,51), tick =F, cex.axis = 1.4, font =2)
			}
	}
}

par(xpd=NA)
################
legend(-280, -2.5,legend= x_lables , col = cols, lwd = 3, bty = "n",  cex = 2
,   ncol = length,

) # ncol = length)
################



par(mar = mar_ori)



dev.off()
