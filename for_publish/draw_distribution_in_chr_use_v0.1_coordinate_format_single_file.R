# R --slave --vanilla --args   input_list  png_pre < /Users/tang58/Kai_BS/for_publish/draw_distribution_in_chr_use_v0.1_coordinate_format_single_file.R
#   setwd("/Users/tang58/misc/Huiming_Zhang/13_Feb15_rrp6L1_smallRNA/list_exp50_fold2_Mar5_2013")
# read me:
# input is 1:1-3 format and has a header

args = commandArgs( trailingOnly =  T)  
input		 = args[1]
png_pre		 = args[2]
#indir        = args[2]
#pattern_input= args[3]
#setwd(indir)


data_individual = read.table( input, header = T )
#  data_individual = read.table("344_exp50_fold2_loci.txt",head=T)

width_val = 8 # 3.8; 
height_val = 9;
units_val = "in"; res_val = 500; pointsize_val =8;
file_name = paste( png_pre , "_w", width_val, "_h", height_val, "_res", res_val, "_pt", pointsize_val, ".png", sep="")
png(file_name, width = width_val, height = height_val, units = units_val , res = res_val, pointsize = pointsize_val );


#cex_lab_val = 2.5 #cex.lab =  1.5, 
#cex_axis_val = 2; #cex.axis = 1.5, 
#cex.main = 2

las_val = 1; # axis labels horizontal or vertical



#par(mar=c(1,3,1,1)+0.1)

par(oma = c(5, 5, 1, 0) ) #down left up rigth
par(mar= rep(2, 4))

mar_ori = par("mar")

#par(mar = c(5, 6 , 4, 2) + 0.1)

chr_len = c(30427671, 19698289, 23459830, 18585056, 26975502 );

plot.new()
plot.window(xlim=c(0,6),ylim=c(-31000000,500000))

lines(c(1,1),c(0,-30427671),lwd = 10)
lines(c(2,2),c(0,-19698289),lwd = 10)
lines(c(3,3),c(0,-23459830),lwd = 10)
lines(c(4,4),c(0,-18585056),lwd = 10)
lines(c(5,5),c(0,-26975502),lwd = 10)
lines(c(1.2,1.2),c(0,-30427671))
lines(c(2.2,2.2),c(0,-19698289))
lines(c(3.2,3.2),c(0,-23459830))
lines(c(4.2,4.2),c(0,-18585056))
lines(c(5.2,5.2),c(0,-26975502))
symbols(c(1,2,3,4,5),c(-15000000,-3600000,-13500000,-3920000,-11700000),circles=c(1,1,1,1,1),add=T,inches=0.08,fg="red",bg="red")

thin_lwd = 0.5
#x= cbind ( data_individual[,c(1,1,1)] ) 
x = cbind (data_individual[,c(1,1) ], 1:dim(data_individual)[1], 1:dim(data_individual)[1] );

for( i in 1:dim(data_individual)[1] ){
		temp = as.character( x[i,1]); 
		a1 = unlist( strsplit( temp, ":"  )); 
		a2 = unlist( strsplit( a1[2], "-" ));
		chr = as.numeric( a1[1] )
		start = as.numeric( a2[1] )
		x[i,3] = chr
		x[i,4] = start
}

#data_individual = x
for (i in 1:5){
#	chr = paste("chr", j, sep="");
	chr = i
	data = x [x[, 3] == chr, ];
	for( k in 1:dim(data)[1]){
		lines( c( i + 0.2, i+0.4) , c( -data[k,4], -data[k,4]), lwd = thin_lwd, col = "black" )
	}
}


par(xpd=NA)

text_y = 700000;

text( c(1,2,3,4,5), rep(text_y, 5), # c(5000,5000,5000,5000,5000), 
  labels = c( "1","2","3","4","5"  )  ,pos=3, cex = 2);

#par(xpd=NA)
#legend(2, -28000000,legend= x_lables , col = cols_db, lwd = 3, bty = "n",  cex = 1.2
#,   ncol = length
#)

dev.off()

if(FALSE) {
for(i in 1:dim(data)[1])
{
	if(data[i,1] =="chr1")
	{lines(c(1.2,1.4),c(-data[i,2],-data[i,2]),lwd = thin_lwd)}
	
	if(data[i,1] =="chr2")
	{lines(c(2.2,2.4),c(-data[i,2],-data[i,2]),lwd = thin_lwd)}
	
	if(data[i,1] =="chr3")
	{lines(c(3.2,3.4),c(-data[i,2],-data[i,2]),lwd = thin_lwd)}
	
	if(data[i,1] =="chr4")
	{lines(c(4.2,4.4),c(-data[i,2],-data[i,2]),lwd = thin_lwd)}
	
	if(data[i,1] =="chr5")
	{lines(c(5.2,5.4),c(-data[i,2],-data[i,2]),lwd = thin_lwd)}
	
}
} #FALSE
