# R --slave --vanilla --args  png_pre  dir pattern < XXX.r
args = commandArgs( trailingOnly =  T)  
png_pre		 = args[1]
indir        = args[2]
pattern_input= args[3]

setwd(indir)


width_val = 8 # 3.8; 
height_val = 9;
units_val = "in"; res_val = 500; pointsize_val =8;


#cex_lab_val = 2.5 #cex.lab =  1.5, 
#cex_axis_val = 2; #cex.axis = 1.5, 
#cex.main = 2

las_val = 1; # axis labels horizontal or vertical

file_name = paste( png_pre , "_w", width_val, "_h", height_val, "_res", res_val, "_pt", pointsize_val, ".png", sep="")

png(file_name, width = width_val, height = height_val, units = units_val , res = res_val, pointsize = pointsize_val );

#par(mar=c(1,3,1,1)+0.1)

par(oma = c(5, 5, 1, 0) ) #down left up rigth
par(mar= rep(2, 4))

mar_ori = par("mar")

#par(mar = c(5, 6 , 4, 2) + 0.1)

#print("modify label TE/GENE")

files = dir(pattern= pattern_input )

#print (files)

length = length(files)

x_data = list()

for (i in 1:length(files)){ x_data[[i]] = read.delim(files[i], head=T, sep = "\t") }

#x_lables = sub(".*TE_.+_in_(.+)_binNum20_flaking2000bp_binNum20.txt", "\\1", files, perl =T)
#x_lables = sub(".*idm2-1_targets_in_(.+)_binNum20_flaking2000bp_binNum20.txt", "\\1", files, perl =T)

# x_lables = sub(".+_in_(.+)_binNum20_flaking2000bp_binNum20_dep..txt", "\\1", files, perl =T)
# ibm1-4A_vs_colA_hyper

# x_lables = sub("(.+)_vs_colA_.+_P0.05_reduced_boundary_All_dep4_WinSize200_sliding50_gap100_iniCut2_repCut5_edm2Paper_detail_EckerPaper_annotated_TAIR10_len100_wCfold2_fCdiff10.txt", "\\1", files, perl =T)
 x_lables = sub("(.+)_vs_colA_.+", "\\1", files, perl =T)

chr_len = c(30427671, 19698289, 23459830, 18585056, 26975502 );


cols_db =  c( 

# rgb ( 145,145,145 ,maxColorValue = 255 ), #col-0 Gray

rgb ( 222,147,142 ,maxColorValue = 255 ), #edm2 red
#rgb (200,226, 126  ,maxColorValue = 255 ), # ibm1 green
rgb (177,199,112  ,maxColorValue = 255 ), # ibm1 green for ibm2paper
rgb ( 142,120,176  ,maxColorValue = 255 ), # ibm2paper ibm2 purple

#rgb ( 160, 132, 200  ,maxColorValue = 255 ) # rdd purple for edm2Paper
rgb ( 127,184,207 ,maxColorValue = 255 ) # rdd blue for edm2Paper

#rgb( 0, 255, 255, maxColorValue = 255)   #blue

)


plot.new()
plot.window(xlim=c(0,6),ylim=c(-31000000,500000))

lines(c(1,1),c(0,-30427671),lwd = 10)
lines(c(2,2),c(0,-19698289),lwd = 10)
lines(c(3,3),c(0,-23459830),lwd = 10)
lines(c(4,4),c(0,-18585056),lwd = 10)
lines(c(5,5),c(0,-26975502),lwd = 10)

#lines(c(1.2,1.2),c(0,-30427671))
#lines(c(2.2,2.2),c(0,-19698289))
#lines(c(3.2,3.2),c(0,-23459830))
#lines(c(4.2,4.2),c(0,-18585056))
#lines(c(5.2,5.2),c(0,-26975502))

symbols(c(1,2,3,4,5),c(-15000000,-3600000,-13500000,-3920000,-11700000),circles=c(1,1,1,1,1),add=T,inches=0.08,fg="red",bg="red")

thin_lwd = 0.3

for (i in 1:length(files) ){
#	data_individual = x_data[[i]][,1:2];
	data_individual = x_data[[i]][,1:2];
	for (j in 1:5){
		chr = paste("chr", j, sep="");
		data = data_individual [data_individual[, 1] == chr, ];
		
		gap =  j+ (i-1) * 0.2 + 0.15;
		
		lines( c( gap, gap ), c(0, -chr_len[j]) ,  col = cols_db[i]);
		data[,3] =  rep (gap, dim(data)[1] );
		data[,4] = data[,3] + 0.08;
		for(k in 1:dim(data)[1]){
			lines( c( data[k,3], data[k,4]) , c( -data[k,2], -data[k,2]  ) , lwd = thin_lwd , col = cols_db[i] )
	
		}
	}
	
}

par(xpd=NA)

text_y = 700000;

text( c(1,2,3,4,5), rep(text_y, 5), # c(5000,5000,5000,5000,5000), 
  labels = c( "chr1","chr2","chr3","chr4","chr5"  )  ,pos=3, cex = 2);

#par(xpd=NA)
legend(2, -28000000,legend= x_lables , col = cols_db, lwd = 3, bty = "n",  cex = 1.2
,   ncol = length
)

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
