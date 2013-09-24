
 library("e1071")

 etb = read.table( "_essenGene.baciProfile.091208.csv",sep="\t",header=T);	
 ctb = read.table( "_summary.coat.clusters.I3.1.073108.csv",sep="\t",header=T);

 bacillus.specs = names(ctb)[10:20];
 
 ### generate profile for coat genes
 # ctb2 = cbind( rep("coat", length(ctb[,1])), ctb[, c(1,10:20) ]);
 ctb2 = matrix( nrow=length(ctb[,1]), ncol=11);
 colnames(ctb2) = bacillus.specs;
 row.names(ctb2) = as.character(ctb$X);

 for( j in 1:11){
   for( i in 1:length(ctb2[,1]) ) {
   #for( i in 1:2 ) {
	ctb2[i,j] = ifelse( is.na(ctb[i,j+9]), 0, 1 );
   }
 }

 head(ctb2);
 ctb2 = ctb2[, -2]

 library(RColorBrewer);
 hmcol = colorRampPalette(brewer.pal(10,"RdBu"))(256);
 hmcol = colorRampPalette(brewer.pal(5,"RdBu"))(16);
# hmcol = c("red","blue","green","black","yellow");
 #hmcol = rainbow(5);

 ###plot 

 pdf("091208.heat.coat.pdf", height=8, width=8); 
 heatmap( ctb2, col=hmcol, scale="none", margins = c(5,10) );
 # heatmap( ctb2, col=hmcol, main="coat proteins" ); #does this make sense? YES, default Euclidean distance => hamming distance
 #ctb2["BG10497",]
 #Bsu Bmo Bam Bli Bpu Ban Bce Bth Bwe Bcl Bha 
 # 1   1   1   1   1   0   1   0   0   1   1
 dev.off();

 myhm =  heatmap( ctb2, col=hmcol ); # this give the rows

 pdf("091208.hclust.coat.pdf", height=8, width=8); 

 d.ham = as.dist( hamming.distance( ctb2 ) );
 plot( hclust( d.ham, "ave"), main="hamming distance, average linkage" )

 plot( hclust( dist(ctb2), "ave"), main="Euclidean distance, average linkage" )

 dev.off();

 ### generate profile for essential genes
 etb2 = matrix( nrow=length(etb[,1]), ncol=11);
 colnames(etb2) = bacillus.specs;
 row.names(etb2) = as.character(etb$BG);

 for( j in 1:11) {
   for( i in 1:length(etb[,1]) ) {
    etb2[i,j] = ifelse( is.na(etb[i, j+2]), 0, 1 )
   }
 }

 head(etb2);
 etb2 = etb2[,-2];

 pdf("091208.heat.essen.pdf", height=10, width=10); 
 #heatmap(etb2, );
 heatmap( etb2, col=hmcol, scale="none", labRow=NA );

 dev.off();

quit("yes");

##########
 newColHSVV = function(h1, h2, n) {
    nside = (n/2) - 1;
    s = 1;
    ret = c( hsv(h1, s, seq(1, 1/nside, length=nside), gamma=1),
             "#000000",
             hsv(h2, s, seq(1/nside, 1, length=nside), gamma=1) );
 };
 hmcol = newColHSVV(1/3, 0, 30); 

