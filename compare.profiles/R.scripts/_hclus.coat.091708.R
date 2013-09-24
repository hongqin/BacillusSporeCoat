
 rm( list = ls() );
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
 ctb2 = ctb2[, -c(2,10,11)]

 #d.ham = as.dist( hamming.distance( ctb2 ) ); 
 #hm, hamming distance and Euclidean distance is different for 3 genes? 

pdf("_091708.iteration.hclust.coat.pdf", height=10, width=10); 

 mymethods = c("ward","complete","average","single","median","centroid", "mcquitty");
 spec.colors = c("cyan","cyan","cyan","cyan","orange","orange","orange","orange" );
 names(spec.colors) = c( "Bsu", "Bam", "Bli", "Bpu", "Ban", "Bce", "Bth", "Bwe")

for( mymethod in mymethods ) {
 hd =  hclust( dist(ctb2), mymethod); 
 # plot( hd, main="hamming distance, ward linkage" )
 coat.cat = cutree(hd, 4 )
 col.palette = c("red","brown","blue","green");
 coat.color = col.palette[coat.cat]

 library(RColorBrewer);
 #hmcol = colorRampPalette(brewer.pal(10,"RdBu"))(256);
 hmcol = colorRampPalette(brewer.pal(5,"RdBu"))(16);

 #heatmap( ctb2, col=hmcol, scale="none", margins = c(5,10) );
 heatmap( ctb2, col=hmcol, scale="none", margins = c(5,10), 
 RowSideColors=coat.color, ColSideColors = spec.colors,
 hclustfun = function(c) hclust( c, method=mymethod),
 #distfun = function(c) as.dist(hamming.distance(c)) #Hamming is less pleasant than Euclidean 
 main = mymethod
 );
}
dev.off();

##ward plot
pdf("_091708.ward.hclust.coat.pdf", height=10, width=10); 
 mymethod = "ward";
 hd =  hclust( dist(ctb2), mymethod); 
 coat.cat = cutree(hd, 4 )
 col.palette = c("red","brown","blue","green");
 coat.color = col.palette[coat.cat]

 library(RColorBrewer);
 hmcol = colorRampPalette(brewer.pal(5,"RdBu"))(16);

 heatmap( ctb2, col=hmcol, scale="none", margins = c(5,10), 
 RowSideColors=coat.color, ColSideColors = spec.colors,
 hclustfun = function(c) hclust( c, method=mymethod) );
dev.off();

#quit("yes");

 cats = c("Fluid","Subtilis","Cereus","Conserved");
 coat.cat2 = cats[coat.cat]

 #plot( hclust( dist(ctb2), "ave"), main="Euclidean distance, average linkage" )

##############################
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
 etb2 = etb2[,-c(2,10,11)];

 pdf("_091708.heat.essen.pdf", height=10, width=10); 

 mymethod = "ward";
 hd =  hclust( dist(etb2), mymethod); 
 essen.cat = cutree(hd, 4 )
 col.palette = c("red","brown","blue","green");
 essen.color = col.palette[essen.cat]

 library(RColorBrewer);
 hmcol = colorRampPalette(brewer.pal(5,"RdBu"))(16);

 heatmap( etb2, col=hmcol, scale="none", margins = c(5,10), 
 RowSideColors=essen.color, ColSideColors = spec.colors,
 hclustfun = function(c) hclust( c, method=mymethod),
 labRow = NA );

 #heatmap(etb2, );
 #heatmap( etb2, col=hmcol, scale="none", labRow=NA );

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

