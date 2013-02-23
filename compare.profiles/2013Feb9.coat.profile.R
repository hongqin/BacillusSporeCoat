# 2013 Feb 9, based on 111608.hclust.cuttree.R
# re-adjust col. use red for labil and blue for conserved.
# re-adjust Col order and Row order
# cuttree(hd, 4)  for coat, 
# cuttree(hd, 3) for essential genes

 rm( list = ls() );
 library("e1071")

# numclus =4 ; 
 numclus =2 ; #111808 change

 etb = read.table( "_essenGene.baciProfile.091208.csv",sep="\t",header=T);	
 ctb = read.csv( "_coat.profile.110808.csv");

 bacillus.specs = names(ctb)[10:20];
 
 ### generate profile for coat genes
 # ctb2 = cbind( rep("coat", length(ctb[,1])), ctb[, c(1,10:20) ]);
 ctb2 = matrix( nrow=length(ctb[,1]), ncol=11);
 colnames(ctb2) = bacillus.specs;
 rownames(ctb2) = as.character( paste(ctb$id, ctb$Gene) );

 for( j in 1:11){
   for( i in 1:length(ctb2[,1]) ) {
   #for( i in 1:2 ) {
	ctb2[i,j] = ifelse( is.na(ctb[i,j+9]), 0, 1 );
   }
 }

 head(ctb2);
 ctb2 = ctb2[, -c(2,10,11)]

 #check
 head(ctb[,11:20]);
 head(ctb2);
# ctb2["BG11380",]

 #d.ham = as.dist( hamming.distance( ctb2 ) ); 
 #hm, hamming distance and Euclidean distance is different for 3 genes? 

pdf("_111808.iteration.hclust.2.coat.pdf", height=10, width=10); 

 mymethods = c("ward","complete","average","single","median","centroid", "mcquitty");
 spec.colors = c("cyan","cyan","cyan","cyan","orange","orange","orange","orange" );
 names(spec.colors) = c( "Bsu", "Bam", "Bli", "Bpu", "Ban", "Bce", "Bth", "Bwe")

for( mymethod in mymethods ) {
 hd =  hclust( dist(ctb2), mymethod); 
 # plot( hd, main="hamming distance, ward linkage" )
 coat.cat = cutree(hd, numclus )  ###<=== change is here
 col.palette = c("blue", "red","brown","green");
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

# ctb2 = ctb2[,c( "Bpu", "Bli", "Bam",  "Bsu", "Bwe", "Bce", "Ban","Bth") ]; 

#pdf("_111808.average.hclust.2.coat.pdf", height=17, width=17); 
pdf("_20130209.average.hclust.2.coat.pdf", height=5, width=5); 
 library(RColorBrewer);
 hmcol = colorRampPalette(brewer.pal(5,"RdBu"))(16);
 
 mymethod = "average";
 hd =  hclust( dist(ctb2), mymethod); 
 coat.cat = cutree(hd, numclus )
 #col.palette = c("blue","red"); #2013 Feb 9 change
 col.palette = hmcol[c(16,1)]; #2013 Feb 17 change
 coat.color = col.palette[coat.cat]

 hcol =  hclust(dist( t(ctb2)), method=mymethod )
 #hcol =  as.dendrogram( hclust(dist( t(ctb2)), method=mymethod ) )

 heatmap( ctb2, col=hmcol, scale="none", margins = c(5,22),  
 #heatmap( ctb2, col=c("red","blue"), scale="none", margins = c(5,22),  #2013Feb9 change
 #Colv = NA,  ## no dedrogram
 #Colv = FALSE, ## use the input column order,  
 #Colv = as.dendrogram( hclust(dist( t(ctb2)), method=mymethod ) ),
 RowSideColors=coat.color, # ColSideColors = spec.colors,  #2013Feb9, remove spec.colors
 hclustfun = function(c) hclust( c, method=mymethod), 
 cexRow=0.4, cexCol=0.6 );

# grid() #not working
 # require(Heatplus); #no grid line
# heatmap_2(ctb2, do.dendro = c(TRUE, TRUE), legend = 0, scale="none",RowSideColors=coat.color,
#           legfrac = 10,
#           hclustfun = function(c) hclust( c, method= mymethod),
#           col = hmcol
# )
 
dev.off();
#colors()[ grep("blue", colors()) ]
 
 # cats = c("Fluid","Subtilis","Cereus","Conserved");
 # coat.cat2 = cats[coat.cat]

##############################
 ### generate profile for essential genes
 numclus = 3 ; 

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

 etb2 = etb2[,c( "Bwe", "Bce", "Ban","Bth","Bpu", "Bli", "Bam",  "Bsu") ]; # fix the order

 # spec.colors = c( "cyan","cyan","cyan","cyan","orange","orange","orange","orange" );
 spec.colors = c( "orange","orange","orange","orange" , "cyan","cyan","cyan","cyan");
 names(spec.colors) = c( "Bsu", "Bam", "Bli", "Bpu", "Ban", "Bce", "Bth", "Bwe")
 # There seems to be a bug here for ColSideColors

 pdf("_111808.average.hclus.essen.pdf", height=10, width=10); 

 mymethod = "ward";
 mymethod = "average"; #111808 
 hd =  hclust( dist(etb2), mymethod); 
 essen.cat = cutree(hd, numclus )
 #col.palette = c("blue","green","red","brown",);
 col.palette = c("blue","red","green" );
 essen.color = col.palette[essen.cat]

 library(RColorBrewer);
 hmcol = colorRampPalette(brewer.pal(5,"RdBu"))(16);

 heatmap( etb2, col=hmcol, scale="none", margins = c(5,10), 
  Colv = NA, ##no dendrogram
  RowSideColors=essen.color, ColSideColors = spec.colors,
  hclustfun = function(c) hclust( c, method=mymethod),
  labRow = NA);

 #heatmap(etb2, );
 #heatmap( etb2, col=hmcol, scale="none", labRow=NA );

 dev.off();

quit("yes");

