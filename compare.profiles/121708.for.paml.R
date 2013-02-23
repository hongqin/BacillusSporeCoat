##121708  parse genes for PAML

 rm( list = ls() );
 library("e1071")

 etb = read.table( "_essenGene.baciProfile.091208.csv",sep="\t",header=T);	
 ctb = read.csv( "_coat.profile.110808.csv");

 bacillus.specs = names(ctb)[10:20];
 
 ### generate profile for coat genes
 ctb2 = matrix( nrow=length(ctb[,1]), ncol=11);
 colnames(ctb2) = bacillus.specs;
 rownames(ctb2) = as.character( paste(ctb$id, ctb$Gene) );

 for( j in 1:11){
   for( i in 1:length(ctb2[,1]) ) {
	ctb2[i,j] = ifelse( is.na(ctb[i,j+9]), 0, 1 );
   }
 }

 head(ctb2);
 ctb2 = ctb2[, -c(2,10,11)]

 #check
 head(ctb[,11:20]);
 head(ctb2);

 #mymethods = c("ward","complete","average","single","median","centroid", "mcquitty");
 spec.colors = c("cyan","cyan","cyan","cyan","orange","orange","orange","orange" );
 names(spec.colors) = c( "Bsu", "Bam", "Bli", "Bpu", "Ban", "Bce", "Bth", "Bwe")

 mymethod = "single";
 hd =  hclust( dist(ctb2), method=mymethod); 
 coat.cat = cutree(hd, h=0 ); ##121708 This is the classificaiton results

 col.palette = rainbow( length(table(coat.cat)) );
 coat.color = col.palette[coat.cat]

 library(RColorBrewer);
 hmcol = colorRampPalette(brewer.pal(5,"RdBu"))(16);

 heatmap( ctb2, scale="none", margins = c(5,22),  
 #col=RGBColVec(64),
 col=hmcol, 
 #Colv = NA,  ## no dedrogram
 #Colv = FALSE, ## use the input column order,  
 #Colv = as.dendrogram( hclust(dist( t(ctb2)), method=mymethod ) ),
 #RowSideColors=coat.color, #ColSideColors = spec.colors,
 hclustfun = function(c) hclust( c, method=mymethod), 
 cexRow=0.4, cexCol=0.6 );

 #retrieve the largest cluster, the 4sub4cereus groups  
 conserved.coat = coat.cat[coat.cat==3]
 write.table( conserved.coat, "_4subtilis4cereus.orth.group.121708.tab", quote=F,sep="\t", col.names=F);


#q("yes");


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

 etb2 = etb2[,c( "Bwe", "Bce", "Ban","Bth","Bpu", "Bli", "Bam",  "Bsu") ]; # fix the order

 # spec.colors = c( "cyan","cyan","cyan","cyan","orange","orange","orange","orange" );
 spec.colors = c( "orange","orange","orange","orange" , "cyan","cyan","cyan","cyan");
 names(spec.colors) = c( "Bsu", "Bam", "Bli", "Bpu", "Ban", "Bce", "Bth", "Bwe")

 mymethod = "average"; #111808 
 hd =  hclust( dist(etb2), mymethod); 
 essen.cat = cutree(hd, h=0); 
 table(essen.cat);

 conserved.essen = essen.cat[essen.cat==1]
 #there are many multiple hits. 
 #I should simple best hits
 
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


quit("yes");

