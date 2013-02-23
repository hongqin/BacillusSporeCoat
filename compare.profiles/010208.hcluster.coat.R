# phylogeny-weighted clustering, not pretty, unused for publication. 

# re-adjust Col order and Row order
# cuttree(hd, 4)  for coat, 
# cuttree(hd, 3) for essential genes

 rm( list = ls() );
 library("e1071")

 numclus =4 ; 
 numclus =2 ; #111808 change

 etb = read.table( "_essenGene.baciProfile.091208.csv",sep="\t",header=T);	
 ctb = read.csv( "_coat.profile.122908.csv", sep="\t");

 bacillus.specs = names(ctb)[10:20];
 
 ### generate profile for coat genes
 #ctb2 = ctb[, c(10,12:20)];
 ctb2 = ctb[, c(10,12:18)];
 rownames(ctb2) = as.character( paste(ctb$id, ctb$Gene) );

 ctb3 = ctb2;

 ctb3 = ifelse( is.na(ctb3), 0, 1);
 ctb3 = data.frame(ctb3);
 head(ctb3); #passed
 
 #check
 head(ctb[,11:20]);
 head(ctb2);
# ctb2["BG11380",]

 ctb4 = ctb3;
 ctb4$Bce = ctb4$Bce * 3;
 ctb4$Bwe = ctb4$Bwe * 3;
 ctb4$Bth = ctb4$Bth * 3;
 ctb4$Ban = ctb4$Ban * 3;
 #ctb4$Bcl = ctb4$Bcl * 3;
 #ctb4$Bha = ctb4$Bha * 3;

 ctb4 = as.matrix(ctb4);
 #ctb4 = ctb4[,5:8]
 library(Heatplus);

 heatmap_2( ctb4, scale='none', legend=2, col=c("red","blue","white","green" ), hclustfun=function(x) hclust(x,method="sing") );

