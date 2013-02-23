library(Heatplus)

myfiles = c("pmat.dN-coat-nonCE-wilcox-greater.tab",
#"pmat.dN-essen-nonCE-wilcox-greater.tab",
"pmat.dS-coat-nonCE-wilcox-greater.tab",
#"pmat.dS-essen-nonCE-wilcox-greater.tab", #All '1's, eror in heatmap
"pmat.omega-coat-nonCE-wilcox-greater.tab",
"pmat.omega-essen-nonCE-wilcox-greater.tab");

strains = c("Ban","Bth","Bce","Bce3L", "Bce87","Bwe","Bsu","Bam","Bli","Bpu","Bcl","Bha");

pdf ("sandbox/all.071310.pdf", width=5, height=5)
for( fl in myfiles ) {
 #ptb = read.table("pmat-coat-nonCE-wilcox-greater.tab")
 ptb = read.table(fl)
 colnames(ptb) = strains;
 rownames(ptb) = strains; 
 pmat = as.matrix(ptb)
 heatmap_2( pmat, do.dendro = c(FALSE, FALSE), legend = 4, legfrac = 6, scale="none",
  #hclustfun = function(x) hclust( as.dist(x), method= "average"),
  Rowv = NA, Colv=NA,  col = RGBColVec(64),  
  main = fl ); 

 ptb = ptb[ -c(4,5), -c(4,5)]
 pmat = as.matrix(ptb)
 heatmap_2( pmat, do.dendro = c(FALSE, FALSE), legend = 4, legfrac = 6, scale="none",
  #hclustfun = function(x) hclust( as.dist(x), method= "average"),
  Rowv = NA, Colv=NA,  col = RGBColVec(64),  
  main = fl ); 

}
dev.off()

quit("no")

 strains = c("Ban","Bce","Bth","Bwe","Bsu","Bam","Bli","Bpu","Bcl","Bha");
 coat = read.table("_omega.coat.summary.121508.csv",row.name=1, header=T);
 cbind(strains, names(coat) ); #visual check, passed

 colnames(coat) = strains;
 rownames(coat) = strains;

 dc = as.dist( coat );
 plot( hclust(dc, method="ave"), main="omega coat" );

 ###do heat map
 library(Heatplus);
 dat = as.matrix(coat);
 pdf("omega.coat.average.121508.pdf",width=5,height=5);
 heatmap_2( dat, do.dendro = c(TRUE, FALSE), legend = 1, legfrac = 6, scale="none",
  hclustfun = function(c) hclust( dc, method= "average"),
  #col = RGBColVec(64),  main = "coat omega average-linkage" ); 
  col = rev(gray.colors(20))) #  main = "coat omega average-linkage" ); 

 heatmap_2( dat, do.dendro = c(FALSE, FALSE), legend = 1, legfrac = 5, scale="none",
  hclustfun = function(c) hclust( dc, method= "average"),
  #col = RGBColVec(64),  main = "coat omega average-linkage" ); 
  col = rev(gray.colors(20)) )#  main = "coat omega average-linkage" ); 

  dc = as.dist( coat );
  plot( hclust(dc, method="ave") ,  main="" );

 dev.off();

 heatmap_2( dat, do.dendro = c(TRUE, FALSE), legend = 1, legfrac = 8, scale="none",
  hclustfun = function(c) hclust( dc, method= "complete"),
  col = RGBColVec(64),  main = "coat omega complete-linkage" ); 

 heatmap_2( dat, do.dendro = c(TRUE, FALSE), legend = 1, legfrac = 8, scale="none",
  hclustfun = function(c) hclust( dc, method= "ward"),
  col = RGBColVec(64),  main = "coat omega ward-linkage" ); 

 #dat = as.matrix(coat);
 #dat[upper.tri(dat)] <- NA;
 #image(dat, col=RGBColVec(64));
 
 #library(LDheatmap);
 #LDheatmap(dat,dat); #not very useful here

######### essential genes omega

 essen = read.table("_omega.essential.summary.121508.csv",row.name=1,header=T);
 cbind(strains, names(essen) ); #visual check, passed
 colnames(essen) = strains;
 rownames(essen) = strains;

 de = as.dist( essen );
 plot( hclust(de, method="ave"), main="omega essential genes") ;

 dat = as.matrix(essen);
 pdf("omega.essential.average.121508.pdf",width=5,height=5);
 heatmap_2( dat, do.dendro = c(TRUE, FALSE), legend = 1, legfrac = 6, scale="none",
  hclustfun = function(c) hclust( de, method= "average"),
  #col = RGBColVec(64),  main = "essential genes omega average-linkage" ); 
  col = rev(gray.colors(20))) #  main = "essential genes omega average-linkage" ); 

 heatmap_2( dat, do.dendro = c(FALSE, FALSE), legend = 1, legfrac = 5, scale="none",
  hclustfun = function(c) hclust( de, method= "average"),
  #col = RGBColVec(64),  main = "essential genes omega average-linkage" ); 
  col = rev(gray.colors(20)) )#  main = "essential genes omega average-linkage" ); 

  de = as.dist( essen );
  plot( hclust(de, method="ave") ,  main="" );

 dev.off();


##### noncoat omega

 noncoat = read.table("_omega.noncoat.summary.121508.csv",row.name=1,header=T);
 cbind(strains, names(noncoat) ); #visual check, passed
 colnames(noncoat) = strains;
 rownames(noncoat) = strains;

 dnc = as.dist( noncoat );
 plot( hclust(dnc, method="ave") );

 plot( hclust(dc/de, method="ave"), main="omega coat/ omega essential)" );

 dat = as.matrix(noncoat);
 pdf("omega.noncoat.average.121508.pdf",width=5,height=5);
 heatmap_2( dat, do.dendro = c(TRUE, FALSE), legend = 1, legfrac = 6, scale="none",
  hclustfun = function(c) hclust( dnc, method= "average"),
  #col = RGBColVec(64),  main = "noncoat omega average-linkage" ); 
  col = rev(gray.colors(20))) #  main = "noncoat omega average-linkage" ); 

 heatmap_2( dat, do.dendro = c(FALSE, FALSE), legend = 1, legfrac = 5, scale="none",
  hclustfun = function(c) hclust( dnc, method= "average"),
  #col = RGBColVec(64),  main = "noncoat omega average-linkage" ); 
  col = rev(gray.colors(20)) )#  main = "noncoat omega average-linkage" ); 

  dnc = as.dist( noncoat );
  plot( hclust(dnc, method="ave") ,  main="" );

 dev.off();

##### pomega
 strains = c("Ban","Bce","Bth","Bwe","Bsu","Bam","Bli","Bpu","Bcl","Bha");
 pomega = read.table("_pomega.coat-noncoat.summary.121508.csv",row.name=1, header=T);
 cbind(strains, names(pomega) ); #visual check, passed

 #for( row in 1:length(pomega[1,]) ) {
 #  for( col in 1:length(pomega[,1]) ) {
 #     pomega[row,col] = ifelse( pomega[row,col]==0, 1E-9, pomega[row,col] );
 #  }
 #}
 
 dpomega = -log10(pomega); 
 dp = as.dist( dpomega );
 plot( hclust(dp, method="ave") );
 
 dp2 = dist( pomega );
 plot( hclust(dp2, method="ave") ); #This seems to be best tree

####
  plot( hclust(dc, method="ave") ,  main=" coat omega" );
  plot( hclust(de, method="ave") ,  main=" essential omega" );
  plot( hclust(dc/de, method="ave") ,  main=" coat omega/essen omega" );

  plot( hclust(dnc, method="ave") ,  main=" noncoat omega" );
  plot( hclust(dc/dnc, method="ave") ,  main=" coat omega /noncoat omega" );

