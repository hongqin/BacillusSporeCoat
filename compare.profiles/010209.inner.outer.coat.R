
 rm(list=ls());

 ### read the coat profile 
 ctb = read.csv( "_coat.profile.122908.csv",sep="\t");
 ctb$id = as.character( ctb$id );
 row.names( ctb ) = as.character( ctb$id );

 names(ctb)[ -c(2:9,21:38) ]
 ctb2 = ctb[,-c(2:9,21:38)];
 names(ctb2) =  c( "NCBI", names(ctb2)[2:12] );
 bacillus.specs = names(ctb)[10:20];

 plot( match( ctb$id, ctb$NCBI) ); #check passed 

 profile = ifelse( is.na(ctb2[,2:12]), 0, 1 );
 profile = data.frame( profile );

 ctb3 = ctb[,-c(10:20)]; 
 profile$id = rownames(profile);

 all2 = merge( ctb3, profile); #

 table( all2$CoatLocation );

### remove the cereus conserved ones
all2$hits = all2$Bce + all2$Ban + all2$Bth + all2$Bwe

labile = all2[all2$hits<4, ] #31 genes
cereus4 = all2[all2$hits==4, ]#42 genes

table( labile$CoatLocation );
table( cereus4$CoatLocation );

### Table 2
 #c4  = c(17,8,5,12);
 #a73 = c(23,20,5,25);
 c4  = c(17,8);
 a73 = c(23,20);
 d = a73 - c4;
sub2 = rbind( d, c4 );
fisher.test(sub2); #Table 2A
fisher.test(sub2, alternative="less"); #Table 2A


### fisher.test, Table 2B
subtb = table( all2[,c("CoatLocation","Bce")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2, alternative="less" ); #p=0.5206

subtb = table( all2[,c("CoatLocation","Bth")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2,alternative="less" ); #p=0.0964

subtb = table( all2[,c("CoatLocation","Bwe")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2, alternative="less" ); #p=0.01187

subtb = table( all2[,c("CoatLocation","Ban")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2, alternative="less" ); #p = 0.0972

subtb = table( all2[,c("CoatLocation","Bha")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2 ); #p = 0.33

subtb = table( all2[,c("CoatLocation","Bam")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2 ); #p = 1

subtb = table( all2[,c("CoatLocation","Bcl")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2 ); #p = 0.2

#this produces erros
#subtb = table( all2[,c("CoatLocation","Bwe", "Bha")] ); #
#sub2 = subtb[c("InnerCoat","OuterCoat"),]
#fisher.test( sub2 ); #p = 0.52


#now try hits
for( row in 1:length(all2[,1]) ) {
 all2$hits.cereus[row] = sum(all2[ row, c("Bce","Ban","Bwe","Bth")]);
# all2$hits.cereus[row] = sum(all2[ row, c("Bce","Ban","Bwe","Bth", "Bha","Bcl","Bpu","Bli")]);
}

all3 = all2[ all2$CoatLocation=="InnerCoat" | all2$CoatLocation=="OuterCoat",  ]

summary( lm( all3$hits.cereus ~ all3$CoatLocation) ); #p=0.08 for sum(all2[ row, c("Bce","Ban","Bwe","Bth", "Bha")]);

summary( all3$hits.cereus[all3$CoatLocation=="InnerCoat"]);
summary( all3$hits.cereus[all3$CoatLocation=="OuterCoat"]);

#Table 2B
wilcox.test( all3$hits.cereus[all3$CoatLocation=="InnerCoat"], all3$hits.cereus[all3$CoatLocation=="OuterCoat"], alternative="greater"); #p=0.045
#This p-value changes = 0.096 if c("Bce","Ban","Bwe","Bth", "Bha","Bcl","Bpu","Bli" are used 
t.test( all3$hits.cereus[all3$CoatLocation=="InnerCoat"], all3$hits.cereus[all3$CoatLocation=="OuterCoat"], alternative="greater"); #p=0.045

##Change to histogram comparison, 2013 Jan 22
hIn =  hist( all3$hits.cereus[all3$CoatLocation=="InnerCoat"], br=c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5)) 
hOut = hist( all3$hits.cereus[all3$CoatLocation=="OuterCoat"], br=c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5)  ) 

source("barplot3D.Rfunction.R")
# lifespan/barplot-example-1.r: barplot( as.matrix(tb), beside=T, col=c("red","blue"), ylab="Relative Density", xlab="Repl
barInOut = data.frame( rbind(hIn$counts, hOut$counts))
names(barInOut) = c(0,1,2,3,4) 
barplot( as.matrix(barInOut), beside=T)
 
#pdf( "barplotInOut2013Jan22.pdf" );
barplot( as.matrix(barInOut), beside=T, col=c("black","gray"), ylab="Counts", xlab="Number of Bce species with orthologus hits"
         , legend = c( "Inner Coat", "Outer Coat" ), xlim = c(15,0)
         );
#legend(0.5,15, c( "Inner Coat", "Outer Coat" ))
# title(main="Comparison of orthologous distributions for inner and outer coat genes" )
 #dev.off();  
 
### gene length, xnu length

 ## all coat protein
 boxplot( ctb$len ~ ctb$CoatLocation );

 wilcox.test( ctb$len[ctb$CoatLocation=="OuterCoat"], ctb$len[ctb$CoatLocation=="InnerCoat"]) #p=0.67

 boxplot( ctb$xnuLen ~ ctb$CoatLocation );

 ## only xnu proteins
  sub = all2[all2$xnuLen>0,  ]
  boxplot( sub$xnuLen ~ sub$CoatLocation );
  boxplot( sub$xnuFraction ~ sub$CoatLocation );
  boxplot( sub$len ~ sub$CoatLocation );

  t.test( sub$len[sub$CoatLocation=="OuterCoat"], sub$len[sub$CoatLocation=="InnerCoat"], alt="less") #p=0.05797

  t.test( sub$xnuLen[sub$CoatLocation=="OuterCoat"], sub$xnuLen[sub$CoatLocation=="InnerCoat"], alt="less") 



q("yes");
