
 rm(list=ls());

 ### read the 110808 profile for coat genes
 ctb = read.table( "_coat.profile.110808.csv",sep=",",header=T);
 ctb$id = as.character( ctb$id );
 row.names( ctb ) = as.character( ctb$id );

 names(ctb)[ -c(2:9,21:38) ]
 ctb2 = ctb[,-c(2:9,21:38)];
 names(ctb2) =  c( "NCBI", names(ctb2)[2:12] );
 bacillus.specs = names(ctb)[10:20];

 plot( match( ctb$id, ctb$NCBI) ); #check passed 

 profile = ifelse( is.na(ctb2[,2:12]), 0, 1 );

 profile = data.frame( profile );

### now conver hit id to 1 and 0 in ctb
 ctb3 = ctb[,-c(10:20)]; 
 ctb3 = cbind( ctb3, profile );

 all2 = ctb3; #just for convenience

 table( all2$CoatLocation );

### fisher.test
subtb = table( all2[,c("CoatLocation","Bce")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2 ); #p=0.50

subtb = table( all2[,c("CoatLocation","Bth")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2,alternative="less" ); #p=0.035

subtb = table( all2[,c("CoatLocation","Bwe")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2, alternative="less" ); #p=0.0065

subtb = table( all2[,c("CoatLocation","Ban")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2, alternative="less" ); #p = 0.029

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
wilcox.test( all3$hits.cereus[all3$CoatLocation=="InnerCoat"], all3$hits.cereus[all3$CoatLocation=="OuterCoat"], alternative="greater"); #p=0.045
#This p-value changes = 0.096 if c("Bce","Ban","Bwe","Bth", "Bha","Bcl","Bpu","Bli" are used 

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
