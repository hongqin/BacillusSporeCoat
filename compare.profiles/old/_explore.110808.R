
 rm(list=ls());

 #library("e1071")

 ### read the profile for coat genes
 ctb = read.table( "_coat.profile.110608.csv",sep=",",header=T);
 ctb$id = as.character( ctb$id );
 row.names( ctb ) = as.character( ctb$id );

 ctb2 = ctb[, -c(2:9, 21)];
 names(ctb2) =  c( "NCBI", names(ctb2)[2:12] );
 bacillus.specs = names(ctb)[10:20];

### read in the coat annotation table
 ### read the annotation table
 tb = read.table( "/home/hqin/coat.protein07/key.data/Bacillus_subtilis_coat_proteins.092008.csv",
header=T,sep="\t");

 str(tb);
 tb$Gene = as.character( tb$Gene );
 tb$NCBI = as.character( tb$NCBI );

 ids = sort( tb$NCBI );
 tb2 = tb[ match( ids, tb$NCBI),  ]
 tb[tb$Gene=="CotS",]
 tb2[tb2$Gene=="CotS",] #check, OK
 tb = tb2; 

### consistence check
table( tb$NCBI %in% ctb2$NCBI ); #73 True, passed. 
cbind( tb$NCBI, ctb2$NCBI );

### merge tb and ctb2
 all = data.frame( cbind(tb[1:11],matrix(nrow=73,ncol=length(ctb2[1,]))) )
 names(all) = c( names(tb[1:11]), names(ctb2) );

 all[ ,12] = ctb2$NCBI[ match(tb$NCBI,ctb2$NCBI) ]

head(all$Bmo);
head(ctb2$Bmo)
head(all);
all[, c(2, 12)] #check, OK

# generate numeric profile
for( r in 1:73) {
  for( c in 13:23 ) {
    all[r,c] = ifelse( is.na(ctb2[r,c-11]), 0, 1); 
  }
}

### regression
  all2 = all;
  #all2[,12:22] = ifelse ( is.na(all2[,12:22]), 0, 1 )

#  summary( lm(all2$CoatLocation ~ all2$Ban+ all2$Bam+all2$Bce) );
  summary( lm( all2$Ban ~ all2$CoatLocation ) );
  summary( lm( all2$Bce ~ all2$CoatLocation ) );
  summary( lm( all2$Bth ~ all2$CoatLocation ) );
  summary( lm( all2$Bwe ~ all2$CoatLocation ) );
  summary( lm( all2$Bam ~ all2$CoatLocation ) );
  summary( lm( all2$Bpu ~ all2$CoatLocation ) );

table( all2[,c("CoatLocation","Bmo")] );

subtb = table( all2[,c("CoatLocation","Bce")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2 ); #p=0.73

subtb = table( all2[,c("CoatLocation","Bth")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2,alternative="less" ); #p=0.32

subtb = table( all2[,c("CoatLocation","Bwe")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2, alternative="less" ); #p=0.13

subtb = table( all2[,c("CoatLocation","Ban")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2, alternative="less" ); #p = 0.13

subtb = table( all2[,c("CoatLocation","Bha")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2 ); #p = 0.75

subtb = table( all2[,c("CoatLocation","Bam")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2 ); #p = 0.71

subtb = table( all2[,c("CoatLocation","Bcl")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2 ); #p = 0.52

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



q("yes");

#######################
 
 # ctb2 = cbind( rep("coat", length(ctb[,1])), ctb[, c(1,10:20) ]);
 ctb2 = data.frame( matrix( nrow=length(ctb[,1]), ncol=11) ); 
 colnames(ctb2) = bacillus.specs;
 row.names(ctb2) = as.character(ctb$X);

 for( j in 1:11){
   for( i in 1:length(ctb2[,1]) ) {
   #for( i in 1:2 ) {
	ctb2[i,j] = ifelse( is.na(ctb[i,j+9]), 0, 1 );
   }
 }
