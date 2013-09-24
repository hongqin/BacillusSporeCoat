
 rm(list=ls());

 library("e1071")

 ### read the profile for coat genes
 ctb = read.table( "_summary.coat.clusters.I3.1.073108.csv",sep="\t",header=T);
 ctb2 = ctb[, -c(2:9, 21)];
 names(ctb2) =  c( "NCBI", names(ctb2)[2:12] );
 bacillus.specs = names(ctb)[10:20];

 ##=> rename the duplicates
 tmp = ctb2[c(14,33,61),]
 ctb3 = rbind (ctb2, tmp );
 ctb3$NCBI = as.character( ctb3$NCBI );
 ctb3$NCBI[c(14,33,61, 71:73)] = c("BG10499","BG11791","BG13471","BG10498","BG13097","BG10492");
 row.names( ctb3) = ctb3$NCBI;

### read in the coat annotation table
 tb = read.table("BGcoat3.csv", header=T, sep=",")
 row.names(tb) = as.character( tb$NCBI );

 ids = sort( ctb3$NCBI );
 tb2 = tb[ match( ids, tb$NCBI),  ]
 tb[tb$Gene=="CotS",]
 tb2[tb2$Gene=="CotS",] #check, OK
 tb = tb2; 

### merge tb and ctb3
 all = data.frame( cbind(tb,matrix(nrow=73,ncol=length(ctb3[1,]))) )
 names(all) = c( names(tb), names(ctb3) );

# all[ , bacillus.specs] = ctb3[ match(tb$NCBI,ctb3$NCBI), bacillus.specs]
 all[ , 11:22] = ctb3[ match(tb$NCBI,ctb3$NCBI), ]

head(all$Bmo);
head(ctb3$Bmo)
head(all);
all[, c(2, 11)] #check, OK

### regression
  all2 = all;
  all2[,12:22] = ifelse ( is.na(all2[,12:22]), 0, 1 )

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
fisher.test( sub2 ); #p=0.50

subtb = table( all2[,c("CoatLocation","Bth")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2,alternative="less" ); #p=0.13

subtb = table( all2[,c("CoatLocation","Bwe")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2, alternative="less" ); #p=0.072

subtb = table( all2[,c("CoatLocation","Ban")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2, alternative="less" ); #p = 0.072

subtb = table( all2[,c("CoatLocation","Bha")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2 ); #p = 0.19

subtb = table( all2[,c("CoatLocation","Bam")] ); #
sub2 = subtb[c("InnerCoat","OuterCoat"),]
fisher.test( sub2 ); #p = 0.19

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
t.test( all3$hits.cereus[all3$CoatLocation=="InnerCoat"], all3$hits.cereus[all3$CoatLocation=="OuterCoat"], alternative="greater"); #p=0.056
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
