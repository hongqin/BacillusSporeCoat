
 rm(list=ls());

 ### read in Henriques07 coat profile
 Hcoat = read.csv( 'Henriques07.coat.profile.122108.csv',header=T);
 Hcoat[,1] = as.character(Hcoat[,1]);
 Hcoat[,2] = as.character(Hcoat[,2]);

 old.col = names(Hcoat);

 HBG = sort( Hcoat[,2] );

 Hcoat = Hcoat[ match(HBG, Hcoat$NCBI), ];
 cols = c("B.s.","B.li.","B.a.","B.c.","B.t.","B.cl.","B.h.");
 Hcoat3 = Hcoat[, match(cols, old.col) ]
 row.names(Hcoat3) = HBG;
 summary(Hcoat3); 

 ### read the profile for coat genes 
 ctb = read.table( "_coat.profile.121708.csv",sep="\t",header=T);
 ctb$id = as.character( ctb$id );
 row.names( ctb ) = as.character( ctb$id );
 bacillus.specs = names(ctb)[10:20];

 ctb2 = ctb[, c(1,10,13,15:17,19:20)];
 names(ctb2) =  c( "NCBI", names(ctb2)[2:8] );

 ctb3 = ctb2[, 2:8] 

 #row.names(ctb4) = as.character( ctb4$Bsu );
 ctb3 = ifelse( is.na(ctb3), 0, 1);
 ctb3 = data.frame(ctb3);
 head(ctb3); #passed

###############
 ### compare ctb3 and Hcoat3 
 table( HBG %in% ctb$id )
#FALSE  TRUE 
#   10    60

 HC  = Hcoat3[HBG %in% ctb$id, ] 
 MyC = ctb3[ctb$id %in% HBG, ] 
 names(HC) = names(MyC);

 #reorder them
 HC  = HC[match(sort(row.names(HC)),row.names(HC)), ]
 MyC = MyC[match(sort(row.names(MyC)),row.names(MyC)), ]

 #check the orders
 hrow = row.names( HC );
 myrow = row.names( MyC );
 plot( match( hrow, myrow), 1:60) ; #OK, a straight line. 

 ### calculate the Hamming distance
 d = numeric(length(HC[,1]));
 names(d) = row.names( HC );

 for( row in 1:length(HC[,1]) ){
    #d[row] = hamming.distance( rbind(HC[row,], MyC[row,]) );
    d[row] = dist( rbind(HC[row,], MyC[row,]) ); #Euclidean distance
 }
 d = d^2; #convert to hamming distance
 
 table(d);
 hist(d);

 #store the updated hamming distances
 ctb$hamming2 = d[ match(ctb$id, names(d) ) ];
 ctb$hamming2 = as.integer(floor(ctb$hamming2 + 0.5) );

 #consistency check 
 x = ctb[,c("hamming","hamming2")]
 x$diff = x[,1] - x[,2];
 x[x$diff>0,]

 ### which genes have different profiles?
 diff = ctb[ ctb$hamming2>0 & (! is.na(ctb$hamming2) ), ]

 ### which genes in Henriques07 are not in ours?
 Hcoat$hamming = d[ match(Hcoat$NCBI, names(d) ) ]
 Hcoat$HenriUniq = ifelse( is.na(Hcoat$hamming), 1, 0); 
 Hcoat.uniq = Hcoat[Hcoat$HenriUniq==1,]

 ### output the new files
 write.csv(Hcoat.uniq, "uniq.Henriques07.coat.122208.csv",row.names=F)
 write.csv(ctb, "Bacillus_subtilis_coat_proteins.122208.csv", row.names=F)

quit('yes');














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
