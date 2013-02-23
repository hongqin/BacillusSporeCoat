
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

 # now conver hit id to 1 and 0 in ctb
 ctb3 = ctb[,-c(10:20)]; 
 ctb3 = cbind( ctb3, profile );

 coatall2 = ctb3; #just for convenience

 table( coatall2$CoatLocation );

 coatall2$coatflag = 1;

### read the coordinates
 coor = read.csv("Bsu168.features.coor.csv", sep="\t",header=F,as.is=T,  colClasses=c("character","character","character","numeric","numeric") );

colnames(coor) = c("bg","type","strand","start","end");
rownames(coor) = coor$bg;

 coor$pos = 1:length(coor[,1])
 coor$pos = ifelse(coor$strand=='+', coor$pos, -coor$pos);

 #sort by $start, make sure it is correct
 sort.start = sort( coor$start );
 coor = coor[match(sort.start, coor$start), ]

 #visual check
 x = 1:length(coor[,1])
 plot( coor$pos ~ x );

### find consecutive tandem clusters 
 pos = coor$pos
 names(pos) = coor$bg
 
 h =  hclust(dist(pos), method="single"); 
 clusters = cutree(h, h=1);
 
 table(clusters); #check, looks good, 1100 consecutive tandem clusters

 coor$clusters = clusters[coor$bg]

 #map back to coatall2
 coatall2$clusters = clusters[coatall2$id]

###visual check of tandem clusters
 coatall2[,c("id","CoatLocation","clusters")]
 myids = paste( "BG10",490:500,sep='');
 coatall2[myids,]
 coor[myids,]

##noncoat clusters
 coor$coatflag = coatall2$coatflag[match(coor$bg, coatall2$id)]
 noncoat = coor[ is.na(coor$coatflag), ]
 
## coat clusters
 coatall2$strand   = coor$strand[ match(coatall2$id, coor$bg)]
 coatall2$start    = coor$start[ match(coatall2$id, coor$bg)]
 coatall2$end      = coor$end[ match(coatall2$id, coor$bg)]
 
 x = sort( coatall2$start );
 
 # coatall3 = coatall2[match(x, coatall2$start), ]
 coatall3 = coatall2[match(x, coatall2$start), c("id", "Gene", "CoatLocation","clusters", "strand")]

 coatall3

## do coat proteins tend to cluster in the genome? 
### test homogeneity of coat location among tandem clusters by permutation
 coor$coatflag = ifelse( is.na(coor$coatflag), 0, 1)
 coor.original = coor;
 coatclusters = coor$clusters[coor$coatflag==1]
 uniq.coat.clusters = length(unique(coatclusters))

 N=10000;
 #random.coat.cluster.numbers = uniq.coat.clusters;
 random.coat.cluster.numbers = NA;

 for( i in 1:N) {
   coor$coatflag2 = sample(coor$coatflag) ;
   coatclusters.tmp = coor$clusters[coor$coatflag2==1]
   uniq.coat.clusters.tmp = length(unique(coatclusters.tmp))
   random.coat.cluster.numbers = c(random.coat.cluster.numbers,uniq.coat.clusters.tmp);
 }

 mean(random.coat.cluster.numbers, na.rm=T);
 sd(random.coat.cluster.numbers, na.rm=T);
 x = random.coat.cluster.numbers;
 y = x[x<=uniq.coat.clusters  ]
 y = y[! is.na(y)]
 p = length(y) / N
 #So, the coat proteins are indeed highly clustered together, p=0.0015


 MyBox <- function( x0,y0,h,w,col) {
   xin = c( x0, x0,   x0+w, x0+w );
   yin = c( y0, y0+h, y0+h, y0   );
   polygon(xin,yin, col=col, border=NA);
 }

his = hist(random.coat.cluster.numbers);

pdf("_123008.coat.cluster.pdf",width=7,height=7);
plot(his,xlim=c(50,80), main="");
arrows(57,500,57,50,lwd=1,col="red", length=0.1) 
#MyBox( 57,0,500, 0.1, "red") #ugly
text(57,600, "coat prot clusters")
dev.off();



######################################################
## given the current coat protein clusters, do inner and outer coat protein
## have preference of other members in clusters? 
## permutation to generate a Z-score matrix

 myOutCats = c("InnerCoat","OuterCoat","Other","Unkown");

 myInCats  = unique(coatall3$CoatLocation)

 #re-partition the CoatLocation to myOutCats
 coatall3$CoatLocation = as.character( coatall3$CoatLocation ) 
 coatall4 = coatall3;
 for( row in 1:length(coatall4[,1]) ) {
   if ( ( coatall4$CoatLocation[row] !="InnerCoat") & ( coatall4$CoatLocation[row] !="OuterCoat") ) {
      if ( coatall4$CoatLocation[row] == myInCats[1] ) {
         coatall4$CoatLocation[row] = "Unkown"
      } else {
         coatall4$CoatLocation[row] = "Other"
      }    
   }
 }

 myclusters = unique( coatall4$clusters ); 

 observation = matrix(0,nrow=length(myOutCats), ncol=length(myOutCats));
 colnames(observation) = myOutCats;
 rownames(observation) = myOutCats;

 binomialcoef <- function(n) { n*(n-1)/2; }

 #count the observations
 for( mycluster in myclusters ) {
   sub = coatall4[coatall4$clusters==mycluster, ]
   counts = table(sub$CoatLocation);
   for( i in 1:length(counts) ) {
      cat = names(counts)[i]
      observation[cat,cat] = observation[cat,cat]  + binomialcoef(counts[cat])
   }
   
   if (length(counts)>=2) { 
   len = length(counts); 
   for( i in 1:(len-1) ) {
      for( j in i:len ) {
       cati = names(counts)[i]
       catj = names(counts)[j]
       observation[cati,catj] = observation[cati,catj]+ counts[cati]*counts[catj] 
       observation[catj,cati] = observation[catj,cati]+ counts[catj]*counts[cati]       
      }
   }
   }
 }
observation

 #bookeeping 
 coatall4.original = coatall4;
 observation.original = observation; 

 ##do permutations
 cat1s = NA;
 cat2s = NA;
 edges = NA;
 iterations = NA;

 N=1000;
rep=1;
while( rep<=N ) {
 coatall4$CoatLocation = sample( coatall4$CoatLocation); #permutation

 observation = matrix(0,nrow=length(myOutCats), ncol=length(myOutCats));
 colnames(observation) = myOutCats;
 rownames(observation) = myOutCats;

 for( mycluster in myclusters ) {
   sub = coatall4[coatall4$clusters==mycluster, ]
   counts = table(sub$CoatLocation);
   for( i in 1:length(counts) ) {
      cat = names(counts)[i]
      observation[cat,cat] = observation[cat,cat]  + binomialcoef(counts[cat])
   }
   
   if (length(counts)>=2) { 
   len = length(counts); 
   for( i in 1:(len-1) ) {
      for( j in i:len ) {
       cati = names(counts)[i]
       catj = names(counts)[j]
       observation[cati,catj] = observation[cati,catj]+ counts[cati]*counts[catj] 
       observation[catj,cati] = observation[catj,cati]+ counts[catj]*counts[cati]       
      }
   }
   }
 }

 #bookkeeping
   for( i in 1:3 ) {
      for( j in (i+1):4 ) {
       cat1s = c(cat1s, myOutCats[i]);
       cat2s = c(cat2s, myOutCats[j]);
       edges= c(edges,observation[i,j]);
       iterations = c(iterations, rep);
      }
   }
   for( i in 1:4 ) {
       cat1s = c(cat1s, myOutCats[i]);
       cat2s = c(cat2s, myOutCats[i]);
       edges= c(edges,observation[i,i]);
       iterations = c(iterations, rep);
   }
rep = rep + 1;
}

permOut = data.frame(cbind(cat1s,cat2s));
permOut$edges = edges;
permOut$iterations = iterations;
permOut = permOut[-1,]
permOut$cat1 = as.character(permOut$cat1s);
permOut$cat2 = as.character(permOut$cat2s);

 meanM= matrix(0,nrow=length(myOutCats), ncol=length(myOutCats));
 colnames(meanM) = myOutCats;
 rownames(meanM) = myOutCats;
 sdM = meanM;
 pM = meanM;
 ZM = meanM;

 for( i in 1:4 ) {
      for( j in i:4 ) {
       cat1 = myOutCats[i];
       cat2 = myOutCats[j];
       sub = permOut[permOut$cat1==cat1 & permOut$cat2==cat2, ]
       meanM[i,j] = mean( sub$edges );
       meanM[j,i] = meanM[i,j]
       sdM[i,j] = sd( sub$edges ); sdM[j,i] = sdM[i,j];
       pM[i,j]  = length(sub$edges[sub$edges>= observation.original[i,j]]) / length(sub[,1]);
       pM[j,i] = pM[i,j];
       if(sdM[i,j]>0) {
        ZM[i,j]  = (observation.original[i,j] - meanM[i,j] ) / sdM[i,j]; 
        ZM[j,i]  = ZM[i,j]
       } else {
        ZM[i,j] = NA; ZM[j,i] = NA;
       } 
      }
 }

library(Heatplus)
colnames(pM) = rep(NA,4);

 pdf("_gene.neighborhood.123008.pdf",width=8,height=5);
heatmap_2(pM,scale='none',do.dendro=c(F,F), Rowv=NA, Colv=NA, legend=2, legfrac=6, 
 col = RGBColVec(64)
)
dev.off();


q("yes");

###################old codes

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

