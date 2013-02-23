
#060511 add inner and outer coat partition, but the results is not satifactory. 

# 051411 add coat profile
# I will partion Bsu, Bce-clade omegas, wilcox.test in coat proteins and between coat and nonCE protein

##########################################
rm(list=ls());
 
## some test to demonstrate the syntax for alternative hypothesis
#x = rnorm(100)+5
#y = x + 2
#summary(x); 
#summary(y)
#wilcox.test( x, y, alt='gr') #accept null, large p 
#wilcox.test( x, y, alt='less') #reject null, accept alternative, small p

## coat gene profile
 ctb = read.table( "_coat.profile.122908.csv",sep="\t",header=T);
 ctb$id = as.character( ctb$id );
 row.names( ctb ) = as.character( ctb$id );
 bacillus.specs = names(ctb)[10:20];

 ctb2 = ctb[,c(1,10:20)]
 ctb3 = ctb2[, 2:12] #profile matrix
 ctb3 = ifelse( is.na(ctb3), 0, 1);
 ctb3 = data.frame(ctb3);
 #head(ctb3); #passed

 ctb3$cereushits = ctb3$Ban + ctb3$Bth + ctb3$Bce + ctb3$Bwe
 ctb3$subtilishits = ctb3$Bam + ctb3$Bli + ctb3$Bpu 
 ctb3$outhits = ctb3$Bha + ctb3$Bcl
 ctb3$hits = ctb3$cereushits + ctb3$subtilishits + ctb3$outhits;  
 #hist(ctb3$cereushits)

### coat free omega, and match to coat orthologous hits in the profile table
coat.tb = read.delim("_free.omega.coat.paml.Jan20,2011.txt", header=F);
coat.tb = coat.tb[ ! is.na(coat.tb$V4), ]
coat.tb$V2 = as.character(coat.tb$V2) 

positions = match( coat.tb$V1, rownames(ctb3))
#ctb3[positions[1:5], ] #check, passed
coat.tb = cbind( coat.tb, ctb3[positions, c('cereushits','subtilishits', 'outhits', 'hits')])
#head(coat.tb)

positions = match( coat.tb$V1, ctb$id)
coat.tb = cbind( coat.tb, ctb$CoatLocation[positions])
tmp = names(coat.tb)
names(coat.tb)= c(tmp[-23], 'CoatLocation')

### nonCE free omega
nonCE.tb = read.delim("_nonCE.omega.7March2011.tab", header=F);
nonCE.tb = nonCE.tb[nonCE.tb$V2=='H1C',] # H1C is where free omega are calculated 
nonCE.tb$V3 = as.character(nonCE.tb$V3)

#remove 999s Many labile coat genes have omega=999?? 
#Wcutoff = 998
#coat.tb$V10[coat.tb$V10> Wcutoff ] = NA
#nonCE.tb$V11[nonCE.tb$V11> Wcutoff ] = NA


#### parse the omega results 
##need to add more
clades = unique(c( nonCE.tb$V3, coat.tb$V2))
# "(Bha,Bcl)"             "Bha"                   "Bcl"                   
#"(Bcl,Bli,Bam,Bsu)"    

#BsuClade = c("(Bpu,Bli,Bam,Bsu)","(Bli,Bam,Bsu)","(Bpu,Bli,Bsu)", "Bpu", "Bli", "(Bam,Bsu)", "Bam", "Bsu", "(Bpu,Bam,Bsu)", "(Bpu,Bli,Bam,Bmo,Bsu)", "(Bli,Bam,Bmo,Bsu)", "(Bam,Bmo,Bsu)","(Bmo,Bsu)","Bmo","(Bli,Bam,Bsu,Bmo)", "(Bam,Bsu,Bmo)", "(Bsu,Bmo)", "(Bpu,Bam,Bmo,Bsu)", "(Bli,Bmo,Bsu)", "(Bli,Bsu)" );   #p=0.0305   

#leaf branches for publication     
BsuClade = c("Bpu", "Bli", "Bam", "Bsu");   #p-value = 0.00823,

# internal branches for publication
#BsuClade = c("(Bpu,Bli,Bam,Bsu)","(Bli,Bam,Bsu)","(Bpu,Bli,Bsu)", "Bpu", "Bli", "(Bam,Bsu)", "Bam", "(Bpu,Bam,Bsu)", "(Bpu,Bli,Bam,Bmo,Bsu)", "(Bli,Bam,Bmo,Bsu)", "(Bam,Bmo,Bsu)","(Bmo,Bsu)","Bmo","(Bli,Bam,Bsu,Bmo)", "(Bam,Bsu,Bmo)", "(Bsu,Bmo)", "(Bpu,Bam,Bmo,Bsu)", "(Bli,Bmo,Bsu)", "(Bli,Bsu)" );   #no Bsu leaf, p=0.13

#BsuClade = c("(Bli,Bam,Bsu)","Bli","(Bam,Bsu)","Bam","Bsu","(Bpu,Bam,Bsu)", "(Bli,Bam,Bmo,Bsu)", "(Bam,Bmo,Bsu)", "(Bmo,Bsu)","Bmo","(Bli,Bam,Bsu,Bmo)", "(Bam,Bsu,Bmo)","(Bsu,Bmo)","(Bli,Bmo,Bsu)","(Bli,Bsu)" );  # p-value = 0.01123, no Bpu

#BsuClade = c("(Bpu,Bli,Bam,Bsu)","(Bli,Bam,Bsu)","(Bpu,Bli,Bsu)", "(Bam,Bsu)", "(Bpu,Bam,Bsu)", "(Bpu,Bli,Bam,Bmo,Bsu)", "(Bli,Bam,Bmo,Bsu)", "(Bam,Bmo,Bsu)","(Bmo,Bsu)","(Bli,Bam,Bsu,Bmo)", "(Bam,Bsu,Bmo)", "(Bsu,Bmo)", "(Bpu,Bam,Bmo,Bsu)", "(Bli,Bmo,Bsu)", "(Bli,Bsu)" );   #no leaf nodes, high omega all over the places, large p-value
#so, these suggest coat gene only recently evolves faster than the nonCE genes, suggest adpation to niches



#BsuClade = c("Bli", "Bam", "Bsu");   #0.01082
#BsuClade = c("Bpu", "Bam", "Bsu");   # p-value = 0.002779
BsuClade = c("Bsu") # p-value = 0.002939
#**BsuClade = c("Bam") # p=0.0858
#BsuClade = c("Bsu", "Bmo", "(Bsu,Bmo)","(Bmo,Bsu)") # p-value = 3.483e-08
#BsuClade = c("Bsu", "Bmo") #  p-value = 4.654e-06
#BsuClade = c("Bsu", "(Bsu,Bmo)","(Bmo,Bsu)") # p-value = 1.899e-05

#BsuClade = c("Bli") # p=0.5 ??
#BsuClade = c("Bpu") # p=0.19
#BsuClade = c("(Bam,Bmo,Bsu)", "(Bsu,Bmo)", "(Bam,Bsu)") #large p
#BsuClade = c( "(Bam,Bsu)") #p=0.1568
     
BceClade = c("(Bwe,Bce,Ban,Bth)","(Bce,Ban,Bth)","(Ban,Bth)","(Bwe,Bce,Bth)","(Bce,Bth)", "(Bwe,Ban,Bth)","Bwe","Ban","Bth","Bce"); # p=0.1694 cCW>nCW
#BceClade = c("(Bwe,Bce,Ban,Bth)","(Bwe,Bce,Bth)", "(Bwe,Ban,Bth)"); # ?? ambiguous

#BceClade = c("(Bwe,Bce,Ban,Bth)","(Bce,Ban,Bth)","(Ban,Bth)","(Bwe,Bce,Bth)","(Bce,Bth)", "(Bwe,Ban,Bth)"); #p=0.214 cCW > nCW
#BceClade = c("Ban","Bth","Bce"); #p=0.25 cCW < nCW, 
#BceClade = c("(Bce,Ban,Bth)","(Ban,Bth)","(Bce,Bth)"); #p=0.1886
#BceClade = c("(Bce,Ban,Bth)","(Ban,Bth)","(Bce,Bth)", "Ban","Bth","Bce"); # p = 0.157

nonCE.tb$flag = NA;
coat.tb$flag  = NA; 
for( j in 1:length(BsuClade)){
  nonCE.tb$flag[ nonCE.tb$V3 == BsuClade[j] ] = 'bsu';
  coat.tb$flag[ coat.tb$V2==BsuClade[j] ]     = 'bsu';
}

for( j in 1:length(BceClade)){
  nonCE.tb$flag[ nonCE.tb$V3 == BceClade[j] ] = 'bce';
  coat.tb$flag[ coat.tb$V2==BceClade[j] ]     = 'bce';
}

table( coat.tb$V2, coat.tb$flag) #check flags, passed. 051411

summary(lm( coat.tb$V10 ~ coat.tb$CoatLocation )) #not significant

nBW =  nonCE.tb$V11[nonCE.tb$flag=='bsu'];
nCW =  nonCE.tb$V11[nonCE.tb$flag=='bce'];
cBW =  coat.tb$V10[coat.tb$flag=='bsu']; 
cCW =  coat.tb$V10[coat.tb$flag=='bce']; 

cBW.outer =  coat.tb$V10[coat.tb$flag=='bsu' & coat.tb$CoatLocation == 'OuterCoat']; 
cCW.outer =  coat.tb$V10[coat.tb$flag=='bce' & coat.tb$CoatLocation == 'OuterCoat']; 
cBW.inner =  coat.tb$V10[coat.tb$flag=='bsu' & coat.tb$CoatLocation == 'InnerCoat' ]; 
cCW.inner =  coat.tb$V10[coat.tb$flag=='bce' & coat.tb$CoatLocation == 'InnerCoat' ]; 

####################### summary


summary(nBW)
summary(cBW)
summary(cBW.outer)
summary(cBW.inner)

summary(nCW)
summary(cCW)
summary(cCW.outer)
summary(cCW.inner)
