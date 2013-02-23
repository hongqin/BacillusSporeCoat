#2013 Jan 13, compare free omega bw coat and nonCE after GBlocks

# Nov 29, 2011,  partition between <1, >1, 

# 051411 add coat profile
# I will partion Bsu, Bce-clade omegas, wilcox.test in coat proteins and between coat and nonCE protein

##########################################
rm(list=ls());
 
## some test to demonstrate the syntax for alternative hypothesis
x = rnorm(100)+5
y = x + 2
#summary(x); 
#summary(y)
wilcox.test( x, y, alt='gr') #accept null, large p 
wilcox.test( y, x, alt='gr') #accept alter, small p

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
coat11.tb = read.delim("_free.omega.coat.paml.Jan20,2011.txt", header=F);
coat.tb = read.delim("_free.omega.GBlocked.coat.paml.Jan17,2013.txt", header=F); #20130123
coat.tb = coat.tb[ ! is.na(coat.tb$V4), ]
coat.tb$V2 = as.character(coat.tb$V2) 
coat.tb$V1 = as.character(coat.tb$V1) 

positions = match( coat.tb$V1, rownames(ctb3))
#ctb3[positions[1:5], ] #check, passed
coat.tb = cbind( coat.tb, ctb3[positions, c('cereushits','subtilishits', 'outhits', 'hits')])
#head(coat.tb)

positions = match( coat.tb$V1, ctb$id)
coat.tb$CoatLocation = ctb$CoatLocation[positions]
coat.tb$node = "internal"
coat.tb$node[grep("^B", coat.tb$V2)] = "leaf"
head(coat.tb)
unique(coat.tb$CoatLocation)


### nonCE free omega
nonCE11.tb = read.delim("_nonCE.omega.7March2011.tab", header=F); #2011 results, without GBlocks
nonCE11.tb = nonCE11.tb[nonCE11.tb$V2=='H1C',] # H1C is where free omega are calculated 
head(nonCE11.tb)

nonCE.tb = read.delim("_nonCE.omega.2013Jan23.tab", header=F); #20130123, after GBlocks
nonCE.tb$V1 = as.character( nonCE.tb$V1 )
nonCE.tb$V2 = as.character( nonCE.tb$V2 )
nonCE.tb = nonCE.tb[nonCE.tb$V2=='H1C',] # H1C is where free omega are calculated 
#nonCE.tb = nonCE.tb[- grep("#", nonCE.tb$V1) , ]
nonCE.tb$V3 = as.character(nonCE.tb$V3)
nonCE.tb$node = "internal"
nonCE.tb$node[grep("^B", nonCE.tb$V3)] = "leaf"
head(nonCE.tb)
str(nonCE.tb)

#compare nonCE w be and afer GBlocks
summary(nonCE11.tb$V11)
summary(nonCE.tb$V11)
t.test( nonCE11.tb$V11, nonCE.tb$V11)
hist(nonCE11.tb$V11)
hist(nonCE.tb$V11)

#how much w > 1, 999 
quantile( nonCE.tb$V11, prob=c(0.2, 0.5, 0.8,0.83, 0.85, 0.86, 0.87, 0.88, 0.9))
quantile( coat.tb$V10, prob=c(0.2, 0.5, 0.8,0.815, 0.82, 0.83, 0.85, 0.86, 0.87, 0.88, 0.9))

#### parse the omega results 
clades = unique(c( nonCE.tb$V3, coat.tb$V2))  #35 clades

BsuClade0 = clades[ grep( "Bsu|Bpu|Bam|Bli|Bmo", clades)]; #20130123
BsuClade = BsuClade0[- grep("Bce|Bwe|Ban|Bth", BsuClade0)]

BceClade0 = clades[grep("Bce|Bwe|Ban|Bth", clades)]
BceClade = BceClade0[ - grep( "Bsu|Bpu|Bam|Bli|Bmo", BceClade0) ]

BsuClade11 = c("(Bpu,Bli,Bam,Bsu)","(Bli,Bam,Bsu)","(Bpu,Bli,Bsu)", "Bpu", "Bli", "(Bam,Bsu)", "Bam", "Bsu", "(Bpu,Bam,Bsu)", 
             "(Bpu,Bli,Bam,Bmo,Bsu)", "(Bli,Bam,Bmo,Bsu)", "(Bam,Bmo,Bsu)","(Bmo,Bsu)","Bmo","(Bli,Bam,Bsu,Bmo)", "(Bam,Bsu,Bmo)", "(Bsu,Bmo)", 
             "(Bpu,Bam,Bmo,Bsu)", "(Bli,Bmo,Bsu)", "(Bli,Bsu)" );   
BceClade11 = c("(Bwe,Bce,Ban,Bth)","(Bce,Ban,Bth)","(Ban,Bth)","(Bwe,Bce,Bth)","(Bce,Bth)", "(Bwe,Ban,Bth)",
               "Bwe","Ban","Bth","Bce");
match(sort(BceClade), sort(BceClade11)) #20130123, 2013 and 2011 gave exactly the same branchs. YES.
match(sort(BsuClade), sort(BsuClade11))


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

table( coat.tb$V2, coat.tb$flag) #check flags
table(coat.tb$node, coat.tb$flag) #check flags

#2013 for updated table 
head(coat.tb)
head(nonCE.tb)
length(unique(nonCE.tb$V1)) #2085 genes
length(unique(coat.tb$V1)) #47 genes

# remove overly high omegas
# how many?? 
Wcutoff = 10
coat.tb$V10[coat.tb$V10> Wcutoff ] = NA
nonCE.tb$V11[nonCE.tb$V11> Wcutoff ] = NA

### leaf node in Bsu clade
BsuClade
summary( coat.tb$V10[coat.tb$node=='leaf' & coat.tb$flag=='bsu'] )
summary( nonCE.tb$V11[nonCE.tb$node=='leaf' & nonCE.tb$flag=='bsu'] )
wilcox.test( coat.tb$V10[coat.tb$node=='leaf' & coat.tb$flag=='bsu'], 
             nonCE.tb$V11[nonCE.tb$node=='leaf' & nonCE.tb$flag=='bsu'],
             alternative='gr') #

summary( coat.tb$V10[coat.tb$node=='leaf' & coat.tb$flag=='bsu' & coat.tb$V10<1 ] )
summary( nonCE.tb$V11[nonCE.tb$node=='leaf' & nonCE.tb$flag=='bsu' & nonCE.tb$V11<1 ] )
wilcox.test( coat.tb$V10[coat.tb$node=='leaf' & coat.tb$flag=='bsu' & coat.tb$V10<1 ], 
             nonCE.tb$V11[nonCE.tb$node=='leaf' & nonCE.tb$flag=='bsu' & nonCE.tb$V11<1],
             alternative='gr') #
head(coat.tb[(coat.tb$node=='leaf' & coat.tb$flag=='bsu' & coat.tb$V10<1 ), ] )
head(nonCE.tb[nonCE.tb$node=='leaf' & nonCE.tb$flag=='bsu' & nonCE.tb$V11<1,])


summary( coat.tb$V10[coat.tb$node=='leaf' & coat.tb$flag=='bsu' & coat.tb$V10 > 1 ] )
summary( nonCE.tb$V11[nonCE.tb$node=='leaf' & nonCE.tb$flag=='bsu' & nonCE.tb$V11 > 1 ] )
wilcox.test( coat.tb$V10[coat.tb$node=='leaf' & coat.tb$flag=='bsu' & coat.tb$V10 > 1 ], 
             nonCE.tb$V11[nonCE.tb$node=='leaf' & nonCE.tb$flag=='bsu' & nonCE.tb$V11 > 1],
             alternative='gr') #
head(coat.tb[(coat.tb$node=='leaf' & coat.tb$flag=='bsu' & coat.tb$V10 > 1 ), ] )
head(nonCE.tb[nonCE.tb$node=='leaf' & nonCE.tb$flag=='bsu' & nonCE.tb$V11 > 1,])

##################
### internal node in BsuClade
nodechoice = 'internal'; cladechoice = 'bsu';
BsuClade
head(coat.tb[(coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10<1 ), ] )
head(nonCE.tb[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11<1,])

summary( nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice] )
summary( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice] )
wilcox.test( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice], 
             nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice],
             alternative='gr') #

summary( nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice  & nonCE.tb$V11<1 ] )
summary( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10<1 ] )
wilcox.test( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice  & coat.tb$V10<1 ], 
             nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11<1],
             alternative='gr') #


head(nonCE.tb[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11 > 1,])
head(coat.tb[(coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10 > 1 ), ] )

summary( nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11 > 1 ] )
summary( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10 > 1 ] )
wilcox.test( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10 > 1 ], 
             nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11 > 1],
             alternative='gr') 

##################
### leaf node in BceClade
nodechoice = 'leaf'; cladechoice = 'bce';
BceClade
head(coat.tb[(coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10<1 ), ] )
head(nonCE.tb[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11<1,])

summary( nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice] )
summary( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice] )
wilcox.test( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice], 
             nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice],
             alternative='gr') #

summary( nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice  & nonCE.tb$V11<1 ] )
summary( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10<1 ] )
wilcox.test( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice  & coat.tb$V10<1 ], 
             nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11<1],
             alternative='gr') #


head(nonCE.tb[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11 > 1,])
head(coat.tb[(coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10 > 1 ), ] )

summary( nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11 > 1 ] )
summary( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10 > 1 ] )
wilcox.test( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10 > 1 ], 
             nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11 > 1],
             alternative='gr') 


##################
### internal node in BceClade
nodechoice = 'internal'; cladechoice = 'bce';
BceClade
head(coat.tb[(coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10<1 ), ] )
head(nonCE.tb[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11<1,])

summary( nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice] )
summary( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice] )
wilcox.test( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice], 
             nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice],
             alternative='gr') #

summary( nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice  & nonCE.tb$V11<1 ] )
summary( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10<1 ] )
wilcox.test( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice  & coat.tb$V10<1 ], 
             nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11<1],
             alternative='gr') #

head(nonCE.tb[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11 > 1,], 10)
head(coat.tb[(coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10 > 1 ), ], 20 ) #none!
hist(coat.tb$V10[(coat.tb$node==nodechoice & coat.tb$flag==cladechoice)])

summary( nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11 > 1 ] )
summary( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10 > 1 ] )
#wilcox.test( coat.tb$V10[coat.tb$node==nodechoice & coat.tb$flag==cladechoice & coat.tb$V10 > 1 ], 
#             nonCE.tb$V11[nonCE.tb$node==nodechoice & nonCE.tb$flag==cladechoice & nonCE.tb$V11 > 1],
#             alternative='gr') 


quit("yes")

#########
#
# END of 2013 Jan 23 analysis
#
########

BceClade
wilcox.test( coat.tb$V10[coat.tb$node=='leaf' & coat.tb$flag=='bce'], 
             nonCE.tb$V11[nonCE.tb$node=='leaf' & nonCE.tb$flag=='bce'],
             alternative='gr') #
summary( coat.tb$V10[coat.tb$node=='leaf' & coat.tb$flag=='bce'] )
summary( nonCE.tb$V11[nonCE.tb$node=='leaf' & nonCE.tb$flag=='bce'] )

coat.tb[coat.tb$node=='internal' & coat.tb$flag=='bsu' & coat.tb$V10<1, ]


######end 2013Jan23

#leaf branches for publication  ********  
#BsuClade = c("Bpu", "Bli", "Bam", "Bsu");   #2011 p-value = 0.00823

# internal branches for publication  ********
### BsuClade = c("(Bpu,Bli,Bam,Bsu)","(Bli,Bam,Bsu)","(Bpu,Bli,Bsu)", "(Bam,Bsu)", "(Bpu,Bam,Bsu)", "(Bpu,Bli,Bam,Bmo,Bsu)", "(Bli,Bam,Bmo,Bsu)", "(Bam,Bmo,Bsu)","(Bmo,Bsu)","(Bli,Bam,Bsu,Bmo)", "(Bam,Bsu,Bmo)", "(Bsu,Bmo)", "(Bpu,Bam,Bmo,Bsu)", "(Bli,Bmo,Bsu)", "(Bli,Bsu)" );   
#no leaf nodes, high omega all over the places, large p-value
#so, these suggest coat gene only recently evolves faster than the nonCE genes, suggest adpation to niches

#BsuClade = c("(Bpu,Bli,Bam,Bsu)","(Bli,Bam,Bsu)","(Bpu,Bli,Bsu)", "Bpu", "Bli", "(Bam,Bsu)", "Bam", "(Bpu,Bam,Bsu)", "(Bpu,Bli,Bam,Bmo,Bsu)", "(Bli,Bam,Bmo,Bsu)", "(Bam,Bmo,Bsu)","(Bmo,Bsu)","Bmo","(Bli,Bam,Bsu,Bmo)", "(Bam,Bsu,Bmo)", "(Bsu,Bmo)", "(Bpu,Bam,Bmo,Bsu)", "(Bli,Bmo,Bsu)", "(Bli,Bsu)" );   
#no Bsu leaf, p=0.13

#BsuClade = c("(Bli,Bam,Bsu)","Bli","(Bam,Bsu)","Bam","Bsu","(Bpu,Bam,Bsu)", "(Bli,Bam,Bmo,Bsu)", "(Bam,Bmo,Bsu)", "(Bmo,Bsu)","Bmo","(Bli,Bam,Bsu,Bmo)", "(Bam,Bsu,Bmo)","(Bsu,Bmo)","(Bli,Bmo,Bsu)","(Bli,Bsu)" );  
# p-value = 0.01123, no Bpu


#BsuClade = c("Bli", "Bam", "Bsu");   #0.01082
#BsuClade = c("Bpu", "Bam", "Bsu");   # p-value = 0.002779
#** BsuClade = c("Bsu") # p-value = 0.002939
#**BsuClade = c("Bam") # p=0.0858
#BsuClade = c("Bsu", "Bmo", "(Bsu,Bmo)","(Bmo,Bsu)") # p-value = 3.483e-08
#BsuClade = c("Bsu", "Bmo") #  p-value = 4.654e-06
#BsuClade = c("Bsu", "(Bsu,Bmo)","(Bmo,Bsu)") # p-value = 1.899e-05

#BsuClade = c("Bli") # p=0.5 ??
#BsuClade = c("Bpu") # p=0.19
#BsuClade = c("(Bam,Bmo,Bsu)", "(Bsu,Bmo)", "(Bam,Bsu)") #large p
#BsuClade = c( "(Bam,Bsu)") #p=0.1568
     

# internal branches for publication  ********
BceClade11 = c("(Bwe,Bce,Ban,Bth)","(Bce,Ban,Bth)","(Ban,Bth)","(Bwe,Bce,Bth)","(Bce,Bth)", "(Bwe,Ban,Bth)"); # p=?  cCW>nCW

# leaf brnaches for publication, 20130123
BceClade11 = c("Bwe","Ban","Bth","Bce"); # 2011, p=0.4, 0.79 for w>1, 0.11 for w< 1 cCW>nCW

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

table( coat.tb$V2, coat.tb$flag) #check flags
table(coat.tb$node, coat.tb$flag) #check flags

nonCE.tb$omegaFlag = ifelse(nonCE.tb$V11>1, 'high','low' )
coat.tb$omegaFlag = ifelse(coat.tb$V10>1, 'high','low' )

table( nonCE.tb$flag, nonCE.tb$omegaFlag )
table( coat.tb$flag,  coat.tb$omegaFlag )
table( coat.tb$V2, coat.tb$omegaFlag)
table( coat.tb$flag, coat.tb$V2, coat.tb$omegaFlag)


tmptb = table( coat.tb$V1, coat.tb$omegaFlag)
tmptb = tmptb[-1, ]
coat.summary.tb = as.data.frame( tmptb )
coat.summary.tb$highPer = coat.summary.tb$high / ( coat.summary.tb$low + coat.summary.tb$high) 

positions = match( rownames(ctb3), rownames(coat.summary.tb) )
ctb3 = cbind( ctb3, coat.summary.tb[positions, ])
plot(ctb3$highPer ~ jitter(ctb3$cereushits) )
summary(lm(ctb3$highPer ~ ctb3$cereushits)) #significant, conserved coat genes evolve have more high omegas? how about bsu- bce- clades? 

#summary(lm( nonCE.tb$V11 ~ nonCE.tb$flag ))
#summary(lm( coat.tb$V10 ~ coat.tb$flag + coat.tb$hits ))
#summary(lm( coat.tb$V10 ~ coat.tb$hits )) #not good
#summary(lm( coat.tb$V10 ~ coat.tb$subtilishits )) # terrible
#summary(lm( log(coat.tb$V10) ~ coat.tb$cereushits )) #p=0.3488
summary(lm( coat.tb$V10 ~ coat.tb$cereushits )) #p=0.0037
plot( coat.tb$V10 ~ jitter(coat.tb$cereushits) )

#partition by omegas, Nov 29, 2011
#n.high = nonCE.tb[nonCE.tb$V11>1 & nonCE.tb$V11<100, ]
n.high = nonCE.tb[nonCE.tb$V11>1 , ]
n.low = nonCE.tb[nonCE.tb$V11<1, ]
#c.high = coat.tb[coat.tb$V10>1 & coat.tb$V10<100, ]
c.high = coat.tb[coat.tb$V10>1 , ]
c.low = coat.tb[coat.tb$V10<1, ]

table(n.high$V3)
table(n.low$V3)
table(c.high$V2)
table(c.low$V2)

FlagPartitionOmega = 2; #1 for w > 1, 2 for w<1
if( FlagPartitionOmega ==1 ){
  nonCE.tb = n.high;
  coat.tb = c.high; 
}
if( FlagPartitionOmega ==2 ){
  nonCE.tb = n.low;
  coat.tb = c.low; 
}


nBW =  nonCE.tb$V11[nonCE.tb$flag=='bsu'];
nCW =  nonCE.tb$V11[nonCE.tb$flag=='bce'];
cBW =  coat.tb$V10[coat.tb$flag=='bsu']; 
cCW =  coat.tb$V10[coat.tb$flag=='bce']; 
cBW.conserved =  coat.tb$V10[coat.tb$flag=='bsu' & coat.tb$cereushits>=4]; 
cCW.conserved =  coat.tb$V10[coat.tb$flag=='bce' & coat.tb$cereushits>=4]; 
cBW.labile =  coat.tb$V10[coat.tb$flag=='bsu' & coat.tb$cereushits <=3 ]; 
cCW.labile =  coat.tb$V10[coat.tb$flag=='bce' & coat.tb$cereushits <=3 ]; 

####################### summary

summary(nBW)
summary(cBW)
summary(nCW)
summary(cCW)

summary(cBW.conserved)
summary(cBW.labile)
summary(cCW.conserved)
summary(cCW.labile)


### analysis of Bsu clade (Because I can change Bsu clade, this gave me opportunity to compare internal branches and leap nodes. )
BsuClade
summary(nBW)
summary(cBW)
summary(cBW.conserved)
summary(cBW.labile)

BsuClade
wilcox.test( cBW, nBW, al='gr') #
t.test( log(cBW), log(nBW), al='gr') 
#ks.test( cBW, nBW, al='gr')
#wilcox.test( cBW.conserved, nBW, al='less') #
#wilcox.test( cBW.labile, nBW, al='gr') #

### analysis of Bce clade
summary(nCW)
summary(cCW)
summary(cCW.conserved)
summary(cCW.labile)

wilcox.test( cCW, nCW, al='gr')
t.test( log(cCW), log(nCW), al='gr'  )
#t.test( cCW, nCW, al='gr' )
#ks.test( cCW, nCW, al="less" )

#wilcox.test( cCW.conserved, nCW, al='less')
#wilcox.test( cCW.labile, nCW, al='gr')

#wilcox.test( log(cCW), log(nCW), al='gr'  ) #for log transformation make no difference for wilcox test
#q("no")

unique(n.high$V1[n.high$V11>100]) #2013 Jan 16

