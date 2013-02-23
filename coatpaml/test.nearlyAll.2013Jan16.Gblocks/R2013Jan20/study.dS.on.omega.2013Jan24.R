# 2013 Jan 24, how does dS affect omega? 

# 2013 Jan 20, How does GBLOCKs affect omega, does it remove extremely high omegas, such as 9999?
rm(list=ls())

newtb = read.table("_free.omega.GBlocked.coat.paml.Jan17,2013.txt", sep="\t", header=F)
oldtb = read.table("_free.omega.coat.paml.Jan20,2011.txt", sep="\t", header=F)
head(newtb)
head(oldtb)
summary(newtb$V10)
summary(oldtb$V10) #omega

summary(newtb$V12) #dN
newtb$V14[newtb$V12==0] = NA

summary(newtb$V14) #dS
newtb$V14[newtb$V14==0] = NA

#omega ~ dS
summary(lm(newtb$V10 ~ newtb$V14))
plot( newtb$V10 ~ newtb$V14)

#omega ~ dS
summary(lm(log10(newtb$V10) ~ log10(newtb$V14)))  
plot( log10(newtb$V10) ~ log10(newtb$V14))

#summary(lm(log10(oldtb$V10) ~ log10(oldtb$V14)))  
plot( log10(oldtb$V10) ~ log10(oldtb$V14))

#omega ~ dN
summary(lm(log10(newtb$V10) ~ log10(newtb$V12)))
plot( log10(newtb$V10) ~ log10(newtb$V12))
plot( log2(newtb$V10) ~ log2(newtb$V12))

plot( log2(oldtb$V10) ~ log2(oldtb$V12))
#old and new data show similar patterns, removing gapp did not have a major effect on omega estimation? 


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

head(nonCE.tb)
#nonCE omega ~ dS, negative correlation
plot( log2(nonCE.tb$V11) ~ log2(nonCE.tb$V15))

#nonCE omega ~ dN, ?correlation
plot( log2(nonCE.tb$V11) ~ log2(nonCE.tb$V13))



#
# end of dS study on 2013 Jan 24
#

summary(lm(newtb$V10 ~ oldtb$V10))
#R=0.158, weak correlation? p<<0.001
plot( newtb$V10 ~ oldtb$V10 )

summary(lm(log10(newtb$V10) ~ log10(oldtb$V10)))
#R=0.307, weak correlation? p<<0.001
plot( log10(newtb$V10) ~ log10(oldtb$V10) )

#merge the two tables
newtb$BGbranch = paste( newtb[,1], newtb[,2] )
oldtb$BGbranch = paste( oldtb[,1], oldtb[,2] )
names1 = names(newtb)
names1[10]= "w.new"
names2 = names(oldtb)
names2[10]= "w.old"
names(newtb) = names1
names(oldtb) = names2


########################
#how many branch-omega are influenced by alignment, as indicated by GBlocks? 
new.old.tb = merge(newtb[,c("V1","V2","w.new","BGbranch")], oldtb[,c("V1","V2","w.old","BGbranch")])
new.old.tb$new.vs.old = new.old.tb$w.new / new.old.tb$w.old
summary(new.old.tb)
hist( new.old.tb$new.vs.old[ new.old.tb$new.vs.old<3  ], br=100 )
myquants = quantile(new.old.tb$new.vs.old, prob= c(0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))

new.old.tb$flag = 'Different'
for( i in 1:length( new.old.tb[,1] )) {
  if (( new.old.tb$new.vs.old[i] > 1/1.25) & (new.old.tb$new.vs.old[i]< 1.25)) {
    new.old.tb$flag[i] = "Consistent"
  }
}

table( new.old.tb$flag )
table( new.old.tb$V1, new.old.tb$flag )
table( new.old.tb$V2, new.old.tb$flag )

consistent.tb = new.old.tb[ new.old.tb$flag=='Consistent',  ]
different.tb = new.old.tb[ new.old.tb$flag=='Different',  ]
summary(consisten.tb)
summary(different.tb)
plot( consistent.tb$w.new ~ consistent.tb$w.old )
plot( different.tb$w.new ~ different.tb$w.old )
plot( log10(consistent.tb$w.new) ~ log10(consistent.tb$w.old ))
plot( log10(different.tb$w.new) ~ log10(different.tb$w.old ) )
quantile(consistent.tb$w.new, prob= c(0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9, 0.91, 0.92, 0.95, 0.99))
length( consistent.tb$w.new[ consistent.tb$w.new >1  ] ) / length( consistent.tb[,1])
#8.2% has omega.new > 1



                   