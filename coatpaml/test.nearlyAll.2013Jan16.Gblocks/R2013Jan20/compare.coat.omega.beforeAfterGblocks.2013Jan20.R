# 2013 Jan 20, How does GBLOCKs affect omega, does it remove extremely high omegas, such as 9999?

newtb = read.table("_free.omega.GBlocked.coat.paml.Jan17,2013.txt", sep="\t", header=F)
oldtb = read.table("_free.omega.coat.paml.Jan20,2011.txt", sep="\t", header=F)
head(newtb)
head(oldtb)
summary(newtb$V10)
summary(oldtb$V10)

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

#how many w > 10? 
head(newtb)
tmp = newtb$w.new
length( tmp[tmp>10]) / length(tmp)

  
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



                   