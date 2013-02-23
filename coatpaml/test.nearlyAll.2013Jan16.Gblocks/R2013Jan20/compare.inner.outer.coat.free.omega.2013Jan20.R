rm(list=ls());

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
#coat.tb = read.delim("_free.omega.coat.paml.Jan20,2011.txt", header=F);
coat.tb = read.table("_free.omega.GBlocked.coat.paml.Jan17,2013.txt", header=F);
coat.tb = coat.tb[ ! is.na(coat.tb$V4), ]
coat.tb$V2 = as.character(coat.tb$V2) 

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
wilcox.test( coat.tb$V10[coat.tb$CoatLocation=="InnerCoat"], coat.tb$V10[coat.tb$CoatLocation=="OuterCoat"],  )
#p=0.81
t.test( coat.tb$V10[coat.tb$CoatLocation=="InnerCoat"], coat.tb$V10[coat.tb$CoatLocation=="OuterCoat"]  )
#p=0.16
t.test( log10(coat.tb$V10[coat.tb$CoatLocation=="InnerCoat"]), log10(coat.tb$V10[coat.tb$CoatLocation=="OuterCoat"])  )
#p=0.46

t.test( log10(coat.tb$V10[coat.tb$CoatLocation=="InnerCoat" &coat.tb$node=='leaf' ]), 
        log10(coat.tb$V10[coat.tb$CoatLocation=="OuterCoat"& coat.tb$node=='leaf' ])  )

t.test( log10(coat.tb$V10[coat.tb$CoatLocation=="InnerCoat" &coat.tb$node=='internal' ]), 
        log10(coat.tb$V10[coat.tb$CoatLocation=="OuterCoat"& coat.tb$node=='internal' ])  )

wilcox.test( log10(coat.tb$V10[coat.tb$CoatLocation=="InnerCoat" &coat.tb$node=='leaf' ]), 
        log10(coat.tb$V10[coat.tb$CoatLocation=="OuterCoat"& coat.tb$node=='leaf' ])  )

t.test( log10(coat.tb$V10[coat.tb$CoatLocation=="InnerCoat" &coat.tb$node=='internal' ]), 
        log10(coat.tb$V10[coat.tb$CoatLocation=="OuterCoat"& coat.tb$node=='internal' ])  )



