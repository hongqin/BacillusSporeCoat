
new = read.csv("_profile.032410.csv", sep="\t")
old = read.csv("_profile.073108b.csv", sep="\t")
new[,1] = as.character(new[,1])
old[,1] = as.character(old[,1])

tmp = cbind( new[,1], old[,1])
table( table(tmp) ) #all matched

new=new[,c("BGcoat","Bce79", "Bce87", "Bce3L")]
old=old[,c("id","Bce")]

table( new$Bce79 == old$Bce ) #47 hits

table( table(new[,2])) #47 hits
table( table(new[,3])) #46 hits
table( table(new[,4])) #45 hits

new2 = new
new[,2] = ifelse(is.na(new[,2]), 0, 1)
new[,3] = ifelse(is.na(new[,3]), 0, 1)
new[,4] = ifelse(is.na(new[,4]), 0, 1)
new$bce.sum = apply( new[,2:4], 1, sum)

profile = read.csv("_coat.profile.summary.March22,2010.csv");
profile[,1] = as.character( profile[,1] )

profile$id == new$BGcoat #not matched due to duplications

new$CoatLocation = profile$CoatLocation[  match(new$BGcoat, profile$id) ]
tb1 = table( new$CoatLocation, new$bce.sum )
chisq.test( tb1[c(2,6), c(1,4)]  )
fisher.test( tb1[c(2,6), c(1,4)]  )

new = data.frame( cbind( new2, new ) );
write.csv(new, "_new.temp.hits.csv")




