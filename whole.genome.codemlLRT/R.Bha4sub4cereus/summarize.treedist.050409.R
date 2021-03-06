
 tretb = read.table( "_BGids.treedist.result.csv", sep="\t" );
 names(tretb) = c("id","d");
 table(tretb$d);
#  0   2   4   6   8  10  12 
#978 292 106  27   7   1   2 

 out0 = tretb[tretb$d==0,]
 write.table( out0, "_BGids.treedist.is.0.050409.csv", row.names=F, quote=F, col.names=F );

 out02 = tretb[tretb$d<=2,]
 write.table( out02, "_BGids.treedist.lessthan2.050409.csv", row.names=F, quote=F, col.names=F );

 ###coat bg
 tbc = read.csv( "_coat.profile.122908.csv", sep="\t");
 tbc$id = as.character(tbc$id);
 tbc$type = "coat";
 rownames(tbc) = tbc$id;

 intersect( out0$id, tbc$id); #15 genes
 intersect( out02$id, tbc$id); #19 genes





