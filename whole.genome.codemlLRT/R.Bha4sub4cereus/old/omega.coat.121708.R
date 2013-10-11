
 rm(list=ls());

 tb = read.csv("_LRT.coat.4sub4cer.121708.all.csv", colClasses= c("character","numeric") )
 sig = tb[tb[,2]<0.05,]

 omega = read.table( "__omega.All.Ban.121708.csv", sep="\t" );
 omega[,1] = as.character(omega[,1]);

 wilcox.test( omega$V4, omega$V6); # cereus actually have higher omega

 sub = omega[match(sig$bg, omega[,1]),]
 wilcox.test(sub$V4,sub$V6); #cereus still have higher omega


q("yes");

###############################
### END  ###
###############################

 labels2  = c("bg",paste( "p", c("H1b","H1c","H1d","H2b","H2c"),"H3a",sep=".") );
 labels2c = c("bg",paste( "p", c("H1b","H1c","H1d","H2b","H2c"),"H3c",sep=".") );
 out2  = data.frame( matrix(nrow=length(bgs), ncol=length(labels2) )); 
 out2c = data.frame( matrix(nrow=length(bgs), ncol=length(labels2c) )); 
 names( out2)  = labels2;
 names( out2c) = labels2c;
 out2$bg  = bgs;
 out2c$bg = bgs;
 
 labels3 = c("bg",paste( "p", c("H3a","H3c"), "H4a",sep=".") );
 out3 = data.frame( matrix(nrow=length(bgs), ncol=length(labels3) )); 
 names( out3) = labels3;
 out3$bg = bgs;





out2 = out;

for( row in 1:length(out2[,1]) ){
  for ( col in 2:16 ) {
    out2[row,col] = ifelse( out[row,col]< 0.05, out[row,col], NA);
  }
  tmp1 = ifelse( (out[row,"p.H0.H1c"] < 0.05) & (out[row,"p.H1c.H3a"] < 0.05 ), 1, 0)
  tmp2 = ifelse( (out[row,"p.H0.H2c"] < 0.05) & (out[row,"p.H2c.H3a"] < 0.05 ), 1, 0)
  tmp3 = ifelse( (out[row,"p.H0.H2c"] < 0.05) & (out[row,"p.H2c.H3b"] < 0.05 ), 1, 0)
  tmp = c(tmp1, tmp2, tmp3);
  out2$flag[row] = sum( tmp );
}

out3 = out2[out2$flag>0,]
write.csv(out3, "chisq.coat.Bha5sub4cer.091508.small.csv", quote=F, row.names=F)

out2 = out;
for( row in 1:length(out2[,1]) ){
  for ( col in 2:16 ) {
    out2[row,col] = ifelse( out[row,col]< 0.05, out[row,col], '');
  }
}

out2;

#idtab = read.table("/home/hqin/coat.protein07/key.data/BGcoat.csv", header=T, sep="\t");
idtab = read.csv("/home/hqin/coat.protein07/key.data/BGcoat3.csv", header=T);

idtab$NCBI = as.character( idtab$NCBI );
idtab$Gene = as.character( idtab$Gene );
out2$name = idtab$Gene[ match( out2$bg, idtab$NCBI) ];
out2$CoatLocation = idtab$CoatLocation[ match( out2$bg, idtab$NCBI) ];

write.csv(out2, "chisq.coat.Bha5sub4cer.091508.csv")

for( row in 1:length(out.all[,1]) ){
  for ( col in 2:(length(out1[1,])-1) ) {
    #out1[row,col] = ifelse( out1[row,col]< 0.05, out1[row,col], NA);
    out1[row,col] = ifelse( out1[row,col]< 0.05, out1[row,col], '');
  }
  for ( col in 2:length(out2[1,]) ) {
    out2[row,col]  = ifelse( out2 [row,col]< 0.05, out2 [row,col], '');
    out2c[row,col] = ifelse( out2c[row,col]< 0.05, out2c[row,col], '');
  }
  for ( col in 2:length(out3[1,]) ) {
    out3[row,col] = ifelse( out3[row,col]< 0.05, out3[row,col], '');
  }
}

