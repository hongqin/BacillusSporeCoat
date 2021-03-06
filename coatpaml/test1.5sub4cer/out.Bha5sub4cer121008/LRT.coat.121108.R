#LRT use 2 *(loglikelihood1 - loglikelihood2)

#my @Htres = ( 
#1#"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,(Ban,Bth))));", #H0
#2#"(Bha,(Bpu,(Bli,(Bam ,(Bmo #1,Bsu #1)#1))),(Bwe,(Bce,(Ban,Bth))));", #H1a
#3#"(Bha,(Bpu,(Bli,(Bam #1,(Bmo #1,Bsu #1)#1)#1)),(Bwe,(Bce,(Ban,Bth))));", #H1b
#4#"(Bha,(Bpu,(Bli #1,(Bam #1,(Bmo #1,Bsu #1)#1)#1)#1),(Bwe,(Bce,(Ban,Bth))));", #H1c
#5#"(Bha,(Bpu #1,(Bli #1,(Bam #1,(Bmo #1,Bsu #1)#1)#1)#1)#1,(Bwe,(Bce,(Ban,Bth))));", #H1d
#6#"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,(Ban #1,Bth #1)#1)));", #H2a
#7#"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce #1,(Ban #1,Bth #1)#1)#1));", #H2b
#8#"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #1,(Bce #1,(Ban #1,Bth #1 )#1)#1));", #H2c
#9#"(Bha,(Bpu #2,(Bli #2,(Bam #2,(Bmo #2,Bsu #2)#2)#2)#2)#2,(Bwe #1,(Bce #1,(Ban #1,Bth #1)#1)#1)#1);", #H3a
#10#"(Bha,(Bpu ,(Bli,(Bam ,(Bmo ,Bsu )) ) ),                  (Bwe #2,(Bce #2,(Ban #1,Bth #2)#2)#2)#2 );", #H3c
#11#"(Bha,(Bpu #3,(Bli#3, (Bam #3,(Bmo #3 ,Bsu #3)#3)#3)#3)#3,(Bwe #2,(Bce #2,(Ban #1,Bth #2)#2)#2)#2 );" #H4a
#);

 tb = read.table( "out.Bha5sub4cer121108.csv", header=F
 #   colClasses=c("chr","chr","int","chr","chr","num","num")
  );
 tb[,1] = as.character( tb[,1] );

 bgs = unique( tb[,1]);

 models = c("H0","H1a","H1b","H1c","H1d", "H2a","H2b","H2c","H3a","H3c","H4a");
 # labels1 = c("bg","p.H0.H1a", "p.H0.H1b","p.H0.H1c","p.H0.H1d","p.H0.H2a","p.H0.H2b","p.H0.H2c"); 
 labels1 = c("bg",paste( "p","H0", models[2:8], sep=".") ); labels1;
 out1 = data.frame( matrix(nrow=length(bgs), ncol=length(labels1) )); 
 names( out1 ) = labels1;
 out1$bg = bgs;

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
 
 for( j in 1:length(bgs) ) {
  bg = bgs[j];
  sub = tb[ tb[,1] == bg,]
  for( i in 2:8 ) { #i is column in out1, row in sub
    out1[ out1$bg==bg,i] = round( pchisq( 2*(sub[i,6] - sub[1,6]), df=1,low=F), 4); ##??
  }
  sub[,c(2,6)]
  sub[,c(2,6)]

  col=2;
  for( i in c(3,4,5,7,8) ) {  #c("H1b","H1c","H1d","H2b","H2c")
    out2[ out2$bg==bg, col ] = round( pchisq(  2*(-sub[i,6] + sub[9, 6]), df=2,low=F), 4); ##H3a
    out2c[ out2c$bg==bg, col ] = round( pchisq( 2*(-sub[i,6] + sub[10,6]), df=2,low=F), 4); ##H3c
    col= col+1;
  }
  sub[,c(2,6)]
  out2[j,]

  col=2;
  for( i in 9:10){ # "p.H3c.H4a"
   out3[ out3$bg==bg, col] = round( pchisq( 2*(-sub[i,6] + sub[11,6]), df=1,low=F), 4); ## H4a
   #out3[ out3$bg==bg, col] =  2*(-sub[i,6] + sub[11,6])  ## H4a
   col = col+1;
  }

  cutoff = 0.05
  ps = (out1[j,2:8] < cutoff)
  x = ps[ps==TRUE]
  out1$sig1[j] = length(x);
 # out1$sig[j] = ifelse( min(out1[j,2:8])<cutoff, 1 ,0 )
}

out.all = cbind( out1, out2, out2c, out3);
write.csv(out.all, "_LRT.coat.Bha5sub4cer.121108.all.csv", quote=F, row.names=F)

cutoff = 0.05
#mycol = 1:length(out.all[1,]);
#mycol2 = mycol[-c(1,9,10,16,22) ];
#names(out.all)[mycol2];

out.all = cbind( out1[,2:9], out2[,2:length(out2[1,])], out2c[,2:length(out2c[1,])], out3[,2:length(out3[1,])]);
rownames(out.all) = out1$bg;

mycols = 1:length(out.all[1,]);

for ( mycol in mycols[-8] ) {
  for( myrow in 1:length(out.all[,1]) ){
    out.all[myrow,mycol] = ifelse( out.all[myrow,mycol]< cutoff, out.all[myrow,mycol], NA);
  }
}


for ( mycol in mycols ) {
    tmp = out.all[,mycol];
    tmp[is.na(tmp)] = '';
    out.all[,mycol] = tmp;
}

write.csv(out.all, "_LRT.coat.Bha5sub4cer.121108.p0.05.csv", quote=F, row.names=T)

q("yes");
### END
###############################
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

