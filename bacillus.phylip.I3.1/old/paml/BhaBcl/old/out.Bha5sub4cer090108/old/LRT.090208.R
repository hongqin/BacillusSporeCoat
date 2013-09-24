#LRT use 2 *(loglikelihood1 - loglikelihood2)

 tb = read.table( "out.Bha5sub4cer090108.csv", header=F
 #   colClasses=c("chr","chr","int","chr","chr","num","num")
  );
 tb[,1] = as.character( tb[,1] );

 bgs = unique( tb[,1]);

 labels = c("bg","p.H0.H1a", "p.H0.H1b","p.H0.H1c","p.H0.H1d","p.H0.H2a","p.H0.H2b","p.H0.H2c", "p.H1c.H3a","pH2c.H3a","p.H2c.H3b" );
 out = data.frame( matrix(nrow=length(bgs), ncol=length(labels) )); 
 names( out ) = labels;
 out$bg = bgs;

 for( j in 1:length(bgs) ) {
  bg = bgs[j];
  sub = tb[ tb[,1] == bg,]
  for( i in 2:8 ) {
   # out[ out$bg==bg,i] = round( pchisq( (sub[i,6] - sub[1,6]), df=1,low=F), 3); ##error here?
    out[ out$bg==bg,i] = round( pchisq( 2*(sub[i,6] - sub[1,6]), df=1,low=F), 3); ##??
  }

  i=9 #"p.H1c.H3a",
  out[ out$bg==bg,i] = round( pchisq( 2*(sub[i,6] - sub[4,6]), df=1,low=F), 3); ##??
  i=9  #"pH2c.H3a",
  out[ out$bg==bg,i+1] = round( pchisq( 2*(sub[i,6] - sub[8,6]), df=1,low=F), 3); ##??
  i=10; # "p.H2c.H3b"
  out[ out$bg==bg,i+1] = round( pchisq( 2*(sub[i,6] - sub[8,6]), df=1,low=F), 3); ##??
  ps = (out[j,2:11] < 0.05)
  x = ps[ps==TRUE]
  out$sig[j] = length(x);
  #out$sig[j] = ifelse( min(out[j,2:8])<0.05, 1 ,0 )
}

out2 = out;
for( row in 1:length(out2[,1]) ){
  for ( col in 2:11 ) {
    out2[row,col] = ifelse( out[row,col]< 0.05, out[row,col], '');
  }
}

out2;

idtab = read.table("/home/hqin/coat.protein07/key.data/BGcoat.csv", header=T, sep="\t");
idtab$NCBI = as.character( idtab$NCBI );
idtab$Gene = as.character( idtab$Gene );

out2$name = idtab$Gene[ match( out2$bg, idtab$NCBI) ];

write.csv(out2, "chisq.Bha5sub4cer.090208.csv")


q("yes");

 #out$s1 = ifelse( out$p.H0to1 < 0.05, 1,0);
 #out$s2 = ifelse( out$p.H0to2 < 0.05, 1,0);

  #out[ out$bg==bg,3] = pchisq( (sub[3,6] - sub[1,6]), df=1,low=F);
  #out[ out$bg==bg,4] = pchisq( (sub[4,6] - sub[1,6]), df=1,low=F);
  #out[ out$bg==bg,5] = pchisq( (sub[5,6] - sub[1,6]), df=1,low=F);
  #out[ out$bg==bg,6] = pchisq( (sub[6,6] - sub[1,6]), df=1,low=F);
  #out[ out$bg==bg,6] = pchisq( (sub[6,6] - sub[1,6]), df=1,low=F);
  #out[ out$bg==bg,6] = pchisq( (sub[6,6] - sub[1,6]), df=1,low=F);


