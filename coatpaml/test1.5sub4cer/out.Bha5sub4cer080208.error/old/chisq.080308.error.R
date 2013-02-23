
 tb = read.table( "out.Bha5sub4cer080208.csv", header=F);
 tb[,1] = as.character( tb[,1] );

 bgs = unique( tb[,1]);

 labels = c("bg","p.H0.H1a", "p.H0.H1b","p.H0.H1c","p.H0.H1d","p.H0.H2a","p.H0.H2b","p.H0.H2c" );
 out = data.frame( matrix(nrow=length(bgs), ncol=length(labels) )); 
 names( out ) = labels;
 out$bg = bgs;

 for( j in 1:length(bgs) ) {
  bg = bgs[j];
  sub = tb[ tb[,1] == bg,]
  for( i in 2:8 ) {
    out[ out$bg==bg,i] = round( pchisq( (sub[i,6] - sub[1,6]), df=1,low=F), 3);
  }
  ps = (out[j,2:8] < 0.05)
  x = ps[ps==TRUE]
  out$sig[j] = length(x);
  #out$sig[j] = ifelse( min(out[j,2:8])<0.05, 1 ,0 )
}

out2 = out;
for( row in 1:length(out2[,1]) ){
  for ( col in 2:8 ) {
    out2[row,col] = ifelse( out[row,col]< 0.1, out[row,col], '');
  }
}

out2;

write.csv(out, "chisq.Bha5sub4cer.080208.csv")

 #out$s1 = ifelse( out$p.H0to1 < 0.05, 1,0);
 #out$s2 = ifelse( out$p.H0to2 < 0.05, 1,0);

  #out[ out$bg==bg,3] = pchisq( (sub[3,6] - sub[1,6]), df=1,low=F);
  #out[ out$bg==bg,4] = pchisq( (sub[4,6] - sub[1,6]), df=1,low=F);
  #out[ out$bg==bg,5] = pchisq( (sub[5,6] - sub[1,6]), df=1,low=F);
  #out[ out$bg==bg,6] = pchisq( (sub[6,6] - sub[1,6]), df=1,low=F);
  #out[ out$bg==bg,6] = pchisq( (sub[6,6] - sub[1,6]), df=1,low=F);
  #out[ out$bg==bg,6] = pchisq( (sub[6,6] - sub[1,6]), df=1,low=F);


