
pretb = read.csv("_coat.profile.comparisons.Feb11,2010.csv", sep="\t");
psitb = read.csv("_profile.from.psi.mcl.I31.clus.021510.csv", sep="\t");

oldprofile = pretb[,c(1,14:24)]
psiprofile = psitb[,-c(15,16)]

for( c in 1:length(oldprofile[1,]) ) {
 oldprofile[,c] = as.character( oldprofile[,c]); 
}

for( c in 1:length(psiprofile[1,]) ) {
 psiprofile[,c] = as.character( psiprofile[,c]); 
}

#check the orders
x = match(psiprofile$coat, oldprofile$id );
x - 1:73 #good

comp = psiprofile; #store comparison results
oldnames = names(oldprofile)[2:12]
oldnames = oldnames[c(1:6, 7,7,7:11)]
psinames = names( psiprofile ) [2:14]

cbind( oldnames, psinames) #double check
oldnames == psinames

#compare profiles
for( r in 1 : length(oldprofile[,1]) ) {
  for( c in 1 : length(oldnames) ) {
    answer = '';
    oldhits = oldprofile[ r, oldnames[c] ];
    psihits = psiprofile[ r, psinames[c] ];
    if ( (is.na(oldhits)) & (is.na(psihits) ) ) { 
       answer = TRUE; 
    } else if (is.na(oldhits) ) {
       answer = psihits; 
    } else if ( is.na(psihits) ) {
       answer = oldhits;
    }  else {
       els1 = strsplit( oldhits, ' ' )[[1]];
       els2 = strsplit( psihits, ' ')[[1]];
       #answer = (oldhits == psihits )[1]
       inter = intersect( els1, els2 );
       if (is.na(inter[1]) ) { answer =FALSE;
       } else { answer = TRUE; }
    }  
    comp[r,(c+1)] = answer; 
  }
}

table( comp$Bsu );
table( comp$Bmo );
table( comp$Bam );
table( comp$Bli );
table( comp$Bpu );

table( comp$Ban );
table( comp$Bth );
table( comp$Bwe );
table( comp$Bcl );
table( comp$Bha );

table( comp$Bce87 );
table( comp$Bce79 );
table( comp$Bce3L );

#summarize Ban Bth Bce79 Bwe
for( r in 1 : length(oldprofile[,1]) ) {
    subcer4 = comp[r, c("Ban", "Bth", "Bce79", "Bwe")]; 
    tmp = table( t( subcer4 ) );
    comp$cer4[r] = tmp["TRUE"];
    if( is.na( tmp["TRUE"] ) ) { 
      comp$cer4[r] = 0;
    }
}

#Find inconsistent ones
inconsistent.coat = comp$coat[ comp$cer4 < 4 ]

inconsistent.oldprofile = oldprofile[ match(inconsistent.coat, oldprofile$id),  ]
inconsistent.psiprofile = psiprofile[ match(inconsistent.coat, psiprofile$coat), ]

inconsistent.oldprofile[ , c("id","Ban","Bth","Bce","Bwe")]
inconsistent.psiprofile[ , c("coat","Ban","Bth","Bce79","Bwe")] 

write.csv( inconsistent.oldprofile, "_inconsistent.oldprofile.022410.csv", row.names=F );
write.csv( inconsistent.psiprofile, "_inconsistent.psiprofile.022410.csv", row.names=F );

