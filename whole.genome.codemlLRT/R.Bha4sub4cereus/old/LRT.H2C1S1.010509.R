#LRT use 2 *(loglikelihood1 - loglikelihood2)

rm(list=ls());
 debug = 9;
 pickSig = 'NO'; qcutoff=0.05; #'YES' will use qcutoff 
 my.model = 'H2C1S1';    my.df = 2;  #for H0-H2 pchisq test
 my.models = c('H1C1','H1S1','H2C1S1');    
 my.dfs    = c(1,      1,        2);  

 library(qvalue)

# my @Hs = ( 'H0', 'H1C1','H1S1','H2C1S1','H2C2bwe','H3C2bweS1','H3C2S1');

#my @Htres = ( 
#"(Bha,(Bpu,(Bli,(Bam,Bsu))),(Bwe,(Bce,(Ban,Bth))));", #H0
#"(Bha,(Bpu,(Bli,(Bam,Bsu))),(Bwe #1,(Bce #1,(Ban #1,Bth #1)#1)#1)#1);", #H1C1
#"(Bha,(Bpu #1,(Bli #1,(Bam #1,Bsu #1)#1)#1)#1,(Bwe,(Bce,(Ban,Bth))));", #H1S1
#"(Bha,(Bpu #2,(Bli #2,(Bam #2,Bsu #2)#2)#2)#2,(Bwe #1,(Bce #1,(Ban #1,Bth #1)#1)#1)#1);", #H2C1S1
#"(Bha,(Bpu ,(Bli ,(Bam ,Bsu ))),(Bwe #2,(Bce #1,(Ban #1,Bth #1)#1)#1)#1);", #H2C2bwe
#"(Bha,(Bpu #2,(Bli #2,(Bam #2,Bsu #2)#2)#2)#2,(Bwe #1,(Bce #1,(Ban #3,Bth #1)#1)#1)#1);", #H3C2S1
#"(Bha,(Bpu #2,(Bli #2,(Bam #2,Bsu #2)#2)#2)#2,(Bwe #3,(Bce #1,(Ban #1,Bth #1)#1)#1)#1);", #H3C2bweS1
#);

#### first, do lnL 
 tblnL = read.table( "out.all.Bha4sub4cereus.lnL.010509.csv", header=F );
 tblnL[,1] = as.character( tblnL[,1] );
 names(tblnL) = c("bg","model",NA,NA,NA,"lnL",NA);

 # parse the lnL results
 tbH0     = tblnL[grep("H0",     tblnL[,2]), ]
 # tbHprime = tblnL[grep(my.model, tblnL[,2]), ] 
 
myLRT = function( my.df, qcutoff, tag ) {
 # LRT using 2 *(loglikelihood1 - loglikelihood2)
 for( j in 1:length(tbHprime[,1]) ) {
  my.bg = tbHprime[j,1];
  lnLHA     = tbHA$lnL[tbHA$bg==my.bg];
  lnLHprime = tbHprime$lnL[tbHprime$bg==my.bg]; 
  deltalnL  = - lnLHA  + lnLHprime 
  tbHprime$deltalnL[j] = deltalnL ;
  tbHprime$p[j] = round( pchisq( 2* deltalnL, df=my.df,low=F), 4);
 }

 qobj = qvalue( tbHprime$p);
 tbHprime$q = qobj$qvalues; 
 tbHprime$sig = ifelse( tbHprime$q < qcutoff, tag, '' );
 tbHprime;
}

###
 tbH1C1     = tblnL[grep('H1C1',   tblnL[,2]), ] 
 tbH1S1     = tblnL[grep('H1S1',   tblnL[,2]), ] 
 tbH2C1S1   = tblnL[grep('H2C1S1', tblnL[,2]), ] 

 tbHprime = tbH1C1;
 tbHA = tbH0;
 tbH1C1 = myLRT( 1, 0.05, 'C' );

 tbHprime = tbH1S1;
 tbHA = tbH0;
 tbH1S1 = myLRT( 1, 0.05, 'S');

 tbHprime = tbH2C1S1;
 tbHA = tbH0;
 tbH2C1S1 = myLRT( 2, 0.05, '2');

 head(tbH1C1)
 head(tbH1S1);
 table( tbH1C1$sig );
 table( tbH1S1$sig );
 table( tbH2C1S1$sig );

 ## how many are shared significant pairs?
 tbH1C1$sigS1 = tbH1S1$sig[match(tbH1C1$bg, tbH1S1$bg)] 
 tbH1C1$sig2  = tbH2C1S1$sig[match(tbH1C1$bg, tbH2C1S1$bg)] 
 for( j in 1:length(tbH1C1[,1]) ) {
   tbH1C1$summary[j] = paste(tbH1C1$sig[j], tbH1C1$sigS1[j], tbH1C1$sig2[j], sep='');  
 }
 table( tbH1C1$summary );


 ##now, do H1-H2 test and update summary
 for( j in 1:length(tbH1C1[,1]) ) {
   my.bg = tbH1C1$bg[j];
   tag = tbH1C1$summary[j]; 
   if (tag == 'C2' | tag=='CS2' ) { #H1C1 vs H2C1S1
    lnLHA     = tbH1C1$lnL[tbH1C1$bg==my.bg];
    lnLHprime = tbH2C1S1$lnL[tbH2C1S1$bg==my.bg]; 
    deltalnL  = - lnLHA  + lnLHprime 
    tbH1C1$delCS2[j] = deltalnL[1] ;
    tbH1C1$pCS2[j]   = round( pchisq( 2* deltalnL, df=1,low=F), 4);
   } else if ( tag == 'S2' ) { #H1S1 vs H2C1S1
    lnLHA     = tbH1S1$lnL[tbH1S1$bg==my.bg];
    lnLHprime = tbH2C1S1$lnL[tbH2C1S1$bg==my.bg]; 
    deltalnL  = - lnLHA  + lnLHprime 
    tbH1C1$delCS2[j] = deltalnL ;
    tbH1C1$pCS2[j]   = round( pchisq( 2* deltalnL, df=1,low=F), 4);
   } else {
    tbH1C1$delCS2[j] = NA;
    tbH1C1$pCS2[j]   = NA;
   } 
 }

 x = tbH1C1[ ! is.na(tbH1C1$pCS2), ] 
 qobj = qvalue( x$pCS2);
 x$qCS2 = qobj$qvalues; 
 tbH1C1$qCS2 = x$qCS2[match(tbH1C1$bg, x$bg)]
 summary( x$qC2 );
 summary( tbH1C1$qCS2 );
 tbH1C1[!is.na(tbH1C1$qCS2),][1:10,]
 
 qcutoff = 0.05; 
 tbH1C1$sig2 = ifelse( tbH1C1$qCS2 < qcutoff  , 'cs2',  tbH1C1$sig2 );
 tbH1C1[!is.na(tbH1C1$sig2), ][1:10,]

 table( tbH1C1$sig2 );
 
 out = tbH1C1[-c(2:4,6)]
 # outfile = paste( "_", my.model, "LRT.all.Bha4sub4cer.010609.csv", sep='.');
 outfile = paste( "_nest.LRT.Bha4sub4cer.010609.csv", sep='.');
 write.csv(out, outfile , quote=F, row.names=F)
 #I need put these summary report to manuscript and ppt. 
 
###second merge with omega
 tbomega = read.table("__omega.all.Bha4sub4cereus.010509.csv",header=F,fill=T,sep="\t");
 for( i in c(1:3,5,7) ) {
  tbomega[,i] = as.character(tbomega[,i])
 }
 
 my.model ='H2C1S1';
 omega = tbomega[grep(my.model, tbomega$V2), ]
 labels= c("bg","model",NA,"wbsu",NA,"wbwe",NA,"wban");
 names(omega) = labels;

 omega$summary   = out$summary[match(omega$bg, out$bg)]
 omega$sig2   = out$sig2[match(omega$bg, out$bg)]

 omega.original = omega; #keep a copy 

 pickSig='NO';
 ##pick only significant ones
 #pickSig = 'NO'; qcutoff=0.1; #'YES' will use qcutoff 
 if (pickSig == 'YES' ) {
  omega=omega[ omega$sig2=='cs2', ]
 }

 head(omega) #visual check the first 2 genes, passed.

###Thirdaly, partition coat, noncoat, essential genes
 ###essential bg
 tbe = read.table( "BG.essential.Kobayashi03.tab", header=F  );
 tbe[,1] = as.character( tbe[,1] );
 bge = tbe[,1];
 tbe$type = "essential";
 rownames(tbe) = tbe[,1];

 ###coat bg
 tbc = read.csv( "_coat.profile.122908.csv", sep="\t");
 tbc$id = as.character(tbc$id);
 tbc$type = "coat";
 rownames(tbc) = tbc$id;

 omega$coatflag  = tbc$type[match(omega$bg, tbc$id)]
 omega$essenflag = tbe$type[match(omega$bg, tbe[,1])]

 tbH1C1$coatflag  = tbc$type[match(tbH1C1$bg, tbc$id)]
 tbH1C1$essenflag = tbe$type[match(tbH1C1$bg, tbe[,1])]

 ## report for publications
 table(tbH1C1$summary)


####test by parition
 wcoat    = omega[!is.na(omega$coatflag),  ] 
 wessen   = omega[!is.na(omega$essenflag), ] 
 wnoncoat = omega[ is.na(omega$coatflag),  ] 

 #wbsu 
 summary(wcoat$wbsu); 
 summary(wessen$wbsu); 
 summary(wnoncoat$wbsu);
 dn = density( wnoncoat$wbsu );
 dc = density( wcoat$wbsu );
 de = density( wessen$wbsu );
 myY = max( dc$y, de$y, dn$y );
 plot( dn, xlim=c(0,0.5), ylim=c(0,myY) );
 lines( dc, col="red" );
 lines( de, col="blue");
 wilcox.test( wcoat$wbsu, wnoncoat$wbsu, alt="gr");


 #wbwe 
 summary(wcoat$wbwe); 
 summary(wessen$wbwe); 
 summary(wnoncoat$wbwe);
 dn = density( wnoncoat$wbwe );
 dc = density( wcoat$wbwe );
 de = density( wessen$wbwe );
 myY = max( dc$y, de$y, dn$y );
 plot( dn, xlim=c(0,0.5), ylim=c(0,myY) );
 lines( dc, col="red" );
 lines( de, col="blue");
 wilcox.test( wcoat$wbwe, wnoncoat$wbwe, alt="gr");


 #wban 
 summary(wcoat$wban); 
 summary(wessen$wban); 
 summary(wnoncoat$wban);
 dn = density( wnoncoat$wban );
 dc = density( wcoat$wban );
 de = density( wessen$wban );
 myY = max( dc$y, de$y, dn$y );
 plot( dn, xlim=c(0,0.5), ylim=c(0,myY) );
 lines( dc, col="red" );
 lines( de, col="blue");
 wilcox.test( wcoat$wban, wnoncoat$wban, alt="gr");




quit("yes");

#
 
###second merge with omega
 tbomega = read.table("__omega.all.Bha4sub4cereus.010509.csv",header=F,fill=T,sep="\t");
 for( i in c(1:3,5,7) ) {
  tbomega[,i] = as.character(tbomega[,i])
 }
 
 omega = tbomega[grep(my.model, tbomega$V2), ]
 labels= c("bg","model",NA,"wbsu",NA,"wbwe",NA,"wban");
 names(omega) = labels;

 omega$lnL = tbHprime$deltalnL[match(omega$bg, tbHprime$bg)]
 omega$p   = tbHprime$p[match(omega$bg, tbHprime$bg)]
 omega$q   = tbHprime$q[match(omega$bg, tbHprime$bg)]

 omega.original = omega; #keep a copy 


q("yes");
### END
###############################
