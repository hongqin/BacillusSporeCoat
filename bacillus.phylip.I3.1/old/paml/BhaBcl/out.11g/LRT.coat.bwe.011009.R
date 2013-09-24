# focus on models related to bwe lineage
#LRT use 2 *(loglikelihood1 - loglikelihood2)

rm(list=ls());
 debug = 1;
 pickSig = 'NO'; qcutoff=0.05; #'YES' will use qcutoff 
 my.model = 'H1bwe';    my.df = 1;  #for H0-H2 pchisq test

 library(qvalue)

#### first, do lnL 
 tblnL = read.table( "out.coat.11glnL.csv", header=F );
 tblnL[,1] = as.character( tblnL[,1] );
 names(tblnL) = c("bg","model",NA,NA,NA,NA,"lnL",NA);

 # parse the lnL results
 tbH0     = tblnL[grep("H0",     tblnL[,2]), ]
 # tbHprime = tblnL[grep(my.model, tblnL[,2]), ] 
 
myLRT = function( my.df, cutoff, tag ) {
 # LRT using 2 *(loglikelihood1 - loglikelihood2)
 for( j in 1:length(tbHprime[,1]) ) {
  my.bg = tbHprime[j,1];
  lnLHA     = tbHA$lnL[tbHA$bg==my.bg];
  lnLHprime = tbHprime$lnL[tbHprime$bg==my.bg]; 
  deltalnL  = - lnLHA  + lnLHprime 
  tbHprime$deltalnL[j] = deltalnL ;
  tbHprime$p[j] = round( pchisq( 2* deltalnL, df=my.df,low=F), 4);
 }

 #qobj = qvalue( tbHprime$p);
 #tbHprime$q = qobj$qvalues; 
 #tbHprime$sig = ifelse( tbHprime$q < qcutoff, tag, '' );

 tbHprime$sig = ifelse( tbHprime$p < cutoff, tag, '' );
 tbHprime;
}

###
# my.models = c('H1C1','H2C2bwe','H3C2bweS1');    
# my.dfs    = c(1,      2,        3);  

 tbH1w     = tblnL[grep('H1w',   tblnL[,2]), ] 
 tbH2C2w   = tblnL[grep('H2C2w',   tblnL[,2]), ] 

 tbHprime = tbH1w;
 tbHA = tbH0;
 tbH1w = myLRT( 1, 0.05, '1b' );

 tbHprime = tbH2C2w;
 tbHA = tbH0;
 tbH2C2w = myLRT( 2, 0.05, '2b');

 table( tbH1w$sig );
 table( tbH2C2w$sig );

  out = tbH1w[-c(2:8)]
 
###second merge with omega
 tbomega = read.table("__omega.coat.11g.011009.csv",header=F,fill=T,sep="\t");
 for( i in c(1:3,5,7,9,11) ) {
  tbomega[,i] = as.character(tbomega[,i])
 }
 
 #my.model ='H1w';
 my.model ='H2C2w';
 omega = tbomega[grep(my.model, tbomega$V2), ]
 labels= c("bg","model",NA,"wbsu",NA,"wbwe",NA,"wbce",NA,"wbth","wban");
 names(omega) = labels;

 omega$sig    = out$sig[match(omega$bg, out$bg)]
 omega$sig2   = out$sig2[match(omega$bg, out$bg)]

q("yes");

#####################

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

 tbH1bwe$coatflag  = tbc$type[match(tbH1bwe$bg, tbc$id)]
 tbH1bwe$essenflag = tbe$type[match(tbH1bwe$bg, tbe[,1])]

 ## report for publications
 table(tbH1bwe$summary)


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
