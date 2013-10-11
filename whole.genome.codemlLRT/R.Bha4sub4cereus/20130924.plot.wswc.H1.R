# plot ws ~ wc for 2013Sep revision

# use only genes with treedist<= 2 for publication

# LRT use 2 *(loglikelihood1 - loglikelihood2)

rm(list=ls());
 debug = 9;
 pickSig = 'NO'; 
 qcutoff=0.05; #'YES' will use qcutoff 
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

### First, filter the old results
 ## read the list of ids with treedist <= 2
 filter = read.table( "_BGids.treedist.lessthan2.050409.csv", header=F);
 names(filter) = c("id", "treedist");
 filter$id = as.character( filter$id );

 ## read old lnL 
 tblnL = read.table( "out.all.Bha4sub4cereus.lnL.010509.csv", header=F );
 tblnL[,1] = as.character( tblnL[,1] );
 names(tblnL) = c("bg","model",NA,NA,NA,"lnL",NA);
 tblnL.old = tblnL

 ## filter old lnL table
 #mt = match(filter$id, tblnL$bg)  #error, this gives only 1 record
 #tblnL.filtered = tblnL[ mt, ]

 tblnL$selected = 0;
 for( i in 1:length(filter$id) ) {
   currentID = filter$id[i]; 
   tblnL$selected[ grep(currentID, tblnL$bg)] = 1
 } #match does not work here.

 tblnL = tblnL[ tblnL$selected == 1, ]
 #tblnL = tblnL.old #2013Jan30 change

 x= intersect( filter$id, tblnL.old$bg );
 y= intersect( filter$id, tblnL$bg );
 identical(x,y) #[1] TRUE
 head(tblnL)
 head(tblnL.old)
 unique(tblnL$bg) #1174
 unique(tblnL.old$bg)

#### Now, do LRT
 # parse the lnL results
 tbH0     = tblnL[grep("H0",     tblnL[,2]), ]
 tbH1C1     = tblnL[grep('H1C1',   tblnL[,2]), ] 
 tbH1S1     = tblnL[grep('H1S1',   tblnL[,2]), ] 
 tbH2C1S1   = tblnL[grep('H2C1S1', tblnL[,2]), ] 
 
myLRT2 = function( tbHprime, tbHA, my.df, qcutoff, tag ) {
 # LRT using 2 *(loglikelihood1 - loglikelihood2)
 for( j in 1:length(tbHprime[,1]) ) {
  # print(j);
  my.bg = tbHprime[j,1];
  lnLHA     = tbHA$lnL[tbHA$bg==my.bg];
  lnLHprime = tbHprime$lnL[tbHprime$bg==my.bg]; 
  deltalnL  = - lnLHA  + lnLHprime 
  tbHprime$deltalnL[j] = deltalnL ;
  tbHprime$p[j] = round( pchisq( 2* deltalnL, df=my.df,low=F), 4);
 }

 qobj = qvalue( tbHprime$p );
 tbHprime$q = qobj$qvalues; 
 tbHprime$sig = ifelse( tbHprime$q < qcutoff, tag, '' );
 tbHprime;
 return(tbHprime)
}

 # tbHprime$sig = ifelse( tbHprime$q < qcutoff, tag, '' );
 tbH1C1   = myLRT2( tbH1C1,   tbH0, 1, 0.05, 'C' );
 tbH1S1   = myLRT2( tbH1S1,   tbH0, 1, 0.05, 'S' );
 tbH2C1S1 = myLRT2( tbH2C1S1, tbH0, 2, 0.05, 'CS' );

 head(tbH1C1)
 head(tbH1S1);
 table( tbH1C1$sig );
 table( tbH1S1$sig );
 #table( tbH2C1S1$sig );#not meaningful

 ## summarize the test resutls
 tbH1C1$sigS1 = tbH1S1$sig[match(tbH1C1$bg, tbH1S1$bg)] 
 for( j in 1:length(tbH1C1[,1]) ) {
   tbH1C1$summary[j] = paste(tbH1C1$sig[j], tbH1C1$sigS1[j], sep='');  
 }
 table( tbH1C1$summary );

####now, do H1-H2 test and update summary
 for( j in 1:length(tbH1C1[,1]) ) {
   my.bg = tbH1C1$bg[j];
   tag = tbH1C1$summary[j]; 
   if (tag == 'C' | tag=='CS' ) { #H1C1 vs H2C1S1
    lnLHA     = tbH1C1$lnL[tbH1C1$bg==my.bg];
    lnLHprime = tbH2C1S1$lnL[tbH2C1S1$bg==my.bg]; 
    deltalnL  = - lnLHA  + lnLHprime 
    tbH1C1$delCS2[j] = deltalnL[1] ;
    tbH1C1$pCS2[j]   = round( pchisq( 2* deltalnL, df=1,low=F), 4);
   } else if ( tag == 'S' ) { #H1S1 vs H2C1S1
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
 tbH1C1$sig2 = ifelse( tbH1C1$qCS2 < qcutoff  , 'CS2',  NA ); #CS2 is two step test and CS is one-step.
 tbH1C1[!is.na(tbH1C1$sig2), ][1:10,]

 table( tbH1C1$sig2 ); #86 H2C1S1 <0.05 fdr for publication
 
 out = tbH1C1[-c(2:4,6)]
 # outfile = paste( "_", my.model, "LRT.all.Bha4sub4cer.010609.csv", sep='.');
 #outfile = paste( "_nest.LRT.Bha4sub4cer.050409.csv", sep='.');
 #write.csv(out, outfile , quote=F, row.names=F)
 #I need put these summary report to manuscript and ppt. 
#### End, do H1-H2 test and update summary
 
### partition coat, noncoat, essential genes
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

 #omega$coatflag  = tbc$type[match(omega$bg, tbc$id)]
 #omega$essenflag = tbe$type[match(omega$bg, tbe[,1])]

 tbH1C1$coatflag  = tbc$type[match(tbH1C1$bg, tbc$id)]
 tbH1C1$essenflag = tbe$type[match(tbH1C1$bg, tbe[,1])]

##### report for publications, Venn diagram
 table(tbH1C1$summary) #1174, 86 in C, 310 in CS, 97 in S
 table(tbH1C1$summary, tbH1C1$coatflag) #8 in S, 5 in CS
 table(tbH1C1$summary, tbH1C1$essenflag) #6 in C, 32 in CS, 19 in S

 table(tbH1C1$coatflag) #19
 table(tbH1C1$essenflag) #182

 #H2
 table(tbH1C1$sig2) #86
 coat  = tbH1C1[ ! is.na(tbH1C1$coatflag), ]
 essen = tbH1C1[ ! is.na(tbH1C1$essenflag), ]

 table(coat$sig2) #3 coat genes H2
 table(essen$sig2) #16 essential genes H2
 
 table(coat$summary) #
 table(essen$summary) 


########2013 Sep 19, plot ws ~ wc for the two-step H2C1S1

 head(tbH1C1[ ! is.na(tbH1C1$sig2), ])
 table( tbH1C1$sig2 ); #86 H2C1S1 <0.05 fdr for publication

 tbomega = read.table("__omega.all.Bha4sub4cereus.010509.csv",header=F,fill=T,sep="\t");
 tbomega$coatflag  = tbc$type[match(tbomega$V1, tbc$id)]
 tbomega$CoatLocation  = tbc$CoatLocation[match(tbomega$V1, tbc$id)]
 tbomega$essenflag = tbe$type[match(tbomega$V1, tbe[,1])]
 tbomega$noCEflag = ifelse(is.na(tbomega$coatflag)&is.na(tbomega$essenflag), 'nonCE', NA)
 head(tbomega[tbomega$V2=="H2C1S1",  ])
 tbomegaH2C1S1 = tbomega[tbomega$V2=="H2C1S1",  ] #I have not applied tree-dist yet. 

 #apply treedis<=2 filter
 tbomegaH2C1S1 = tbomegaH2C1S1[match(filter$id, tbomegaH2C1S1$V1), ]
 tbomegaH2C1S1 = tbomegaH2C1S1[ ! is.na(tbomegaH2C1S1[,1]), ]

 #I need to use tbH1C1$sig2 to pick omegas from tbomega[tbomega$V2=="H2C1S1",  ]
 picked.BG.sig2 = as.vector( tbH1C1[ ! is.na(tbH1C1$sig2), "bg"])
 picked.omega = tbomegaH2C1S1[ match(picked.BG.sig2, tbomegaH2C1S1$V1), ]
 picked.omega$color = ifelse( picked.omega$CoatLocation=="OuterCoat", 'red', 
                              ifelse(picked.omega$coatflag=='coat','blue',NA) )
 #picked.omega$color = ifelse( picked.omega$essenflag=="essential", 'black', picked.omega$color)

 #plot wc ~wc for all H2C1S1 model
 tbomegaH2C1S1$CoatLocation = as.character(tbomegaH2C1S1$CoatLocation )
 tbomegaH2C1S1$color = NA;
 for( i in 1:length(tbomegaH2C1S1[,1])) {
# tbomegaH2C1S1$color[i] = ifelse( tbomegaH2C1S1$CoatLocation[i]=="OuterCoat", 'red', tbomegaH2C1S1$color[i] )
# tbomegaH2C1S1$color[i] = ifelse( tbomegaH2C1S1$CoatLocation[i]=="InnerCoat", 'blue', tbomegaH2C1S1$color[i] )
   tbomegaH2C1S1$color[i] = ifelse( tbomegaH2C1S1$coatflag[i]=="coat", 'red', 
            ifelse( tbomegaH2C1S1$essenflag[i]=="essential", 'black', NA )
   )
}
 table( tbomegaH2C1S1$color)
 table( tbomegaH2C1S1$coatflag)
 table( tbomegaH2C1S1$essenflag)
# pdf("ws_wc_20130924_H2C1S1all.pdf",width=5,heigh=5)
 plot( tbomegaH2C1S1$V6, tbomegaH2C1S1$V4, xlab='wc', ylab='ws', col=tbomegaH2C1S1$color ) 
x = seq(0,0.3,by=0.01); y = x; 
lines( y ~ x, lty=2 )
dev.off()

 pdf("ws_wc_20130919_sig2.pdf",width=5,heigh=5)
 with( picked.omega[!is.na(picked.omega$noCEflag),], 
       plot( V6,V4, ylab=expression(paste(omega,s)), xlab=expression(paste(omega,c)),                                           
             xlim=c(0,0.25),ylim=c(0,0.25)   )       
       )
 #with( picked.omega, plot( V4 ~ V6, ylab='ws', xlab='wc', col=picked.omega$color))
 #plot( picked.omega$V4 ~ picked.omega$V6, col=picked.omega$color)
 #with( picked.omega[!is.na(picked.omega$coatflag), ], points( V6, V4, col='red', pch=19))
 with( picked.omega[!is.na(picked.omega$coatflag), ], points( V6, V4, col=color, pch=19))
 with( picked.omega[!is.na(picked.omega$essenflag), ], points( V6, V4, col='black', pch=17))
 x = seq(0,0.3,by=0.01); y = x; 
 lines( y ~ x, lty=2 )
 dev.off()

 #reivewer 2 comment, ws > wc? 2013 Sep 23, using significant ones
 ws = picked.omega$V4[picked.omega$coatflag=='coat']
 wc = picked.omega$V6[picked.omega$coatflag=='coat']
 summary(ws); summary(wc)
 wilcox.test( ws, wc, al='gr') #p=0.65

 ws = picked.omega$V4[picked.omega$essenflag=='essential']
 wc = picked.omega$V6[picked.omega$essenflag=='essential']
 summary(ws); summary(wc)
 wilcox.test( ws, wc, al='gr') #p=0.65
 
#reivewer 2 comment, ws > wc? 2013 Sep 24, using all H2C1S1 coat genes
 wsC = tbomegaH2C1S1$V4[tbomegaH2C1S1$coatflag == 'coat']
 wcC = tbomegaH2C1S1$V6[tbomegaH2C1S1$coatflag == 'coat']
 summary(wsC); summary(wcC)
 ratioWsWcCoat = ws/wc
 summary(ratioWsWcCoat)
 wilcox.test( ratioWsWcCoat, mu=1, al="gr") #p=0.072
 t.test( ratioWsWcCoat, mu=1, al="gr") #p=0.099

plot( ws ~ wc, col='red', xlim=c(0, 0.5), ylim=c(0,0.5))

wsE = tbomegaH2C1S1$V4[tbomegaH2C1S1$essenflag == 'essential']
wcE = tbomegaH2C1S1$V6[tbomegaH2C1S1$essenflag == 'essential']
summary(wsE); summary(wcE)
ratioWsWcEssen = wsE/wcE
summary(ratioWsWcEssen)
wilcox.test( ratioWsWcEssen, mu=1, al="gr") #p=0.081
t.test( ratioWsWcEssen, mu=1, al="gr") #p=0.015
plot( wc, ws, col='blue')

pdf("ws_wc_20130924_H2C1S1all.pdf",width=5,heigh=5)
plot( tbomegaH2C1S1$V6, tbomegaH2C1S1$V4, xlab='wc',ylab='ws', pch=3, col='gray')
points( wcE, wsE, col='blue', pch=2)
x = seq(0,0.3,by=0.01); y = x; 
lines( y ~ x, lty=2 )
points( wcC, wsC, col='red', pch=19)
Mylabels = c("coat", "essential", "nonCE")
Mysymbols = c(19, 2, 3)
Mycolors = c("red", "blue", "gray")
legend(0.0, 0.26, legend=Mylabels, col=Mycolors, pch=Mysymbols)
dev.off()

quit("yes");

