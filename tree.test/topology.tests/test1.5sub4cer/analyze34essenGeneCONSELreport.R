
rm(list=ls())

setwd("~/projects/coat.protein07/ortholog.analysis/essetial.mcl.rbh.090208/bacillus.phylip.I3.1/topology.tests/test1.5sub4cer")
list.files()
tb = read.table( "_conselreport.34essenGene.2012Jan4.tab", sep="\t", header=T)
summary(tb)

sub.top = tb[tb$rank==1, ]
table( sub.top$item )

unique(tb$BG)

sub.accepted = tb[tb$au>0.05, ]
table( sub.accepted$item )
