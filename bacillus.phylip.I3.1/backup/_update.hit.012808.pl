#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

#summarize for mysql 

my $debug = 5;

#my $genomedir = '/home/hqin/projects/coat.protein07/genomes';
my $rootdir = '/home/hqin/coat.protein07/orthog.analysis/phylo012708/mysql';
#my $clusterdir = '/home/hqin/projects/coat.protein07/orthog.analysis/all.vs.all/simple.hits.mcl/single.cluster';

my $infile = "_rbHits.coat.020108.csv";
my $outfile = "_rbHits.coat.020108.v2.csv";

my %spec_freq = ();

chdir $rootdir;
open(IN, "<$rootdir/$infile"); my @lines = <IN>; close(IN);

#count number in each species per bgcoat
foreach my $line ( @lines ) {
  chomp $line; 
  my ($bgname, $bgid, $hit, $spec, @res ) = split (/\t/, $line);
  $spec_freq{"$bgid\t$spec"} += 1;
}


open (OUT, ">$rootdir/$outfile");
foreach my $line ( @lines ) {
	chomp $line;
  	my ($bgname, $bgid, $hit, $spec, @res ) = split (/\t/, $line);
  	#print OUT "$bgname\t$bgid\t$hit\t". $spec_freq{"$bgid\t$spec"}. "\t$spec\n" ;
	print OUT "$bgname\t$bgid\t$hit\t$spec\t". $spec_freq{"$bgid\t$spec"}. "\n" ;
}
close (OUT);

exit;
