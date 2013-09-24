#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

my $debug = 0;

my $coatfl = "_simple.coat.orth.list.122307.versionD.csv";
if ($debug) { $coatfl = "_test.csv"; }

my $root = "/home/hqin/coat.protein07/orthog.analysis/phylo122307";

open (IN, "<$coatfl");
my $count = 1;
while (my $line = <IN>) {
  chomp $line; 
  my ( $bg, @ids) = split ( /\t/, $line );

  chdir $root;
  chdir $bg; 
  if ($debug) {    system( "pwd" );  }

  ### align CDS
  my $cmd1 =  "align_cds_by_protein_bioperl.00.pl -i cds.fasta -o tmpcds.phylip -af phylip ";
  system( $cmd1 );

  ### indicate it is interleaved.
  open(IN2, "<tmpcds.phylip"); my @lines = <IN2>; 
  close (IN2);

  open(OUT, ">cds.phylip");
    my $firstline = shift @lines;
    chomp $firstline;
    print OUT $firstline."  I\n";
    foreach my $line (@lines) { print OUT $line; }
  close(OUT);

  $count ++;
  print "**** count = $count \n\n";
  chdir ".."; 
}
close (IN);

exit;



