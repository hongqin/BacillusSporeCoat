#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

my $debug = 1;

my $allfaa = "_merged.all.122107.faa";
my $allfna = "_merged.all.122107.fna";
my $coatfl = "_simple.coat.orth.pairs.122007.csv";

my $root = "/home/hqin/coat.protein07/orthog.analysis/phylo/";

open (IN, "<$coatfl");
my $count = 0;
while (my $line = <IN>) {
  chomp $line; 
  my ( $bg, @ids) = split ( /\t/, $line );

  chdir $bg;
 
  if ($debug) {    system( "pwd" );  }

  my $cmd1 = "ln -sf $bg.faa prot.fasta";
  my $cmd2 = "ln -sf $bg.fna cds.fasta";
  system( $cmd1 );
  system( $cmd2 );

  

  $count ++;
  print "**** count = $count \n\n";

  chdir $root;
}
close (IN);

exit;


