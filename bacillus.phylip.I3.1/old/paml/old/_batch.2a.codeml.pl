#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

my $debug = 0;

# run codeml in H0, H1 H2

my $root = "/home/hqin/coat.protein07/orthog.analysis/phylo122307";

my $coatfl = "$root/_simple.coat.orth.list.122307.versionD.csv";
if ($debug) { $coatfl = "_test.csv"; }

open (IN, "<$coatfl");
my $count = 0;
while (my $line = <IN>) {
  chomp $line; 
  my ( $bg, @ids) = split ( /\t/, $line );

  #chdir $root; chdir $bg; 
 
  my @Hs = ( "H0", "H1", "H2");
  # my @Hs = ( "H0");

  foreach my $H (@Hs) {
  	chdir "$root/$bg/$H";
  	if ($debug) {    system( "pwd; ls" );  }
  	system( "codeml");
  }

  ### 
  $count ++;
  print "**** count = $count \n\n";
}
close (IN);

exit;


