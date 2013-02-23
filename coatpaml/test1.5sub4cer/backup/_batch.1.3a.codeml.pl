#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

my $debug= 1 ;

my $inclusterfile = '_073108.test.csv';
open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); shift @lines;
my @clusterdirs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($bg, @res) = split( /\t/, $line);
  push (@clusterdirs, $bg);
}

my $root = "/home/hqin/projects/coat.protein07/orthog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1/paml/test1.5sub4cer";

####################################################################

my $count = 0;
foreach my $mydir (@clusterdirs) {
  if (($debug>9) &&($count>1) ) { exit; }
  chdir $root;
  chdir $mydir; 

  print "\n******* count = $count $mydir \n"; #  system( "pwd" );


  my @Hs = ( "H0", "H1", "H2");
  # my @Hs = ( "H0");

  foreach my $H (@Hs) {
        chdir "$root/$mydir/$H";
        if ($debug) {    system( "pwd; ls" );  }
        system( "codeml");
  }

  $count ++;
}
close (IN);

exit;

