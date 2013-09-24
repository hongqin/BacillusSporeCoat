#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

#protein ML tree in paup.

my $debug = 0;

my $coatfl = "_simple.coat.orth.pairs.122007.csv";
if($debug) { $coatfl = "_test.csv"; }

my $root = "/home/hqin/coat.protein07/orthog.analysis/phylo/";

open (IN, "<$coatfl");
my @coatlines = <IN>; 
close (IN);
my $count = 0;
#while (my $line = <IN>) {

foreach my $line (@coatlines) {
  chomp $line; 
  my ( $bg, @ids) = split ( /\t/, $line );

  chdir $root;
  chdir $bg; 
  if ($debug) {    system( "pwd" );  }

  ### generate nj 
  open (IN, "<prot-paup.nex"); my @lines = <IN>; close (IN);
  pop @lines; ### remove the end;
  pop @lines; ### remove charset
  pop @lines; ### remove begin paups;

  open (OUT, ">prot-paup.nexus"); 
  foreach my $line (@lines) { 
    if ( $line =~ "datatype = DNA" ) { 
	$line =~ s/DNA/protein/o;
    } 
    print OUT $line; 
  }

  my $pauplines = "
BEGIN PAUP;
   NJ;
   savetrees file=prot-paup-nj.tre brlens;
END;
  ";
  print OUT "$pauplines";
  close (OUT);
 
  system( "paup prot-paup.nexus -n " );

  ###########boot strap
  open (OUT, ">prot-paup2.nexus"); 
  foreach my $line (@lines) { 
    if ( $line =~ "datatype = DNA" ) { 
	$line =~ s/DNA/protein/o;
    } 
    print OUT $line; 
  }

  my $pauplines2 = "
BEGIN PAUP;
   NJ;
   bootstrap;
   savetrees file=prot-paup-nj-btstrp.tre brlens;
END;
  ";
  print OUT "$pauplines2";
  close (OUT);
 
  system( "paup prot-paup2.nexus -n" );

  $count ++;   print "**** count = $count \n\n";

 # if( ( $debug ) & ( $count > 2 ) ) { exit; }  
}


exit;


