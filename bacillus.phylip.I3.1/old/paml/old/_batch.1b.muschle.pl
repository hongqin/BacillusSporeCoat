#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

### align using muscle

my $debug = 0;

my $allfaa = "_merged.all.122107.faa";
my $allfna = "_merged.all.122107.fna";
my $coatfl = "_simple.coat.orth.pairs.122007.csv";

if( $debug ) {    $coatfl = "_test.csv";   }

my $root = "/home/hqin/coat.protein07/orthog.analysis/phylo/";

open (IN, "<$coatfl");
my $count = 0;
while (my $line = <IN>) {
  chomp $line; 
  my ( $bg, @ids) = split ( /\t/, $line );

  chdir $bg; 
  if ($debug) { system( "pwd" );  }

  ### align protein
  my $cmd1 =  "muscle -in prot.fasta -out prot-aligned.fasta -stable ";
  system( $cmd1 );

  ### convert to nexus
  open  (OUT, ">seqCat.in");   print OUT "prot-aligned.fasta\n";   close (OUT);  
  my $cmd2 = "seqCat.pl -dseqCat.in -if";
  system( $cmd2 );
  system( "mv seqCat_sequences.nex prot-paup.nex "  );

  #system( "rm seqCat.in" );

  $count ++;
  print "**** count = $count \n\n";

  chdir $root;
}
close (IN);

exit;


