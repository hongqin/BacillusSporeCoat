#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

my $debug = 0;

my $allfaa = "_merged.all.122107.faa";
my $allfna = "_merged.all.122107.fna";

my $coatfl = "_simple.coat.orth.list.122307.versionD.csv";
if( $debug ) { $coatfl = "_test.csv"; }

my $root = "/home/hqin/coat.protein07/orthog.analysis/phylo122307/";

open (IN, "<$coatfl");
my $count = 0;
while (my $line = <IN>) {
  chomp $line; 
  my ( $bg, @ids) = split ( /\t/, $line );

  #rmdir $bg;
  mkdir $bg;
  chdir $bg;
 
  if ($debug) {    system( "pwd" );  }

  open (OUT, ">$bg.ortho.txt");
  print OUT $line;
  close (OUT);

  my $outfaa = "$bg.original.faa"; 
  my $outfna = "$bg.original.fna";

  my $cmd  = "pick_fasta_records_by_ids.01.pl -if $root$allfna -id $bg.ortho.txt -o $outfna";
  print $cmd."\n\n";
  system( $cmd );
  system( "clean_fasta_headers.00.pl -if $outfna -of _tmp.fna -pos 0-1 " );
  system( "remove_1stSpace_in_fasta_header.00.pl -if _tmp.fna -of _tmp2.fna " );
  system( "clean_fasta_headers.01.pl -if _tmp2.fna -of $bg.fna -pos 0 " );
  system( "ln -sf $bg.fna cds.fasta" );

  my $cmd2 = "pick_fasta_records_by_ids.01.pl -if $root$allfaa -id $bg.ortho.txt -o $outfaa";
  print $cmd2."\n\n";
  system( $cmd2 );
  system( "clean_fasta_headers.00.pl -if $outfaa -of _tmp.faa -pos 0-1 " );
  system( "remove_1stSpace_in_fasta_header.00.pl -if _tmp.faa -of _tmp2.faa " );
  system( "clean_fasta_headers.01.pl -if _tmp2.faa -of $bg.faa -pos 0 " );
  system( "ln -sf $bg.faa prot.fasta" );

  $count ++;
  print "**** count = $count \n\n";

  chdir $root;
}

close (IN);

exit;


