#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

### align using muscle, for visual check.

my $debug = 1;

my $genomedir = '/home/hqin/projects/coat.protein07/genomes';
my $rootdir = '/home/hqin/coat.protein07/orthog.analysis/phylo012708';

### get the coat ids
system("ls -d BG* > /tmp/_ls-file.txt");
my %coats = ();
open (IN, "</tmp/_ls-file.txt");
while (my $line = <IN> ) {
   chomp $line; 
   my($bg, @els) = split (/-/, $line); 
   $coats{$bg} = $bg;
}
close(IN);
if ($debug) { showHashTable( \%coats); }

my $count = 0; 
foreach my $bg ( sort ( keys %coats ) ) {
  chdir $rootdir;
  chdir $bg; 
  if ($debug) {    system( "pwd" );  }

  ### align protein
  my $cmd1 =  "muscle -in prot.fasta -out prot-aligned.fasta -stable ";
  system( $cmd1 );
  my $cmd1b =  "muscle -in prot.fasta -out prot-aligned.aln -stable -clw ";
  system( $cmd1b );

  ### convert to nexus
  open  (OUT, ">seqCat.in");   print OUT "prot-aligned.fasta\n";   close (OUT);  
  my $cmd2 = "seqCat.pl -dseqCat.in -if";
  system( $cmd2 );
  system( "mv seqCat_sequences.nex prot-paup.nex "  );

  $count ++;
  print "**** count = $count $bg \n\n";
  if (($debug > 9)&&($count > 2)) {  exit; }
}


exit;


