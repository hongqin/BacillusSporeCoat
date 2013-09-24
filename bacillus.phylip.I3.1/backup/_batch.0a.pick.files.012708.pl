#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

my $debug = 0;

#my $allfaa = "_merged.all.122107.faa";
#my $allfna = "_merged.all.122107.fna";

#my $coatfl = "_simple.coat.orth.pairs.122207.versionC.csv";
#if( $debug ) { $coatfl = "_test.csv"; }

my $genomedir = '/home/hqin/projects/coat.protein07/genomes';
my $rootdir = '/home/hqin/coat.protein07/orthog.analysis/phylo012708';
my $clusterdir = '/home/hqin/projects/coat.protein07/orthog.analysis/all.vs.all/simple.hits.mcl/single.cluster';

system("cat $genomedir/*/*faa > /tmp/_merged.faa");
system("cat $genomedir/*/*fna > /tmp/_merged.fna");

### get the cluster files
system("ls $clusterdir > /tmp/_ls-file.txt");
my %clusterfiles = ();
open (IN, "</tmp/_ls-file.txt");
while (my $line = <IN> ) {
   chomp $line; 
   my($bg, @els) = split (/-/, $line); 
   $clusterfiles{$bg} = "$clusterdir/$line";
}
close(IN);
if ($debug) { showHashTable( \%clusterfiles); }

my $count = 0; 
foreach my $bg ( sort ( keys %clusterfiles ) ) {
  chdir $rootdir;
  mkdir $bg; 
  chdir $bg; 
  if ($debug) {    system( "pwd" );  }

  system("ln -s $clusterfiles{$bg} $bg-single-cluster.tab");

  my $outfna = "$bg.original.fna";
  my $outfaa = "$bg.original.faa"; 

  my $cmd  = "pick_fasta_records_by_ids.01.pl -if /tmp/_merged.fna -id $bg-single-cluster.tab -o $outfna";
  print $cmd."\n\n";   system( $cmd );
  system( "clean_fasta_headers.00.pl -if $outfna -of /tmp/_tmp.fna -pos 0-1 " );
  system( "remove_1stSpace_in_fasta_header.00.pl -if /tmp/_tmp.fna -of /tmp/_tmp2.fna " );
  system( "clean_fasta_headers.01.pl -if /tmp/_tmp2.fna -of $bg.fna -pos 0 " );
  system( "ln -sf $bg.fna cds.fasta" );

  $cmd  = "pick_fasta_records_by_ids.01.pl -if /tmp/_merged.faa -id $bg-single-cluster.tab -o $outfaa";
  print $cmd."\n\n";  system( $cmd );
  system( "clean_fasta_headers.00.pl -if $outfaa -of /tmp/_tmp.faa -pos 0-1 " );
  system( "remove_1stSpace_in_fasta_header.00.pl -if /tmp/_tmp.faa -of /tmp/_tmp2.faa " );
  system( "clean_fasta_headers.01.pl -if /tmp/_tmp2.faa -of $bg.faa -pos 0 " );
  system( "ln -sf $bg.faa prot.fasta" );

  $count ++;
  print "**** count = $count $bg \n\n";
  if (($debug > 4)&&($count > 2)) {  exit; }
}


system("rm /tmp/_merged.faa");
system("rm /tmp/_merged.fna");

exit;

