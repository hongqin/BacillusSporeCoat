#073108 pick BG clusters for paml 
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; 

my $debug= 2 ;

my $inclusterfile = '_073108.test.csv';
#my $sourcedir = '/home/hqin/projects/coat.protein07/orthog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1/';
my $allfaa = "/home/hqin/projects/coat.protein07/genomes/_merged.all.062108.faa";
my $allfna = "/home/hqin/projects/coat.protein07/genomes/_merged.all.062108.fna";

open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); shift @lines;
my @clusterdirs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($bg, @res) = split( /\t/, $line);
  push (@clusterdirs, $bg);
}

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 

  system( "clean_fasta_headers.01.pl -if ids.faa -of prot.faa -pos 0 " );

  system( "clean_fasta_headers.01.pl -if ids.fna -of cds.fna  -pos 0 " );

  chdir "..";
}

exit;

#
# END
#

########################

