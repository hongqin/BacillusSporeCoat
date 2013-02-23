use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; 

my $debug= 2 ;
 
my $inprofilefile = '_10g.coat.profile.011009.csv';

system("pwd > __mydir.txt"); open(IN,"<__mydir.txt"); my @tmplines=<IN>; close(IN);
my $wkdir = shift @tmplines; chomp $wkdir; 

open (IN,"<$inprofilefile"); my @lines = <IN>; close (IN);
my $header = shift @lines; 

my $count = 0;
foreach my $line ( @lines ) { 
  chdir $wkdir;
  my @els = split( /\t/, $line );
  my $mydir = $els[0];
  $mydir =~ s/\"//g;
  chdir $mydir;

  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 

  system( "clean_fasta_headers.01.pl -if ids.faa -of prot.faa -pos 0 " );

  system( "clean_fasta_headers.01.pl -if ids.fna -of cds.fna  -pos 0 " );

}

exit;

#
# END
#

########################

