# 062208 align protein sequence using MUSCLE
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; 

my $debug= 2 ;

my $tmpfl = "/tmp/_coatBGclus.062208.txt";
 #   system( "ls -d BG* | head -n 2 > $tmpfl " );
system( "ls -d BG*  > $tmpfl " );
open (IN, "<$tmpfl" ); my @clusterdirs = <IN>; close (IN);
system( "rm $tmpfl" );

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  my $cmd1  =  "muscle -in prot.faa -out prot-aligned.fasta -stable ";
  system( $cmd1 );
  my $cmd1b =  "muscle -in prot.faa -out prot-aligned.aln   -stable -clw ";
  system( $cmd1b );

  chdir "..";
}

exit;

#
# END
#

########################

