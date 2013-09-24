use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 2 ;

my $clusterfl = "infile-BGclusters.txt";

open (IN, "<$clusterfl" ); my @clusterdirs = <IN>; close (IN);

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  my $Rlines = "
  library(ape);
  tr1 = read.nexus(\"prot-bacillus-paup-nj.tre\");
  write.tree(tr1, \"prot-bacillus-paup-nj.nwk\");  ";
  open (RRR, ">_tmp.R"); print RRR $Rlines; close (RRR);

  system( "R --no-save < _tmp.R" );

  chdir "..";
}

exit;

#
# END
#

