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
  tr = read.nexus(\"prot-bacillus-paup-nj.tre\");
  nodes = tr\$tip.label
  Bha = grep( \"Bha\", nodes);
  if ( length(Bha) >= 1 ) {
    trBha1 = root( tr, nodes[Bha] );
    write.tree(trBha1, \"prot-bacillus-paup-nj-BhaRoot1.nwk\");
  }

  Bcl = grep( \"Bcl\", nodes);
  if ( length(Bcl) >= 1 ) {
    trBcl1 = root( tr, nodes[Bcl] );
    write.tree(trBcl1, \"prot-bacillus-paup-nj-BclRoot1.nwk\");
  }
  ";

  open (RRR, ">_tmp.R"); print RRR $Rlines; close (RRR);

  system( "R --no-save < _tmp.R" );

  chdir "..";
}

exit;

#
# END
#

