use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 2 ;

my $clusterfl = "infile-BGclusters.bootstrap.txt";

open (IN, "<$clusterfl" ); my @clusterdirs = <IN>; close (IN);

my $trees="prot-bacillus-paup-nj-brlens.tre
prot-bacillus-paup-njbt.tre
";
#prot-bacillus-paup-nj.tre


my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  open (OUT,">__tree.txt"); print OUT $trees; close (OUT);

  my $Rlines = "
  library(ape);
  files = read.table( \"__tree.txt\")
  files = files[,1];

  pdf (\"_prot.nj.plot.070508.pdf\", width=10,height=8);
  par(mar=c(2,2,2,2)+0.1);

  for( i in 1:length(files) ) {
     file = as.character(files[i]);
     tr = read.nexus( file  );
     plot( tr, main=paste(\"$mydir\",file), show.node.label=T );
     plot( tr, main=paste(\"$mydir\",file), type=\"u\" );
     plot( tr, main=paste(\"$mydir\",file), type=\"c\", show.node.label=T  );
     i;
  }

  tnj = read.nexus( as.character(files[1]) );
  tbs = read.nexus( as.character(files[2]) );

  d = dist.topo( tnj, tbs);
  if( d[1]== 0 ) { 
    write.table( d, \"__nj.unchanged.in.bootstrap.070508.txt\");
  } else {
    write.table( d, \"__nj.changed.in.bootstrap.070508.txt\");
  } 

  dev.off();
  ";

  open (RRR, ">_tmp.plot.compare.R"); print RRR $Rlines; close (RRR);

  system( "R --no-save < _tmp.plot.compare.R > _out.tmp.R.txt" );

  chdir "..";
}

exit;

#
# END
#

