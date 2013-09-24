#071808 add species labels to trees; plot for visual check

use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 0;

#my $idfl   = "/home/hqin/projects/coat.protein07/key.data/_id2genomes.csv";
#my $specfl = "/home/hqin/projects/coat.protein07/key.data/_genome2species.csv";
my $clusterfl = "infile-BGclusters.bootstrap.txt";

#my $tre ="((((((Bmo,Bsu),Bam),Bli),Bpu),((Bce,Ban,Bth),Bwe)),(Bha,Bcl))";
# OR in a more readble format: 
#	( 
#	(
#	 	( ( ( (Bmo,Bsu),Bam), Bli), Bpu ), 
#		( (Bce,Ban, Bth),Bwe)
#	), (Bha, Bcl) 
#	);

open (IN, "<$clusterfl" ); my @clusterdirs = <IN>; close (IN);

if($debug>9) { @clusterdirs = ( $clusterdirs[0], $clusterdirs[1] ); 
};

my $count = 0;
open( REPORT, ">__bacilus.4x.tree.comparison.report.071808.txt");
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n---------------------------------\nCluster:$count I am now working on [$mydir]::\n"; 
 
  system( "R --no-save < ../R.scripts/10char.2.longIDs.071808.R" );

  my $Rlines = "library(ape);
   postscript(\"_tree.071808.ps\",width=8,height=8,horizontal=F);
   t=read.tree(\"protbaci.label.nwk\");
   plot(t, main=\"$mydir final\");
   dev.off();
  ";
  open( RRR, ">__tmp.plot.R"); print RRR $Rlines; close(RRR);

  system("R --no-save < __tmp.plot.R");

  chdir "..";
}
close(REPORT);

exit;

#
# END
#

