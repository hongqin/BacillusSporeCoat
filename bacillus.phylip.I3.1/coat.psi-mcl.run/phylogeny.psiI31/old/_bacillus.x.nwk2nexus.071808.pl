#061708 add plots
#071608 compare nj and bootstrap trees using treedist and consense

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
open( REPORT, ">__bacilus.4a.tree.comparison.report.071708.txt");
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-----------------------------------\nCluster:$count I am now working on [$mydir]::\n"; 
 
  system( "cat protbaci_consense.nwk protbaci_nj.nwk > intree" );
  open(OUT,">promt.1.4a.txt"); print OUT "D\nY\n"; close (OUT);
  if( -e "outfile" ) { system( "rm outfile"); }
  system( "treedist < promt.1.4a.txt" );
  
  open(IN, "<outfile"); my @lines = <IN>; close(IN);
  my $result = pop @lines;
  chomp $result;  
  my @els = split( /:/, $result );
  $result = $els[$#els];
  $result =~ s/\s+//g;
  print REPORT "$mydir treedist result::$result\n";
  system("mv outfile _outfile.treedist");

  if ( $result == 0 ) {
     system( "cp protbaci_nj.nwk $mydir-bacillus.nwk" );
     system("touch __same.topology.btstrap.n.nj.071608.txt");
  } else { #find consensus tree between bootstrap and nj
	open(OUT,">prompt.2.4a.txt"); print OUT "C\nY\n"; close (OUT);	
	if (-e "outtree" ) { system("rm outtree"); } 
	system("consense < prompt.2.4a.txt"); 
	if ( -e "outtree" ) {
	  system ("cat outtree");
	  system ("mv outtree $mydir-bacillus.nwk");
	  system("touch __consensus.of.btstrap.n.nj.071608.txt");
	} else { 
	  system("touch __ERROR.prot-bacillus-nwk.071508.txt");
	}
  }
  
  #plot the trees
  my $Rlines = "library(ape); 
  postscript(\"_plot.trees.ps\", horizontal=F, width=8,height=8);
  t1=read.tree(\"protbaci_nj.nwk\");
  plot(t1,main=\"$mydir nj treedist=$result\");
  t2=read.tree(\"protbaci_consense.nwk\");
  plot(t2,main=\"$mydir bootstrap treedist=$result\");
  t3=read.tree(\"$mydir-bacillus.nwk\");
  plot(t3,main=\"$mydir final treedist=$result\")
  dev.off();
  ";
  open(RRR,">__tmp.plot.R"); print RRR $Rlines; close(OUT);
  
  
  if ( $result > 0 )   {
    system( "R --no-save < __tmp.plot.R" );
  } 
  #system( "ps2pdf _plot.tree.ps" ); #bug ??

  chdir "..";
}
close(REPORT);

exit;

#
# END
#

