# 062208 compare prot nj tree and 16s rRNA tree using NOTUNG
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 2 ;

#my $idfl   = "/home/hqin/projects/coat.protein07/key.data/_id2genomes.csv";
#my $specfl = "/home/hqin/projects/coat.protein07/key.data/_genome2species.csv";
my $clusterfl = "infile-BGclusters.txt";

my $tre ="((((((Bmo,Bsu),Bam),Bli),Bpu),((Bce,Ban,Bth),Bwe)),(Bha,Bcl))";
# OR in a more readble format: 
#	( 
#	(
#	 	( ( ( (Bmo,Bsu),Bam), Bli), Bpu ), 
#		( (Bce,Ban, Bth),Bwe)
#	), (Bha, Bcl) 
#	);




open (IN, "<$clusterfl" ); my @clusterdirs = <IN>; close (IN);

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  open (TRE, ">spec.tre" ); print TRE $tre."\;\n"; close (TRE);
  if ($debug ) { system("cat spec.tre"); }

  open (OUT, ">my.notung.batch" );
  print OUT "spec.tre\n";
  if ( -e "prot-bacillus-paup-nj-BhaRoot1.nwk" ) {
     print OUT "prot-bacillus-paup-nj-BhaRoot1.nwk\n";
  } 
  if ( -e "prot-bacillus-paup-nj-BclRoot1.nwk" ) {
     print OUT "prot-bacillus-paup-nj-BclRoot1.nwk\n";
  }
  close (OUT);
    
  my $cmd1 = "java -jar /home/hqin/tools/Notung-2.5/Notung-2.5.jar -b my.notung.batch --reconcile --speciestag postfix --treeoutput newick";

  system( $cmd1 );

  chdir "..";
}

exit;

#
# END
#

