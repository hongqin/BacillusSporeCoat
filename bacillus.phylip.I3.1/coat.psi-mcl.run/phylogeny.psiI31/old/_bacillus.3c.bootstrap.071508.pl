#071508 seqboot, protdist, neighbor, consensus
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 0;

#my $idfl   = "/home/hqin/projects/coat.protein07/key.data/_id2genomes.csv";
#my $specfl = "/home/hqin/projects/coat.protein07/key.data/_genome2species.csv";
#my $clusterfl = "infile-BGclusters.bootstrap.txt";
my $clusterfl = "infile-BGclusters.txt";

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
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 
 
  system( "rm outfile" );
  my $cmd = " clus2mol.pl prot-bacillus-aligned.10char.aln > _in1";
  system( $cmd );
  $cmd = "mol2phy.pl _in1 > infile  ";
  system( $cmd );
  system( "rm _in1" );
  system( "cp infile prot-bacillus.phylip.infile");

  open  (OUT, ">prompt.1.3c.txt");   
  print OUT "I\nY\n4333\n";   
  close (OUT);  
  system( "seqboot < prompt.1.3c.txt" ); 
  
  system( "mv outfile infile");
  open  (OUT, ">prompt.2.3c.txt");
  print OUT "I\nM\nD\n100\nY\n";
  close (OUT);
  system( "protdist < prompt.2.3c.txt" );

  system ("mv outfile infile");
  system (" rm outtree " );
  open (OUT, ">prompt.3.3c.txt");
  print OUT "M\n100\n433\nY\n";
  close (OUT);
  system( "neighbor < prompt.3.3c.txt");
  system( "mv outfile _outfile.neighbor" );

  system( "mv outtree intree" );
  open (OUT, ">prompt.4.3c.txt");
  print OUT "Y\n";
  close (OUT);
  system( "consense < prompt.4.3c.txt");

  if ( -e "outtree" ) {
    system ("cat outtree");
    system ("mv outtree protbaci_consense.nwk");
    
  } else { 
    system("touch __FAILED.prot-bacillus-consense.071508.txt");
  } 
  
  chdir "..";
}

exit;

#
# END
#

