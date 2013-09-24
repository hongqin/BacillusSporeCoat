# 071208 protdist neighbour
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 0 ;

#my $idfl   = "/home/hqin/projects/coat.protein07/key.data/_id2genomes.csv";
#my $specfl = "/home/hqin/projects/coat.protein07/key.data/_genome2species.csv";
#my $clusterfl = "infile-BGclusters.txt";
my $clusterfl = "_curated.BGclusters.073008.txt";

#my $tre ="((((((Bmo,Bsu),Bam),Bli),Bpu),((Bce,Ban,Bth),Bwe)),(Bha,Bcl))";
# OR in a more readble format: 
#	( 
#	(
#	 	( ( ( (Bmo,Bsu),Bam), Bli), Bpu ), 
#		( (Bce,Ban, Bth),Bwe)
#	), (Bha, Bcl) 
#	);

open (IN, "<$clusterfl" ); 
my @clusterdirs = (); 
while( my $line = <IN> ) {
  my($bg, $fl, @res) = split( /\//, $line);
  push (@clusterdirs, $bg);
}
close (IN);

if($debug>9) { @clusterdirs = ( $clusterdirs[0], $clusterdirs[1] ); 
};

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  my $cmd = " clus2mol.pl prot-bacillus-aligned.10char.2.aln > _in1";
  system( $cmd );
  $cmd = "mol2phy.pl _in1 > infile  ";
  system( $cmd );
  system( "rm _in1" );

  open  (OUT, ">prompt.1.txt");   
  #print OUT "G\nI\n2\nY\n";   
  print OUT "I\nY\n";   
  close (OUT);  
  
  system( "rm outfile" );
  my $cmd2 = "protdist < prompt.1.txt ";
  system( $cmd2 );

  open  (OUT, ">prompt.2.txt");
  print OUT "Y\n";
  close (OUT);
  system( "rm outtree");
  system( "mv outfile infile ");
  system( "neighbor < prompt.2.txt");

  if ( -e "outtree" ) {
    system ("cat outtree");
    system ("mv outtree protbaci_nj.curated.nwk");
    
  } else { 
    system("touch __FAILED.prot-bacillus-neighbor.curated.txt");
  } 
  
  chdir "..";
}

exit;

#
# END
#

