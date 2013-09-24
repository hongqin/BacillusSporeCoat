# 071208 muscle prot-bacillus.faa
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 2 ;

#my $idfl   = "/home/hqin/projects/coat.protein07/key.data/_id2genomes.csv";
#my $specfl = "/home/hqin/projects/coat.protein07/key.data/_genome2species.csv";
#my $clusterfl = "infile-BGclusters.txt";
my $clusterfl = "_curated.BGclusters.073008.txt";

my $tre ="((((((Bmo,Bsu),Bam),Bli),Bpu),((Bce,Ban,Bth),Bwe)),(Bha,Bcl))";
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

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  my $cmd1  =  "muscle -in prot-bacillus.10char.2.faa -out prot-bacillus-aligned.10char.2.fasta -stable ";
  system( $cmd1 );
  my $cmd1b =  "muscle -in prot-bacillus.10char.2.faa -out prot-bacillus-aligned.10char.2.aln   -stable -clw ";
  system( $cmd1b );

  chdir "..";
}

exit;

#
# END
#

