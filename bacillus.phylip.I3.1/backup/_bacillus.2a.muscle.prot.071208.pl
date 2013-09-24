# 071208 muscle prot-bacillus.faa
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

  my $cmd1  =  "muscle -in prot-bacillus.10char.faa -out prot-bacillus-aligned.10char.fasta -stable ";
  system( $cmd1 );
  my $cmd1b =  "muscle -in prot-bacillus.10char.faa -out prot-bacillus-aligned.10char.aln   -stable -clw ";
  system( $cmd1b );

  chdir "..";
}

exit;

#
# END
#

