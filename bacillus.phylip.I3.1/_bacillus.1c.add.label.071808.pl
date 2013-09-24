# 071808  add species labels to _ids.bacillus.csv
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 1 ;

my $idfl   = "/home/hqin/projects/coat.protein07/key.data/_id2genomes.csv";
my $specfl = "/home/hqin/projects/coat.protein07/key.data/_genome2species.csv";
my $clusterfl = "infile-BGclusters.txt";
my $inIDcsv = '_ids.bacillus.csv';
my $outIDcsv = '_ids.bacillus.2.csv';

my $tre ="((((((Bmo,Bsu),Bam),Bli),Bpu),((Bce,Ban,Bth),Bwe)),(Bha,Bcl))";
# OR in a more readble format: 
#	( 
#	(
#	 	( ( ( (Bmo,Bsu),Bam), Bli), Bpu ), 
#		( (Bce,Ban, Bth),Bwe)
#	), (Bha, Bcl) 
#	);

### find out the cluster-directories
#my $tmpfl = "/tmp/_coatBGclus.062208.txt";
#   system( "ls -d BG* | head -n 2 > $tmpfl " );
#system( "ls -d BG*  > $tmpfl " );
open (IN, "<$clusterfl" ); my @clusterdirs = <IN>; close (IN);
#system( "rm $tmpfl" );

if($debug>9) { @clusterdirs = ( $clusterdirs[0], $clusterdirs[1] );  };

my %id2genomes = ();
open (IN, "<$idfl" );
while (my $line = <IN> ) {
  chomp $line;
  my ($id, $g, @res) = split (/\t/, $line);
  $id2genomes{$id} = $g;
}
close (IN);

my %genome2species = ();
open (IN, "<$specfl");
my $count = 1;
while (my $line = <IN> ) {
  if ( $count >= 2 ) {  #skip the header
  chomp $line;
  my ($genome, $spec, $flag, @res) = split (/\t/, $line);
  print "[$genome] [$spec] [$flag]\t";
  $genome2species{ $genome } = { 
   'SPEC' => $spec,
   'FLAG' =>  $flag };
  } else { $count ++;
  } 
}
close (IN);

print "\n*************************\n";
foreach my $key (keys %genome2species ) {
  print $genome2species{$key}->{SPEC} ."::".$genome2species{$key}->{FLAG} . "\t\t";
}
print "\n*************************\n";

$count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  open (IN, "<$inIDcsv");
  open ( OUT, ">$outIDcsv" );
  while (my $line =<IN> ) {
     chomp $line; 
     my ($phylip, $id) = split (/\t/, $line);
     my $myspec = $genome2species{ $id2genomes{$id} }->{SPEC};
     if ($debug) { print "id=[$id] spec =[$myspec]\n"; }
     print OUT "$phylip\t$id\t$id".'_'.$genome2species{ $id2genomes{$id} }->{SPEC}."\n";
  }
  close (OUT);
  chdir "..";
}

exit;

#
# END
#

