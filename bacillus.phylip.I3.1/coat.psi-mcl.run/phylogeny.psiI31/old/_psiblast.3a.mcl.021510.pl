# v1. Feb 15, 2010, MCL on psi-blast output
# 1 take the first run, put into a merged file (for deal with duplicates)
# 2 run mcl

use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 9 ;

my $homedir = "/home/hqin/projects/coat.protein07/ortholog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1";

my $clusterfl = "infile-BGclusters.all.txt";
open (IN, "<$clusterfl" ); my @clusterdirs = <IN>; close (IN);

if ($debug > 9 ) { @clusterdirs = splice( @clusterdirs, 0, 2); }

my $idfl   = "/home/hqin/projects/coat.protein07/key.data/_id2genomes.csv";
my $specfl = "/home/hqin/projects/coat.protein07/key.data/_genome2speciesv2.csv";

my @orders  = ('Bsu','Bmo','Bam','Bli','Bpu','Ban','Bce87', 'Bce79', 'Bce3L','Bth','Bwe','Bcl','Bha');

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
  my ($genome, $spec, $flag, $spec2, $flag2, @res) = split (/\t/, $line);
  #print "[$genome] [$spec] [$flag] [$spec2] [$flag2]\t";
  print "[$genome] [$spec2] [$flag2]\t";
  $genome2species{ $genome } = { 
   'SPEC' =>  $spec2,
   'FLAG' =>  $flag2 };
  } else { $count ++;
  } 
}
close (IN);

print "\n*************************\n";
foreach my $key (sort( keys %genome2species ) ) {
  print $genome2species{$key}->{SPEC} ."::".$genome2species{$key}->{FLAG} . "\n";
}
print "\n*************************\n";

#I parameter for mcl runs
#my @Is = ( 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,21,22,23,24, 25,26,27,28,29, 30,31,32,33,34,35,36,37,38,39,40, 50, 60 ); #done
  my $start = 11; my @Is = ( 60, 70, 80); 
  while ( $start< 50 ) {
    push (@Is, $start);
    $start ++;
  }


#my %RUNs = ();
$count = 0;
my %profile = ();

#my @orders  = ('Bsu','Bmo','Bam','Bli','Bpu','Ban','Bce87', 'Bce79', 'Bce3L','Bth','Bwe','Bcl','Bha');
#open (PROFILE, ">_profile.psi.021010.csv");
#print PROFILE "coat"; foreach my $spec (@orders) { print PROFILE "\t$spec"; } 
#print PROFILE "\tNumOfHits\tNumOfGenomeHits\n";

open (ONE, ">coat.psi-mcl.run/merged.coat.psiblast.run1.csv");

foreach my $mydir ( @clusterdirs ) {
  #%RUNs = (); #021010 change. This cause a bug in previous results. 
  chomp $mydir;    my $myBG = $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  open (IN, "<_$myBG.psi.baci.m9.csv"); my @lines = <IN>; close (IN); 
  my @headers = grep( '# BLASTP 2.2.13 [Nov-27-2005]', @lines );
  my $iteractions = $#headers + 1; 

  my $alllongline = join( '', @lines );
  my @runlongline = split( "# BLASTP", $alllongline );
  shift @runlongline; #first element is empty
  print "There are ". ($#runlongline + 1). " runs\n";
  
  #parse out the first psiblast run
  my $itr = 0; 
  print "####   \$itr =".($itr +1)."  #####\n ";
  my @runlines = split ( "\n", $runlongline[$itr] ); 
  shift @runlines; 
  my @idlines = grep ( /^BG/, @runlines );

  #output to ~/mcl-psi/_**run1.csv
  #if( ! -e 'mcl-psi') { system( "mkdir mcl-psi" ); }
  #chdir 'mcl-psi'; 
  foreach my $idline (@idlines) {
     print ONE $idline."\n"; 
  }  
   
  chdir "..";
}

close (ONE); 


#run mcl   
 chdir $homedir;
 chdir "coat.psi-mcl.run";  
 foreach my $I ( @Is ) {
   #my $cmnd = "cut -f 1,2,11 merged.coat.psiblast.run1.csv | mcl - --abc -I ". $I/10 ." -o _coat.psi.mcl.I$I.clus "; #E-value
   my $cmnd = "cut -f 1,2,12 merged.coat.psiblast.run1.csv | mcl - --abc -I ". $I/10 ." -o _coat.psi.mcl.I$I.clus "; # score
   print $cmnd."\n";
   system( "$cmnd");
 }


exit;
#
# END
#