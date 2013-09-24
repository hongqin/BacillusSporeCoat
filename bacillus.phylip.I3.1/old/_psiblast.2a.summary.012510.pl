# Jan 25, 2010 psiblast summary 
# 1 split run
# 2 collect unique IDs of runs
# 3 compare uniques ID between runs

use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 10 ;

my $clusterfl = "infile-BGclusters.txt";

open (IN, "<$clusterfl" ); my @clusterdirs = <IN>; close (IN);

my %RUNs = ();
my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;    my $myBG = $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  open (IN, "<_$myBG.psi.m9.csv"); my @lines = <IN>; close (IN); 
  my @headers = grep( '# BLASTP 2.2.13 [Nov-27-2005]', @lines );
  my $iteractions = $#headers + 1; 

  my $alllongline = join( '', @lines );
  my @runlongline = split( "# BLASTP", $alllongline );
  shift @runlongline; #first element is empty
  print "There are ". ($#runlongline + 1). " runs\n";
  
  for my $itr  ( 0 .. $#runlongline ) {
    print "####   \$itr = $itr  #####\n ";
    my @runlines = split ( "\n", $runlongline[$itr] ); 
    shift @runlines; 
    my @idlines = grep ( /^BG/, @runlines );
    my %hits = ();
    foreach my $idline (@idlines) {
      my @els = split (/\t/, $idline );
      $hits {$els[1]} ++; 
    }

    if ($debug>9) { showHashTable(\%hits); }
    my @tmpkeys = keys (%hits); 
    $RUNs{ ($itr + 1) } = {
        'HITS' => join( ' ', keys (%hits) ),
        'NUMofHITS' => $#tmpkeys + 1,
        'DIFF' => 0,
      }; 

    if ( $itr > 0 ) { #compare to previous steps
      my $diff = $RUNs{ ($itr+1) } -> {NUMofHITS} - $RUNs{ $itr } -> {NUMofHITS}; 
      $RUNs{ ($itr+1) } -> {DIFF} = $diff;
    }
    if ($debug>9) { showHashTable( $RUNs{($itr+1)} ); }
  }

  #print report
  #system("rm *report* ");
  open(OUT, ">_$myBG.psi.report.csv"); 
     print OUT "iterations\tNumOfHits\tDifference\n";
  foreach my $itr ( sort( keys(%RUNs) ) ) {
     print OUT $itr."\t". $RUNs{ $itr }->{NUMofHITS}; 
     print OUT "\t". $RUNs{ $itr }->{DIFF}; ;
     print OUT "\n"; 
  }
  close(OUT);
  chdir "..";
}

exit;

#
# END
#

