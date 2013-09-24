# v2, Feb 10, 2010, parse profiles
# v1. Jan 25, 2010 psiblast summary 
# 1 split run
# 2 collect unique IDs of runs
# 3 compare uniques ID between runs

use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 10 ;

my $clusterfl = "infile-BGclusters.all.txt";
open (IN, "<$clusterfl" ); my @clusterdirs = <IN>; close (IN);

my $idfl   = "/home/hqin/projects/coat.protein07/key.data/_id2genomes.csv";
my $specfl = "/home/hqin/projects/coat.protein07/key.data/_genome2speciesv2.csv";

my @orders  = ('Bsu','Bmo','Bam','Bli','Bpu','Ban','Bce87', 'Bce79', 'Bce3L','Bth','Bwe','Bcl','Bha');

#system( "ls -d BG* > /tmp/_tmpbg.txt" );
#open(IN, "</tmp/_tmpbg.txt");
#my @clusterdirs = ();
#while(my $line =<IN>) {
# chomp $line;
#  push (@clusterdirs, $line);
#}
#close(IN);

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

my %RUNs = ();
$count = 0;
my %profile = ();

#my @orders  = ('Bsu','Bmo','Bam','Bli','Bpu','Ban','Bce87', 'Bce79', 'Bce3L','Bth','Bwe','Bcl','Bha');
open (PROFILE, ">_profile.psi.021010.csv");
print PROFILE "coat"; foreach my $spec (@orders) { print PROFILE "\t$spec"; } 
print PROFILE "\tNumOfHits\tNumOfGenomeHits\n";

foreach my $mydir ( @clusterdirs ) {
  %RUNs = (); #021010 change. This cause a bug in previous results. 
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
  
  for my $itr  ( 0 .. $#runlongline ) {
    print "####   \$itr =".($itr +1)."  #####\n ";
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
  open(OUT, ">_$myBG.psi.bacilli.report.csv"); 
     print OUT "iterations\tNumOfHits\tDifference\n";
  foreach my $itr ( sort( keys(%RUNs) ) ) {
     print OUT $itr."\t". $RUNs{ $itr }->{NUMofHITS}; 
     print OUT "\t". $RUNs{ $itr }->{DIFF}; ;
     print OUT "\n"; 
  }
  close(OUT);


  ###use $RUNS{1} for new profile, Feb 10, 2010, Ref: _summarize.profile.073108.pl 
  my @hits = split ( /\s+/, $RUNs{1}->{HITS} );
  my %_cur_profiles = ();
  foreach my $hit (@hits ) {
     my $myflag = $genome2species{ $id2genomes{$hit} }->{FLAG};
     if ($myflag eq 'Y' ) {
        $_cur_profiles{$genome2species{ $id2genomes{$hit} }->{SPEC}} .= $hit.' ';
     }
  }
  if ($debug > 5) { 
	print "##\%_cur_profiles::\n";
	showHashTable( \%_cur_profiles); 
  }
  close (IN);

  print PROFILE "$myBG";
  foreach my $spec (@orders) { 
    if (exists $_cur_profiles{$spec} ) {
  	print PROFILE "\t".$_cur_profiles{$spec}; 
    } else {
        print PROFILE "\tNA";
    } 
  } 
  print PROFILE "\t". ($#hits + 1);
  my @genomehits = keys( %_cur_profiles );
  print PROFILE "\t". ($#genomehits + 1) ."\n"; 
  ### END, print PROFILE for $myBD;

  chdir "..";
}

close (PROFILE); 



exit;

#
# END
#

