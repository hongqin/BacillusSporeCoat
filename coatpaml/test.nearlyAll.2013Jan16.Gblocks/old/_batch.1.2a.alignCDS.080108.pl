#073108 pick BG clusters for paml 
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; 

my $debug= 2 ;

my $inclusterfile = '_073108.test.csv';
open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); shift @lines;
my @clusterdirs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($bg, @res) = split( /\t/, $line);
  push (@clusterdirs, $bg);
}

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 

  if ($debug) {    system( "pwd" );  }

  ### align CDS
  my $cmd1 =  "align_cds_by_protein_bioperl.00.pl -i _cds.spec.fna -o tmpcds.phylip -af phylip ";
  system( $cmd1 );

  ### indicate it is interleaved.
  open(IN2, "<tmpcds.phylip"); my @lines = <IN2>; 
  close (IN2);

  open(OUT, ">cds.phylip");
    my $firstline = shift @lines;
    chomp $firstline;
    print OUT $firstline."  I\n";
    foreach my $line (@lines) { print OUT $line; }
  close(OUT);

  chdir "..";
}

exit;

#
# END
#

########################

