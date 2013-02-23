 
use strict; use warnings; 
use lib '/Users/hongqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 20 ;

my $inclusterfile = '_073108.test.csv';

open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); shift @lines;
my @clusterdirs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($bg, @res) = split( /\t/, $line);
  push (@clusterdirs, $bg);
}

my $root = "/home/hqin/projects/coat.protein07/orthog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1/paml/test1.5sub4cer/";
#my $wkdir = "Bha5sub4cer090108";
my $wkdir = "Bha5sub4cer120908";  #120908

my   $long_choices = "Bha Bpu Bli Bam Bmo Bsu Bwe Bce Ban Bth";
 #my $long_choices = "Bha Bpu Bli Bam Bmo Bsu Bwe Bce Ban";
my @choices = split( /\s+/, $long_choices);


my $count = 0;

foreach my $mydir ( @clusterdirs ) {
  if (($debug>9) &&($count>2) ) { exit; }
 
  chdir $root;   chomp $mydir;   chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  if ( ! (-e $wkdir) ) { mkdir $wkdir; }
  chdir $wkdir;
  open (OUT,">_subids.tab");
  foreach my $id (@choices) { print OUT "$id\n"; } close (OUT);
  system( "pick_fasta_records_by_ids.1.pl -if ../_cds.spec.fna -of cdssub.fna -id _subids.tab" );

  ### align CDS
  my $cmd1 =  "align_cds_by_protein_bioperl.00.pl -i cdssub.fna -o tmpcds.phylip -af phylip ";
  system( $cmd1 );

  ### indicate it is interleaved.
  open(IN2, "<tmpcds.phylip"); my @lines = <IN2>;    close (IN2);

  open(OUT, ">cds2.phylip");
    my $firstline = shift @lines;
    chomp $firstline;
    print OUT $firstline."  I\n";
    foreach my $line (@lines) { print OUT $line; }
  close(OUT);

  #chdir "..";
}

exit;

#
# END
#

