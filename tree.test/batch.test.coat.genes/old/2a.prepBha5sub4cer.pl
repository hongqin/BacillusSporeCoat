 
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 1 ;

my $inprofilefile = '_10g.coat.profile.011009.csv';

system("pwd > __mydir.txt"); open(IN,"<__mydir.txt"); my @tmplines=<IN>; close(IN);
my $root = shift @tmplines; chomp $root; 

my $wkdir = 'Bha5sub4cer';

open (IN,"<$inprofilefile"); my @lines = <IN>; close (IN);
my $header = shift @lines; 

my $long_choices = "Bha Bpu Bli Bam Bsu Bmo Bwe Bce Ban Bth";
my @choices = split( /\s+/, $long_choices);

my $count = 0;
foreach my $line ( @lines ) { 
  chdir $root;
  my @els = split( /\t/, $line );
  my $mydir = $els[0];
  $mydir =~ s/\"//g;
  chdir $mydir;

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
  system("rm tmpcds.phylip");
  #chdir "..";
}

exit;

#
# END
#

