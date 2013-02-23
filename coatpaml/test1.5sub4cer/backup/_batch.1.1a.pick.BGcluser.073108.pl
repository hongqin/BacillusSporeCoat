#073108 pick BG clusters for paml 
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; 

my $debug= 2 ;

my $inclusterfile = '_073108.test.csv';
my $sourcedir = '/home/hqin/projects/coat.protein07/orthog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1/';

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
  if ( !( -e $mydir ) ) {  mkdir $mydir; }
  chdir $mydir;
  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 
  system("cp $sourcedir/$mydir/_ids.bacillus.csv .");

  if( -e '_ids.bacillus.curated.csv' ) {
    system(" cut -f 2 _ids.bacillus.curated.csv > _ids.tab" );
  } else {
    system(" cut -f 2 _ids.bacillus.csv > _ids.tab");
  }

  chdir "..";
}

exit;

#
# END
#

########################

