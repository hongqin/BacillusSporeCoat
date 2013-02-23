#090508 pick BG clusters for paml 
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; 

my $debug= 0 ;

my $inclusterfile = 'my.curated.essen.gene.csv.112608b';
my $sourcedir = '/home/hqin/projects/coat.protein07/orthog.analysis/essetial.mcl.rbh.090208/bacillus.phylip.I3.1/wkdir/';

open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); 
my @clusterdirs = ();
foreach my $line (@lines) {
  chomp $line;
  $line =~ s/^\s+//o; 
  my ($nu, $bg, @res) = split( /\s+/, $line);
  push (@clusterdirs, $bg);
}

if( $debug > 9 ) { @clusterdirs = splice(@clusterdirs, 1, 3); } 

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


