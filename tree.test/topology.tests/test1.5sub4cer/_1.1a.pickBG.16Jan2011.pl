#Jan 11,2011
#090508 pick BG clusters for paml 
use strict; use warnings; 
#use lib '/home/hqin/lib/perl'; use Util; 
use lib '/Users/hongqin/lib/perl'; use Util; 

my $debug= 2 ;

my $inclusterfile = '_test.csv';
#my $sourcedir = '/Users/hongqin/projects/coat.protein07/orthog.analysis/essetial.mcl.rbh.090208/bacillus.phylip.I3.1/wkdir';
my $sourcedir = '/Users/hongqin/projects/coat.protein07/ortholog.analysis/essetial.mcl.rbh.090208/bacillus.phylip.I3.1/wkdir';

open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); shift @lines;
my @clusterdirs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($num, $bg, @res) = split( /\t/, $line);
  print "[$num][$bg]------";
  push (@clusterdirs, $bg);
}


my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  if ( !( -e "wkdir/$mydir" ) ) {  mkdir "wkdir/$mydir"; }
  chdir "wkdir/$mydir";
  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 
  system("cp $sourcedir/$mydir/_ids.bacillus.csv .");

  if( -e '_ids.bacillus.curated.csv' ) {
    system(" cut -f 2 _ids.bacillus.curated.csv > _ids.tab" );
  } else {
    system(" cut -f 2 _ids.bacillus.csv > _ids.tab");
  }

  chdir "../..";
}

exit;

#
# END
#

########################

