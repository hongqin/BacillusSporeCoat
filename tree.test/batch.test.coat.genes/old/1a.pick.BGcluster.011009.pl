use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; 

my $debug= 2 ;

my $inprofilefile = '_10g.coat.profile.011009.csv';
my $wkdir = '/home/hqin/projects/coat.protein07/orthog.analysis/coatpaml/Bha10g'; 
 
open (IN,"<$inprofilefile"); my @lines = <IN>; close (IN);
my $header = shift @lines; 
my @els = split( /\t/, $header );
my $count2 = $#els; 
foreach my $i ( 0..$#els ) {
  print "$i $els[$i]\n";
}

my $count = 0;
foreach my $line ( @lines ) { 
  my @els = split( /\t/, $line );
  my $mydir = $els[0];
  $mydir =~ s/\"//g;

  chomp $mydir;
  if ( !( -e $mydir ) ) {  mkdir $mydir; }
  chdir $mydir;

  if($debug) { $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; }
  
  my @els2 = splice( @els, 9, 11);
  if($debug>9) {
   for( my $i=0; $i<=$#els; $i++ ) { 
     print "$i [$els[$i]]\n"; 
   }
  }

   open(OUT,">ids.tab"); 
   foreach my $el (@els2) {
    if ( $el !~ /NA/ )  {
      #if ( $el =~ /.\s+./ ) { print "\nMultiple entries\t[$el]\n";}; 
      $el =~ s/\s+//g;
      $el =~ s/\"//g;
      print OUT "$el\n";
    }
   }
   close(OUT);
  
  chdir $wkdir;
}

exit;

#
# END
#

########################


