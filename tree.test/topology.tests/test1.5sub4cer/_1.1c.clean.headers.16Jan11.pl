use strict; use warnings; 
use lib '/Users/hongqin/lib/perl'; use Util; 

my $debug= 2 ;
my $inclusterfile = '_test.csv';


open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); shift @lines;
my @clusterdirs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($num, $bg, @res) = split( /\t/, $line);
  push (@clusterdirs, $bg);
}

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir "wkdir/$mydir";
  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 

  system( "clean_fasta_headers.01.pl -if ids.faa -of prot.faa -pos 0 " );
  system( "clean_fasta_headers.01.pl -if ids.fna -of cds.fna  -pos 0 " );

  chdir "../..";
}

exit;

#
# END
#

########################

