	# 062208 pick entries in the Bacillus clade, add labels for NOTUNG
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 2 ;

#my $idfl   = "/home/hqin/projects/coat.protein07/key.data/_id2genomes.csv";
#my $specfl = "/home/hqin/projects/coat.protein07/key.data/_genome2species.csv";
#my $clusterfl = "infile-BGclusters.txt";

my $__wrong = "Bmo2454";
my %wrongids = ();
my @els = split( /\s+/, $__wrong);
foreach my $el (@els ) {
   $wrongids{$el} ++;
}
showHashTable( \%wrongids );

### find out the cluster-directories
my $tmpfl = "/tmp/_coatBGclus.062808.txt";
#   system( "ls -d BG* | head -n 2 > $tmpfl " );
system( "ls -d BG*  > $tmpfl " );
open (IN, "<$tmpfl" ); my @clusterdirs = <IN>; close (IN);
system( "rm $tmpfl" );

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  if ( ! ( -e "ids.tab.old" ) ) {
    system( "cp ids.tab ids.tab.old");
  }

  open (IN, "<ids.tab.old"); my @lines=<IN>; close (IN);

  open ( OUT, ">ids.tab" );
  foreach my $line (@lines) {
  	chomp $line;
  	$line =~ s/\s+//g;
  	if ( !(exists($wrongids{$line}) ) ){
  		print OUT $line."\n";
  	} else {
     	    	print "!!! I have removed $line from $mydir\n\n";
        } 
  }
  close (OUT);
  chdir "..";
}

exit;

#
# END
#

