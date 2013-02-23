# build on _1.x _2.x.pl scripts
use strict; use warnings; 
my $debug= 0 ;

BEGIN {    unshift(@INC,"/home/hqin/lib/perl/", "/Users/hongqin/lib/perl"); };
use Util; 

my $HOME = $ENV{'HOME'}; 

#my $inprofilefile = 'RBH.4hits.nonCE.21Jan2011.tab';
#my $inprofilefile = 'specKey-RBH.4hits.nonCE.23Jan2011.tab';
my $inTreefile = '_taxon-tree.hash.b.tab';

my $homedir = "$HOME/coat.protein07/ortholog.analysis/coatpaml/test.nearlyAll.Jan18,2011";
my $idfl   = "$HOME/projects/coat.protein07/key.data/_id2genomes.csv";
my $specfl = "$HOME/projects/coat.protein07/key.data/_genome2species.csv";

#################
my %orf2genomes = ();
open (IN, "<$idfl" );
while (my $line = <IN> ) {
  chomp $line;
  my ($orf, $g, @res) = split (/\t/, $line);
  $orf2genomes{$orf} = $g;
}
close (IN);

my %genome2species = ();
open (IN, "<$specfl");
my $count = 1;
while (my $line = <IN> ) {
  if ( $count >= 2 ) {  #skip the header
  chomp $line;
  my ($genome, $spec, $flag, @res) = split (/\t/, $line);
  if ($debug) { print "[$genome] [$spec] [$flag]\t"; }
  $genome2species{ $genome } = { 
   'SPEC' => $spec,
   'FLAG' =>  $flag };
  } else { $count ++;
  } 
}
close (IN);

####################

### species trees
my %taxons2tree; 

open (IN, "<$inTreefile");
while (my $line = <IN> ) {
  chomp $line; 
  my ($taxons, $tree) = split(/\t/, $line);
  $taxons2tree{$taxons} = $tree; 
}
close(IN);
if($debug) { showHashTable(\%taxons2tree)}

####################
#### read ortholog profile, match them to truncated species trees 
#open (IN, "<$inprofilefile"); my @lines = <IN>; close (IN); 
#shift @lines; #remove header
#my @inSpecies = qw(Bsu Bmo Bam Bli Bpu Ban Bce Bth Bwe Bcl Bha);

my $tmpfile = "/tmp/__coatJan24,2011.asdfat.txt";
system( "ls -d wkdir/BG*/ > $tmpfile");

open (IN, "<$tmpfile"); my @lines = <IN>; close (IN); 
if ($debug>9) { @lines = splice(@lines, 0, 10); }

my @bgs = ();

foreach my $line (@lines) {
  chomp $line;
  my ($tmp, $bg,  @res) = split( /\//, $line);
     push (@bgs, $bg);
}

open (DUMP, ">__dump.b1.2.txt");

my $count = 1; 
my @clusterdirs; 
foreach my $bg (@bgs) {
   chdir $homedir; 
   chdir "wkdir/$bg"; 
   my $BG = $bg; #just lazy. keep old codes. 

   print "\n\nCluster:$count I am now working on [$bg]::\n";   

   my @orthoHits = ();
   open (IN, "<bacillus.ids.tab");
   while (my $line = <IN> ) {
     chomp $line;
     push (@orthoHits, $line); 
   }
   close(IN);

   if ($debug>10) { for(my $i=0; $i<=$#orthoHits; $i++) { print "$i-->$orthoHits[$i]\t"; } print "\n"; }
   if ($debug) { print "count=$count-----$BG:(@orthoHits)\n"; }
   my @currentHitSpecies =();
   my @notNAHits = ();
   for(my $i=0; $i<= $#orthoHits; $i++){
      if ((defined $orthoHits[$i]) & ($orthoHits[$i] !~ m/^\s*$/) & ($orthoHits[$i] !~ m/^\s*NA\s*$/) ) { #neither empty nor NA
	push( @currentHitSpecies, $genome2species{ $orf2genomes{$orthoHits[$i]} }->{SPEC});
	push (@notNAHits, $orthoHits[$i]); 
      } else {
	if ($debug > 5) {	print "Warning: $BG, $i ($orthoHits[$i])\n"; }
      }	
   }
   my @orderedSHits = order_big2small_array(@currentHitSpecies);
   if($debug) {print "$BG hits:(@currentHitSpecies)\n";}
   if($debug) {print "$BG hits:(@orderedSHits)\n";}
   my $key = join( ':', @orderedSHits ); 
   if ( ! defined $taxons2tree{$key} ) { print DUMP "no tree for $BG\t$key\n"}
       
      open (OUT, ">spec.key.tab"); print OUT $key."\n";           close (OUT); 

      #write trees
      #my $key = join( ':', @orderedSHits ); 
      open(OUT, ">bacillus.treeH2C1S1.nwk"); 
      if($debug) { print "Tree: "; }
      if (defined $taxons2tree{$key} ) {
         print OUT $taxons2tree{$key}."\n";
	 if($debug) { print "$BG->$key -> $taxons2tree{$key}\n"; }
       } else {
	 print "Warning: no tree for $BG: $key"; 
       }
      close(OUT);  
       
   $count ++; 
   #paml loop

}

close (DUMP);

exit;

#
# END
#

########################

   #my ($Bsu, $Bmo, $Bam, $Bli, $Bpu) = ($els[15], $els[16], $els[17], $els[18], $els[19]); 
   #my ($Ban, $Bce79, $Bth, $Bwe, $Bcl, $Bha) = ($els[20],$els[21],$els[22],$els[23],$els[24]); 
