use strict; use warnings; 
use lib "/Users/hongqin/lib/perl"; use Util; 

my $HOME = '/Users/hongqin';
my $debug= 1 ;
my $inprofilefile = '_coat.profile.for.truncated.codeml.test.Jan18,2011.tab';
my $inTreefile = '_taxon-tree.hash.tab';

my $sourcedir = "$HOME/coat.protein07/orthog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1/";
my $homedir = "$HOME/coat.protein07/ortholog.analysis/coatpaml/test.nearlyAll.Jan18,2011";

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

#### read ortholog profile, match them to truncated species trees 
open (IN, "<$inprofilefile"); my @lines = <IN>; close (IN); 
shift @lines; #remove header
if ($debug>9) { @lines = splice(@lines, 0, 10); }

my @inSpecies = qw(Bsu Bmo Bam Bli Bpu Ban Bce Bth Bwe Bcl Bha);

my $count = 1; 
my @clusterdirs; 
foreach my $line (@lines) {
   my @els = split( /\t/, $line ); if ($debug>10) { for(my $i=0; $i<=$#els; $i++) { print "$i-->$els[$i]\t"; } print "\n"; }
   my ($BG, $paml) = ($els[0], $els[2] );
   my @orthoHits = splice( @els, 15, 11); 
   if ($debug) { print "$count-----$BG:(@orthoHits)\n"; }
   my @currentHitSpecies =();
   my @notNAHits = ();
   for(my $i=0; $i<= $#orthoHits; $i++){
      if ((defined $orthoHits[$i]) & ($orthoHits[$i] !~ m/^\s*$/) & ($orthoHits[$i] !~ m/^\s*NA\s*$/) ) { #neight empty nor NA
          push (@currentHitSpecies, $inSpecies[$i]);
          push (@notNAHits, $orthoHits[$i]);
      } else {
	if ($debug > 5) {	print "Warning: $BG, $i ($orthoHits[$i])\n"; }
      }	
   }
   my @orderedSHits = order_big2small_array(@currentHitSpecies);
   if($debug) {print "$BG hits:(@currentHitSpecies)\n";}
   if($debug) {print "$BG hits:(@orderedSHits)\n";}

   if ($paml =~ m/\s*YES\s*/ ) { #start here. 
      chdir $homedir; 
      my $mydir = "wkdir/$BG"; 
      if ( !( -e $mydir ) ) {  mkdir $mydir; }
      chdir $mydir;
      print "Cluster:$count I am now working on [$mydir]::\n"; 
      
      #write @orthoHits, pick the first seq if there are multiple ones
      #my $flag = array_2_file(\"bacillus.ids.tab", \@notNAHits, \"\n");
      open (OUT, ">bacillus.ids.tab");
      foreach my $item (@notNAHits) {
	 $item =~ s/^\s+//o; 
	 my @els = split( /\s+/, $item );
	 print OUT $els[0]."\n"; 
      }
      close (OUT); 

      #write trees
      my $key = join( ':', @orderedSHits ); 
      open(OUT, ">bacillus.tree1.nwk"); 
      if (defined $taxons2tree{$key} ) {
         print OUT $taxons2tree{$key}."\n";
	 if($debug) { print "$BG->$key -> $taxons2tree{$key}\n"; }
       } else {
	 print "Warning: no tree for $BG: $key"; 
       }
      close(OUT);  
       
   }

  $count ++; 
}


$count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  if ( !( -e $mydir ) ) {  mkdir $mydir; }
  chdir $mydir;
 
#  if( -e '_ids.bacillus.curated.csv' ) {
#    system(" cut -f 2 _ids.bacillus.curated.csv > _ids.tab" );
#  } else {
#    system(" cut -f 2 _ids.bacillus.csv > _ids.tab");
#  }

  chdir "..";
}

exit;

#
# END
#

########################

   #my ($Bsu, $Bmo, $Bam, $Bli, $Bpu) = ($els[15], $els[16], $els[17], $els[18], $els[19]); 
   #my ($Ban, $Bce79, $Bth, $Bwe, $Bcl, $Bha) = ($els[20],$els[21],$els[22],$els[23],$els[24]); 
