# Jan 25, 2010 psiblast
# plan 
# 1. muscle prot.faa, 
# 2. pick query sequence
# 3. psiblast 
# 4. check ID changes

use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 6 ;

my $clusterfl = "infile-BGclusters.all.txt"; #021010 change, all coat clsuters 

open (IN, "<$clusterfl" ); my @clusterdirs = <IN>; close (IN);
if($debug>10) { @clusterdirs = splice( @clusterdirs, 0, 2); }

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;    my $myBG = $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  my $cmd1  =  "muscle -in prot.faa -out /tmp/_$myBG.aln -clw -stable";
  system( $cmd1 );

  #reformat the alignment
  open(TMPIN, "</tmp/_$myBG.aln"); my @tmplines = <TMPIN>; close (TMPIN); 
  shift @tmplines; #remove the first three lines
  shift @tmplines; 
  shift @tmplines; 

  open(TMPOUT, ">_$myBG.psi.aln"); #This alignment is for blastpgp
  foreach my $line (@tmplines) {
    print TMPOUT $line;
  }
  close (TMPOUT); 

  my $cmd2 = ''; my $cmd3 = "";
  #pick query sequence
  if ($myBG =~ /-/ ) {#more than two IDs
    my ($id1, $id2, @res) = split( /-/, $myBG); 
    $cmd1 = "pick_fasta_records_by_ids.1b.pl -if prot.faa -isd $id1 -o $id1.faa";   
    $cmd3 = "blastpgp -i $id1.faa -B _$myBG.psi.aln -j 5 -d /home/hqin/coat.protein07/genomes/_bacilli.genomes.012710.faa -m 9 -o _$myBG.psi.baci.m9.csv -e 0.0002";
  } else {
    $cmd1 = "pick_fasta_records_by_ids.1b.pl -if prot.faa -isd $myBG -o $myBG.faa";   
    $cmd3 = "blastpgp -i $myBG.faa -B _$myBG.psi.aln -j 5 -d /home/hqin/coat.protein07/genomes/_bacilli.genomes.012710.faa -m 9 -o _$myBG.psi.baci.m9.csv -e 0.0002";
  }
  system( $cmd1 );
  system( $cmd3 );

  #psi-blast
  #my $cmd = "blastpgp -i $myBG.faa -B _$myBG.psi.aln -j 3 -d /home/hqin/coat.protein07/genomes/_merged.all.062108.faa -m 9 -o _$myBG.psi.m9.csv -e 0.0002";

  chdir "..";
}

exit;

#
# END
#

