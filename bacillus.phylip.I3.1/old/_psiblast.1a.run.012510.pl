# Jan 25, 2010 psiblast
# plan 
# 1. muscle prot.faa, 
# 2. pick query sequence
# 3. psiblast 
# 4. check ID changes

use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 2 ;

my $clusterfl = "infile-BGclusters.txt";

open (IN, "<$clusterfl" ); my @clusterdirs = <IN>; close (IN);

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

  #pick query sequence
  my $cmd = "pick_fasta_records_by_ids.1b.pl -if prot.faa -isd $myBG -o $myBG.faa";   
  system( $cmd );

  #psi-blast
  #system("rm __psi.out");
  my $cmd = "blastpgp -i $myBG.faa -B _$myBG.psi.aln -j 3 -d /home/hqin/coat.protein07/genomes/_merged.all.062108.faa -m 9 -o _$myBG.psi.m9.csv -e 0.0002";
  system( $cmd );

  chdir "..";
}

exit;

#
# END
#

