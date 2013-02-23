BEGIN { push(@INC,"/home/hqin/lib/perl/", "/Users/hongqin/lib/perl");};

my $debug= 1 ;

use strict; use warnings; 
use Util; 
use Bio::TreeIO;
my $HOME = $ENV{'HOME'}; 

my $redoalign = 0; #alignment is done; 

my @models = qw(H0 H2C1S1); 
#my $modeldir = "H2C1S1"; #This is the model
my @TreeFiles = qw(bacillus.treeH0.nwk bacillus.treeH2C1S1.nwk); 

my $idfl   = "$HOME/projects/coat.protein07/key.data/_id2genomes.csv";
my $specfl = "$HOME/projects/coat.protein07/key.data/_genome2species.csv";
my $homedir = "$HOME/coat.protein07/ortholog.analysis/coatpaml/test.nearlyAll.Jan18,2011";

my $tmpfile = "/tmp/___tmpncoat.ajdsfl.codeml.txt"; 
system("ls -d wkdir/BG* > $tmpfile");
open (IN,"<$tmpfile"); my @lines = <IN>; close (IN); 
my @bgs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($tmp, $bg, @res) = split( /\//, $line);
  push (@bgs, $bg); 
}

if($debug>9) { @bgs = splice(@bgs, 0, 2); }

my $count = 0;
foreach my $mydir ( @bgs ) {
  chdir $homedir; 
  chomp $mydir;
  chdir "wkdir/$mydir";
  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 
 
  if ($redoalign) {
   ### align CDS
   my $cmd1 =  "align_cds_by_protein_bioperl.00.pl -i _cds.spec.fna -o tmpcds.phylip -af phylip ";
   system( $cmd1 );

   ### indicate the alignment is interleaved. (required by codeml)
   open(IN2, "<tmpcds.phylip"); my @lines2 = <IN2>;    close (IN2);
   open(OUT, ">cds2.phylip");
     my $firstline = shift @lines2;
     chomp $firstline;
     print OUT $firstline."  I\n";
     foreach my $line (@lines2) { print OUT $line; }
   close(OUT);
  }#redo alignment

  foreach my $i ( 0 .. $#models ) {
  my $modeldir = $models[$i]; 
   if ( ! (-e $modeldir) ) { mkdir $modeldir; }
   chdir $modeldir;
   system("cp ../cds2.phylip .");
   system("cp ../bacillus.treeH2C1S1.nwk $TreeFiles[$i]");
   chdir '..'; 
  }
}

exit;

#
# END
#

########################


