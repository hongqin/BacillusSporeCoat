use strict; use warnings; 
use lib "/Users/hongqin/lib/perl"; use Util; 

my $HOME = '/Users/hongqin';
my $debug= 1 ;

my $inprofilefile = '_coat.profile.for.truncated.codeml.test.Jan18,2011.tab';
my $homedir = "$HOME/coat.protein07/ortholog.analysis/coatpaml/test.nearlyAll.Jan18,2011";

my $modeldir = "H1C"; #This is the model

open (IN,"<$inprofilefile"); my @lines = <IN>; close (IN); shift @lines;
my @bgs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($bg, $num, $paml, @res) = split( /\t/, $line);
  if ($paml =~ m/\s*YES\s*/ ) {
     push (@bgs, $bg);
  }
}

if($debug>9) { @bgs = splice(@bgs, 0, 2);  } #if($debug>9) { @bgs = ($bgs[5], $bgs[7], $bgs[11]); }

my $count = 0;
foreach my $mydir ( @bgs ) {
  chdir $homedir; 
  chomp $mydir;
  chdir "wkdir/$mydir";
  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 
 
  ### align CDS
  my $cmd1 =  "align_cds_by_protein_bioperl.00.pl -i _cds.spec.fna -o tmpcds.phylip -af phylip ";
  system( $cmd1 );

  ### indicate the alignment is interleaved. (required by codeml)
  open(IN2, "<tmpcds.phylip"); my @lines = <IN2>;    close (IN2);

  open(OUT, ">cds2.phylip");
    my $firstline = shift @lines;
    chomp $firstline;
    print OUT $firstline."  I\n";
    foreach my $line (@lines) { print OUT $line; }
  close(OUT);

  if ( ! (-e $modeldir) ) { mkdir $modeldir; }
  chdir $modeldir;
  system("cp ../cds2.phylip .");
  system("cp ../bacillus.tree1.nwk .");


}

exit;

#
# END
#

########################


