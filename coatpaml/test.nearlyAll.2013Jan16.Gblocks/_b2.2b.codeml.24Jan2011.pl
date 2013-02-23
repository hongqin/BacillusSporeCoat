# Jan 24, 2011
use strict; use warnings; 
BEGIN { push(@INC,"/home/hqin/lib/perl/", "/Users/hongqin/lib/perl");};

my $debug= 1 ;

use Util; 
#use Bio::TreeIO;
#use Bio::Tools::Phylo::PAML::Result;
#use Bio::Tools::Phylo::PAML;

my $HOME = $ENV{'HOME'}; 

my @HXs = qw(H0 H2C1S1); 
my @codemlmodels = qw(0 2); 
#my $modeldir = "H2C1S1"; #This is the model
#my @TreeFiles = qw(bacillus.treeH0.nwk bacillus.treeH2C1S1.nwk); 

#my $idfl   = "$HOME/projects/coat.protein07/key.data/_id2genomes.csv";
#my $specfl = "$HOME/projects/coat.protein07/key.data/_genome2species.csv";
my $homedir = "$HOME/coat.protein07/ortholog.analysis/coatpaml/test.nearlyAll.Jan18,2011";

my $tmpfile = "/tmp/___tmpcoatJ24,2011,llakjsfi.codeml.txt"; 
system("ls -d wkdir/BG* > $tmpfile");
open (IN,"<$tmpfile"); my @lines = <IN>; close (IN); 
my @bgs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($tmp, $bg, @res) = split( /\//, $line);
  push (@bgs, $bg); 
}

#if($debug>9) { @bgs = ($bgs[5], $bgs[7], $bgs[11]); }
if($debug>9) { @bgs = splice(@bgs, 0, 3); }

my $count = 0;
foreach my $bg ( @bgs ) {
  $count ++; print "Cluster:$count I am now working on [$bg]::\n"; 
 
  foreach my $i ( 0 .. $#HXs ) {
   my $modeldir = $HXs[$i]; 
   chdir $homedir."/wkdir/$bg/$modeldir"; 

   open (OUT, ">codeml.ctl"); 
   my $ctlline = _get_ctl($HXs[$i], $codemlmodels[$i]); 
   print OUT $ctlline; 
   close(OUT);

   system( "codeml codeml.ctl");
  }

}

exit; 

#
# END 
#

#########
sub _get_ctl {
  my ($HX, $model) = @_; 
my $long_ctl_line = "
      seqfile = cds2.phylip   * sequence data filename
     treefile = bacillus.treeHX.nwk      * tree structure file name
      outfile = results.HX.txt   * main result file name

        noisy = 9      * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1      * 1:detailed output
      runmode = 0      * 0:user defined tree

      seqtype = 1      * 1:codons
    CodonFreq = 2      * 0:equal, 1:F1X4, 2:F3X4, 3:F61

        model = X      * 0:one omega ratio for all branches
                       * 1:separate omega for each branch
                       * 2:user specified dN/dS ratios for branches

      NSsites = 0      * 

        icode = 0      * 0:universal code

    fix_kappa = 0      * 1:kappa fixed, 0:kappa to be estimated
        kappa = 2      * initial or fixed kappa

    fix_omega = 0      * 1:omega fixed, 0:omega to be estimated
        omega = 0.2    * initial omega
"; 
 my @lines = split( /\n/, $long_ctl_line ); 
 $lines[2] =~ s/HX/$HX/o; 
 $lines[3] =~ s/HX/$HX/o; 
 $lines[12] =~ s/X/$model/o; 
 
 return join( "\n", @lines); 
} 


#my @Hs = ('H0','H1a', 'H1b', 'H1c', 'H1d','H2a','H2b', 'H2c', 'H3a','H3b', 'H3c','H4a');
#my @Hs = ('H0','H1a', 'H1b', 'H1c', 'H1d','H2a','H2b', 'H2c', 'H3a', 'H3c','H4a');
#my @Hs = ('H1bwe','H1bce');

