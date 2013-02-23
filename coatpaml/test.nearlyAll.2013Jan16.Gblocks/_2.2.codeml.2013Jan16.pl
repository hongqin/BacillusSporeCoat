# Jan 19, 2011 
# I decide to set "model=2" for all trees, including thoses that I did not specify omega. Whill this work?  NO!
#
# I decide to set "model=1" for all trees, including thoses that I did not specify omega. Whill this work?  Yes, for publication


use strict; use warnings; 
use lib "/Users/hongqin/lib/perl"; use Util; 

my $HOME = '/Users/hongqin';
my $debug= 0 ;

my $inprofilefile = '_coat.profile.for.truncated.codeml.test.Jan18,2011.tab';
#my $homedir = "$HOME/coat.protein07/ortholog.analysis/coatpaml/test.nearlyAll.Jan18,2011";
my $homedir = "$HOME/coat.protein07/ortholog.analysis/coatpaml/test.nearlyAll.2013Jan16.Gblocks";

my $modeldir = "H1C"; #This is the model !!!!!!!!!!

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
  chdir "wkdir/$mydir/$modeldir";
  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 
 
   #omega choice
   #if ($h eq 'H0' ) { $lines2[11] =  "       model = 0    \n"; #H0 one omega
   #} else { $lines2[11] =  "       model = 2    \n"; }  #alternatives, 2 omega

   open (OUT, ">codeml.ctl"); 
   my $ctlline = _get_ctl(); 
   print OUT $ctlline; 
   close(OUT);

   system( "codeml codeml.ctl");
}

exit; 

#
# END 
#

#########
sub _get_ctl {
my $long_ctl_line = "
      seqfile = cds3.phylip   * sequence data filename
     treefile = bacillus.tree1.nwk      * tree structure file name
      outfile = results.H1C.txt   * main result file name

        noisy = 9      * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1      * 1:detailed output
      runmode = 0      * 0:user defined tree

      seqtype = 1      * 1:codons
    CodonFreq = 2      * 0:equal, 1:F1X4, 2:F3X4, 3:F61

        model = 1      * 0:one omega ratio for all branches
                       * 1:separate omega for each branch
                       * 2:user specified dN/dS ratios for branches

      NSsites = 0      * 

        icode = 0      * 0:universal code

    fix_kappa = 0      * 1:kappa fixed, 0:kappa to be estimated
        kappa = 2      * initial or fixed kappa

    fix_omega = 0      * 1:omega fixed, 0:omega to be estimated 
        omega = 0.2    * initial omega

                       * comments:
                       * H0 in Table 3: model = 0
                       * H1 in Table 3: model = 2
                       * H2 in Table 3: model = 2
                       * H3 in Table 3: model = 2
"; 
 return $long_ctl_line; 
} 


#my @Hs = ('H0','H1a', 'H1b', 'H1c', 'H1d','H2a','H2b', 'H2c', 'H3a','H3b', 'H3c','H4a');
#my @Hs = ('H0','H1a', 'H1b', 'H1c', 'H1d','H2a','H2b', 'H2c', 'H3a', 'H3c','H4a');
#my @Hs = ('H1bwe','H1bce');

