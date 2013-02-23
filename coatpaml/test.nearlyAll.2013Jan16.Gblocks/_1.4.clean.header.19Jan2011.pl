use strict; use warnings; 
use lib "/Users/hongqin/lib/perl"; use Util; 

my $HOME = '/Users/hongqin';
my $debug= 1 ;

my $inprofilefile = '_coat.profile.for.truncated.codeml.test.Jan18,2011.tab';
#my $sourcedir = "$HOME/coat.protein07/orthog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1/";
my $homedir = "$HOME/coat.protein07/ortholog.analysis/coatpaml/test.nearlyAll.Jan18,2011";

open (IN,"<$inprofilefile"); my @lines = <IN>; close (IN); shift @lines;
my @bgs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($bg, $num, $paml, @res) = split( /\t/, $line);
  if ($paml =~ m/\s*YES\s*/ ) {
     push (@bgs, $bg);
  }
}

if($debug>9) { @bgs = ($bgs[5], $bgs[7], $bgs[11]); }

my $count = 0;
foreach my $mydir ( @bgs ) {
  chdir $homedir; 
  chomp $mydir;
  chdir "wkdir/$mydir";
  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 
  system( "clean_fasta_headers.01.pl -if bacillus.ids.faa -of bacillus.prot.faa -pos 0 " );
  system( "clean_fasta_headers.01.pl -if bacillus.ids.fna -of bacillus.cds.fna  -pos 0 " );
}

exit;

#
# END
#

########################

