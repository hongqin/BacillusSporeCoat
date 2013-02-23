use strict; use warnings; 
use lib "/Users/hongqin/lib/perl"; use Util; 

my $HOME = '/Users/hongqin';
my $debug= 0 ;

my $inprofilefile = '_coat.profile.for.truncated.codeml.test.Jan18,2011.tab';
my $homedir = "$HOME/coat.protein07/ortholog.analysis/coatpaml/test.nearlyAll.2013Jan16.Gblocks";

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

if($debug>9) { @bgs = splice(@bgs, 0, 1);  } #if($debug>9) { @bgs = ($bgs[5], $bgs[7], $bgs[11]); }

my $count = 0;
foreach my $mydir ( @bgs ) {
  chdir $homedir; 
  chomp $mydir;
  chdir "wkdir/$mydir";
  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 
 
  ### align CDS
  #my $cmd1 =  "align_cds_by_protein_bioperl.00.pl -i _cds.spec.fna -o tmpcds.phylip -af phylip ";
  my $cmd1 =  "align_cds_by_protein_bioperl.00.pl -i _cds.spec.fna -o tmpcds.aln -af fasta ";
  system( $cmd1 );

  $cmd1 = "Gblocks tmpcds.aln -t=c -e=.gb";  #2013Jan16
  system( $cmd1 ); 

  #clean headers in gb file
  open(IN2a, "<tmpcds.aln.gb"), my @linesGB = <IN2a>, close (IN2a); 
  open(OUT, ">tmpcds.aln.gb.faa");
  foreach my $line (@linesGB) { 
    if ($line =~ /^>/ ) {
       my @els = split( /\//, $line);
       print OUT "$els[0]\n";
    } else { 
   		print OUT $line; 
   	}
  }
  
  close(OUT);

  $cmd1 =  "align_cds_by_protein_bioperl.00.pl -i tmpcds.aln.gb.faa -o tmpcds.phylip2 -af phylip ";
  system( $cmd1 ); 

  ### indicate the alignment is interleaved. (required by codeml)
  open(IN2, "<tmpcds.phylip2"); my @lines = <IN2>;    close (IN2);  #2013 Jan16

  if ( -e "cds2.phylip" ) { 
     system("mv cds2.phylip cds2.phylip.before.Gblock"); 
  }

  open(OUT, ">cds3.phylip");
    my $firstline = shift @lines;
    chomp $firstline;
    print OUT $firstline."  I\n";
    foreach my $line (@lines) { print OUT $line; }
  close(OUT);

  if ( ! (-e $modeldir) ) { mkdir $modeldir; }
  chdir $modeldir;
  system("cp ../cds3.phylip .");
  system("cp ../bacillus.tree1.nwk .");


}

exit;

#
# END
#

########################


