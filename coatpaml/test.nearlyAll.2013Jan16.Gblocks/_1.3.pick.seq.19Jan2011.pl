use strict; use warnings; 
use lib "/Users/hongqin/lib/perl"; use Util; 

my $HOME = '/Users/hongqin';
my $debug= 1 ;

my $inprofilefile = '_coat.profile.for.truncated.codeml.test.Jan18,2011.tab';
#my $sourcedir = "$HOME/coat.protein07/orthog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1/";
my $homedir = "$HOME/coat.protein07/ortholog.analysis/coatpaml/test.nearlyAll.Jan18,2011";
my $allfaa = "$HOME/projects/coat.protein07/genomes/_merged.all.062108.faa";
my $allfna = "$HOME/projects/coat.protein07/genomes/_merged.all.062108.fna";

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

##########
my $count = 0;
foreach my $mydir ( @bgs ) {
  chdir $homedir; 
  chomp $mydir;
  chdir "wkdir/$mydir";
  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 

  my $cmd = "pick_fasta_records_by_ids.1.pl -if $allfaa -of bacillus.ids.faa -id bacillus.ids.tab";
  if ($debug) { print $cmd . "\n"; }
  system ( $cmd );
  my $cmd2 = "pick_fasta_records_by_ids.1.pl -if $allfna -of bacillus.ids.fna -id bacillus.ids.tab";
  if ($debug) { print $cmd2 . "\n"; }
  system ( $cmd2 );
  
  #check parsing consistency
  my $cmd3 = "wc -l bacillus.ids.tab "; print $cmd3."\n";
  system ( $cmd3 );
  my $cmd4 = "grep \"^>\" bacillus.ids.faa | wc -l "; print $cmd4."\n";
  system ( $cmd4 );
  my $cmd5 = "grep \"^>\" bacillus.ids.fna | wc -l "; print $cmd5."\n"; 
  system ( $cmd5 );

 open (IN, "<bacillus.ids.tab"); my @lines = <IN>; close (IN);
  my $num_of_ids = $#lines + 1;
  print "\$num_of_ids=[$num_of_ids]\n";
    
  open (IN, "<bacillus.ids.faa"); @lines = <IN>; close (IN);
  my @lines2 = grep ( /^>/, @lines );
  my $num_of_faaLines = $#lines2 + 1;
  print "\$num_of_faaLines=[$num_of_faaLines]\n";
  
  open (IN, "<bacillus.ids.fna"); @lines = <IN>; close (IN);
  my @lines3 = grep (/^>/, @lines );
  my $num_of_fnaLines = $#lines3 + 1;
  print "\$num_of_fnaLines=[$num_of_fnaLines]\n";  
  
  if ( ( ($num_of_ids - $num_of_faaLines) == 0 ) && ( ($num_of_ids - $num_of_fnaLines) == 0 ) ) {
     system( "touch __yes.faafna.entries.match.id.numbers.062108.txt" );
  } else {
     system( "touch __no.faafna.entries.NOT.match.id.numbers.062108.txt" );
  } 

}

exit;

#
# END
#

########################

