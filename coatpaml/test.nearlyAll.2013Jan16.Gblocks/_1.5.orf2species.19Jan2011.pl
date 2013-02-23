#convert orf name to 3-letter species name 
use strict; use warnings; 
use lib "/Users/hongqin/lib/perl"; use Util; 

my $HOME = '/Users/hongqin';
my $debug= 1 ;

my $inprofilefile = '_coat.profile.for.truncated.codeml.test.Jan18,2011.tab';
#my $sourcedir = "$HOME/coat.protein07/orthog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1/";
my $homedir = "$HOME/coat.protein07/ortholog.analysis/coatpaml/test.nearlyAll.Jan18,2011";

my $idfl   = "$HOME/projects/coat.protein07/key.data/_id2genomes.csv";
my $specfl = "$HOME/projects/coat.protein07/key.data/_genome2species.csv";

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



my %orf2genomes = ();
open (IN, "<$idfl" );
while (my $line = <IN> ) {
  chomp $line;
  my ($orf, $g, @res) = split (/\t/, $line);
  $orf2genomes{$orf} = $g;
}
close (IN);

my %genome2species = ();
open (IN, "<$specfl");
my $count = 1;
while (my $line = <IN> ) {
  if ( $count >= 2 ) {  #skip the header
  chomp $line;
  my ($genome, $spec, $flag, @res) = split (/\t/, $line);
  print "[$genome] [$spec] [$flag]\t";
  $genome2species{ $genome } = { 
   'SPEC' => $spec,
   'FLAG' =>  $flag };
  } else { $count ++;
  } 
}
close (IN);

print "\n*************************\n";
foreach my $key (keys %genome2species ) {
  print $genome2species{$key}->{SPEC} ."::".$genome2species{$key}->{FLAG} . "\t\t";
}
print "\n*************************\n";

#my $count = 0;
$count = 0; 
foreach my $mydir ( @bgs ) {
  chdir $homedir; 
  chomp $mydir;
  chdir "wkdir/$mydir";
  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 

  open (IN,  "<bacillus.cds.fna");
  open (OUT, ">_cds.spec.fna" );
  open (TAB, ">_orf.spec.csv");
  while (my $line = <IN> ) {
    if( $line =~ /^>/ ) {
   	chomp $line;
   	$line =~ s/^>//o;
   	$line =~ s/\s+//g;
   	print OUT ">". $genome2species{ $orf2genomes{$line} }->{SPEC} . "\n";
   	print TAB "$line\t". $genome2species{ $orf2genomes{$line} }->{SPEC} . "\n";
    }else {
    	print OUT $line; 
    }
  }
  close (IN);
  close (OUT);
  close (TAB);

}

exit;

#
# END
#

########################
