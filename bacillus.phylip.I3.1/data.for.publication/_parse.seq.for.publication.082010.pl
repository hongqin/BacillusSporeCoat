#082010 parse seq for publication

# 073108 pick entries in the Bacillus clade,trim IDs to 10 letters for PHYLIP
# generate id converstion tables, 
 
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; 
use lib '/Users/hongqin/lib/perl';
use Util; use FASTAshort; 

my $debug= 0 ;

my $idfl   = "/Users/hongqin/projects/coat.protein07/key.data/_id2genomes.csv";
my $specfl = "/Users/hongqin/projects/coat.protein07/key.data/_genome2speciesv2.csv";
#my $clusterfl = "infile-BGclusters.txt";
#open (IN, "<$clusterfl" ); my @clusterdirs = <IN>; close (IN);

my @orders  = ('Bsu','Bmo','Bam','Bli','Bpu','Ban','Bce79', 'Bce87','Bce3L','Bth','Bwe','Bcl','Bha');

system( "ls -d BG* > /tmp/_tmpbg.txt" );
open(IN, "</tmp/_tmpbg.txt");
my @clusterdirs = ();
while(my $line =<IN>) {
  chomp $line;
  push (@clusterdirs, $line);
}
close(IN);

my %id2genomes = ();
open (IN, "<$idfl" );
while (my $line = <IN> ) {
  chomp $line;
  my ($id, $g, @res) = split (/\t/, $line);
  $id2genomes{$id} = $g;
}
close (IN);

my %genome2species = ();
open (IN, "<$specfl");
my $count = 1;
while (my $line = <IN> ) {
  if ( $count >= 2 ) {  #skip the header
  chomp $line;
  my ($genome, $spec, $flag, $spec2, $flag2, @res) = split (/\t/, $line);
  print "[$genome] [$spec2] [$flag2]\t";
  $genome2species{ $genome } = { 
   'SPEC' =>  $spec2,
   'FLAG' =>  $flag2 };
  } else { $count ++;
  } 
}
close (IN);

print "\n*************************\n";
foreach my $key (keys %genome2species ) {
  print $genome2species{$key}->{SPEC} ."::".$genome2species{$key}->{FLAG} . "\t\t";
}
print "\n*************************\n";


$count = 0;
my %profile = ();
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  if (($debug >=9 )&&($count>5) ){ exit; }

  my $myfl = '';
  if( -e 'ids.tab.curated') { 
     $myfl = 'ids.tab.curated'; 
     if ($debug) { print "$mydir ids.tab.curated\n" } 
  } else {
     $myfl = 'ids.tab';
     #if ($debug) { print "$mydir ids.tab\n" } 
  } 

 system("rm *faa *fna " );

 system("pick_fasta_records_by_ids.1b.pl -if /Users/hongqin/coat.protein07/genomes/_merged.all.062108.faa -id $myfl -o $mydir.prot.082010.faa");
 system("pick_fasta_records_by_ids.1b.pl -if /Users/hongqin/coat.protein07/genomes/_merged.all.062108.fna -id $myfl -o $mydir.dna.082010.fna ");

   
 chdir "..";
}

exit;

#
# END
#

