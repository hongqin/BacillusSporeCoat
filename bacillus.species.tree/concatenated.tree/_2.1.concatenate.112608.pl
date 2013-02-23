 
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 1;

my $inclusterfile = 'my.curated.essen.gene.csv.112608b';

open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); 
my @clusterdirs = ();
foreach my $line (@lines) {
  chomp $line;
  $line =~ s/^\s+//o; 
  my ($nu, $bg, @res) = split( /\s+/, $line);
  push (@clusterdirs, $bg);
}

if( $debug > 9) { @clusterdirs = splice(@clusterdirs, 1, 3 ); }

my $count = 0;

my %concatenated = ();

foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  FASTAshort::paddle_fasta_file2("_cds.spec.fna", "_cds.spec.oneline.fna" );
  
  open(IN, "<_cds.spec.oneline.fna");  my @lines = <IN>;   close(IN);

  my %lens = ();
  for(my $i=0; $i<=$#lines; $i=$i+2) {
     my $header = $lines[$i];
     my $seq = $lines[$i+1];
     $header =~ s/^>\s*//o;
     my ($id, @res) = split( /\s+/, $header);
     if ($debug>9) { print "\$id = [$id]\n"; }
     $lens{ $id } = length( $seq );
     
     chomp $seq; 
     $concatenated{$id}->{SEQ} .= $seq . ' '; #add a space
     $concatenated{$id}->{COUNTS} += 1 . ' '; #add a space
  }
  
  if($debug) { showHashTable(\%lens); } 
  
  chdir "..";
}

#output and check the quality concanted seq
open(OUT,">_concatenated.cds.112608.fna");
foreach my $id ( sort( keys %concatenated ) ) { 
 my $len =  length( $concatenated{$id}->{SEQ} );
 print "$id $len " . $concatenated{$id}->{COUNTS} . "\n";
 print OUT ">$id\n";
 print OUT $concatenated{$id}->{SEQ} . "\n";
}
close(OUT);


exit;

#
# END
#

