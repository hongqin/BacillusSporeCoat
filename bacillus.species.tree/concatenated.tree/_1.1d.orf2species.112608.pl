#convert orf name to 3-letter species name 
 
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 2 ;

my $inclusterfile = 'my.curated.essen.gene.csv.112608b';
#my $inclusterfile = '_test.csv';
my $idfl   = "/home/hqin/projects/coat.protein07/key.data/_id2genomes.csv";
my $specfl = "/home/hqin/projects/coat.protein07/key.data/_genome2species.csv";

open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); 
my @clusterdirs = ();
foreach my $line (@lines) {
  chomp $line;
  $line =~ s/^\s+//o; 
  my ($nu, $bg, @res) = split( /\s+/, $line);
  push (@clusterdirs, $bg);
}

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

$count = 0;

foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  open (IN,  "<cds.fna");
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

  chdir "..";
}

exit;

#
# END
#

