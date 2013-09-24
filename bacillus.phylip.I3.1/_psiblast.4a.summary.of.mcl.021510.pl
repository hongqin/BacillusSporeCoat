#v1, Feb 15, 2010, profile based mcl-psiblast results
# 

use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 10 ;

my $coatfl = "/home/hqin/projects/coat.protein07/key.data/BGcoat.tab";
open (IN, "<$coatfl" ); my @coatBGs = <IN>; close (IN);

my $idfl   = "/home/hqin/projects/coat.protein07/key.data/_id2genomes.csv";
my $specfl = "/home/hqin/projects/coat.protein07/key.data/_genome2speciesv2.csv";

my @orders  = ('Bsu','Bmo','Bam','Bli','Bpu','Ban','Bce87', 'Bce79', 'Bce3L','Bth','Bwe','Bcl','Bha');

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
  #print "[$genome] [$spec] [$flag] [$spec2] [$flag2]\t";
  print "[$genome] [$spec2] [$flag2]\t";
  $genome2species{ $genome } = { 
   'SPEC' =>  $spec2,
   'FLAG' =>  $flag2 };
  } else { $count ++;
  } 
}
close (IN);

print "\n*************************\n";
foreach my $key (sort( keys %genome2species ) ) {
  print $genome2species{$key}->{SPEC} ."::".$genome2species{$key}->{FLAG} . "\n";
}
print "\n*************************\n";

#my %RUNs = ();
$count = 0;
my %profile = ();

#my @orders  = ('Bsu','Bmo','Bam','Bli','Bpu','Ban','Bce87', 'Bce79', 'Bce3L','Bth','Bwe','Bcl','Bha');
open (PROFILE, ">coat.psi-mcl.run/_profile.from.psi.mcl.I31.clus.021510.csv");
print PROFILE "coat"; foreach my $spec (@orders) { print PROFILE "\t$spec"; } 
print PROFILE "\tNumOfHits\tNumOfGenomeHits\n";

open (IN, "<coat.psi-mcl.run/_coat.psi.mcl.I31.clus"); my @mcllines = <IN>; close (IN); 

my $count = 0; 
foreach my $myBG ( sort @coatBGs ) {
  chomp $myBG; 
  $count ++;
  print "***** $count $myBG ******\n";
  my @hitlines = grep( /$myBG/, @mcllines );
  my @hits = split ( /\s+/, $hitlines[0] );

  my %_cur_profiles = ();
  foreach my $hit (@hits ) {
     my $myflag = $genome2species{ $id2genomes{$hit} }->{FLAG};
     if ($myflag eq 'Y' ) {
        $_cur_profiles{$genome2species{ $id2genomes{$hit} }->{SPEC}} .= $hit.' ';
     }
  }
  if ($debug > 5) { 
	print "##\%_cur_profiles::\n";
	showHashTable( \%_cur_profiles); 
  }

  print PROFILE "$myBG";
  foreach my $spec (@orders) { 
    if (exists $_cur_profiles{$spec} ) {
  	print PROFILE "\t".$_cur_profiles{$spec}; 
    } else {
        print PROFILE "\tNA";
    } 
  } 
  print PROFILE "\t". ($#hits + 1);
  my @genomehits = keys( %_cur_profiles );
  print PROFILE "\t". ($#genomehits + 1) ."\n"; 
  ### END, print PROFILE for $myBD;

}

close (PROFILE); 

exit;

#
# END
#

