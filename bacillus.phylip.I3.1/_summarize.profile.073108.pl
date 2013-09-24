# 073108 pick entries in the Bacillus clade,trim IDs to 10 letters for PHYLIP
# generate id converstion tables, 
 
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; 
use lib '/Users/hongqin/lib/perl'; 
use Util; use FASTAshort; 

my $debug= 1 ;

my $idfl   = "/Users/hongqin/projects/coat.protein07/key.data/_id2genomes.csv";
my $specfl = "/Users/hongqin/projects/coat.protein07/key.data/_genome2species.csv";
#my $clusterfl = "infile-BGclusters.txt";
#open (IN, "<$clusterfl" ); my @clusterdirs = <IN>; close (IN);

my @orders  = ('Bsu','Bmo','Bam','Bli','Bpu','Ban','Bce','Bth','Bwe','Bcl','Bha');

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
  my ($genome, $spec, $flag, @res) = split (/\t/, $line);
  print "[$genome] [$spec] [$flag]\t";
  $genome2species{ $genome } = { 
   'SPEC' =>  $spec,
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

open(PROFILE, ">_profile.073108b.csv");

#first row, header
print PROFILE "id\tsource\t"; #032410 change
foreach my $spec (@orders) { print PROFILE "$spec\t"; } print PROFILE "Profile11\n"; 

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

  my %_cur_profiles = ();

  open (IN, "<$myfl") ; 
  while (my $line = <IN> ) {
    chomp $line; 
    if ($line !~ /^\s*$/ ) {
     my $id = $line; 
     my $myflag = $genome2species{ $id2genomes{$id} }->{FLAG};
     if ($myflag eq 'Y' ) {
        $_cur_profiles{$genome2species{ $id2genomes{$id} }->{SPEC}} .= $id.' ';
     }
    }
  }
  if ($debug) { showHashTable( \%_cur_profiles); }
  close (IN);

  print PROFILE "$mydir\t$myfl\t";

  my $profile_pattern ='Y';
  foreach my $spec (@orders) { 
    if (exists $_cur_profiles{$spec} ) {
  	print PROFILE $_cur_profiles{$spec} . "\t"; 
    } else {
        print PROFILE "NA\t";
        $profile_pattern='';
    } 
  } 
  print PROFILE "$profile_pattern\n"; 
   
  chdir "..";
}
close(PROFILE);

exit;

#
# END
#

