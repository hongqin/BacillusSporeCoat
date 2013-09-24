# 030610 partition psi-MCL clusers by coat BG ids
# output "ids.tab" in each cluster-directory

use strict; use warnings; 
use lib '/home/hqin/lib/perl'; 
use lib '/Users/hongqin/lib/perl';

use Util; 

my $debug= 5 ;
my $coat_fl = "/Users/hongqin/coat.protein07/key.data/BGcoat2.tab";
#my $rbh_fl = "original.mcl.rbh/_simple.coat.orth.pairs.060208.csv";
#my $in_merged_clus_fl = "_rbph.and._coat.I31.clus";
my $in_merged_clus_fl = "_inconsistent.psiprofile.022410.csv";

### read in bg coat ids
my %coat =();

open (IN, "<$coat_fl");
while (my $line = <IN>) {
 chomp $line; 
 my ($name, $bg, @rest) = split( /\t/, $line);
 $coat{ $bg } = $name;
}
close (IN);

if($debug) { showHashTable( \%coat ); }

  open (IN, "<$in_merged_clus_fl");
  my @lines = <IN>;
  close(IN);
  shift (@lines);
  
  my $count = 0;
#  while ( my $line = <IN> ) {
foreach my $line (@lines) {
     chdir "/Users/hongqin/bacillus.phylip.I3.1/coat.psi-mcl.run/phylogeny.psiI31";
     $count ++;
     chomp $line;
     $line =~ s/\"//g;
     my @els = split (/\,|\s+/, $line);
     #@els = sort( @els );
     
     #find out coat BG ids
     my @current_coatBGs = ();
     foreach my $el (@els) {
        #$el =~ s/\"//g; 
        if ( exists $coat{$el} ) {
           push (@current_coatBGs, $el );
           if ($debug>=3) { print "$count \@current_coatBGs are [@current_coatBGs]\n";}
        }
     }
     
     my $mydir = $current_coatBGs[0];
     #if ( $#current_coatBGs >= 1) {
     #   for ( my $i=1; $i <= $#current_coatBGs; $i++) {
     #       $mydir .= '-'.$current_coatBGs[ $i ]; 
     #   }
     #}
     if ($debug>1) { print "\$mydir = [$mydir]\n"}
     
     if ( ! -e $mydir ) { mkdir $mydir ; }
     chdir $mydir; 
     system ( "pwd" );
     
     my %unique_ids = ();
     array2hash( \%unique_ids, \@els, \@els);
     open( OUT, ">ids.tab" );
     foreach my $id (sort (keys %unique_ids) ) {
        if (( $id !~ /^\s*$/ ) & ($id !~/NA/ )) {  #remove empty and NA ones
         print OUT "$id\n";
       }
     }

     chdir "/Users/hongqin/bacillus.phylip.I3.1/coat.psi-mcl.run/phylogeny.psiI31";
     system ( "pwd" );
  }

  close (IN);

exit;
#
# END
#

########################
