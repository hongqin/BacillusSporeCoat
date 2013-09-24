# 030610 partition psi-MCL clusers by coat BG ids
# output "ids.tab" in each cluster-directory

use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; 

my $debug= 2 ;
my $coat_fl = "/home/hqin/coat.protein07/key.data/BGcoat2.tab";
#my $rbh_fl = "original.mcl.rbh/_simple.coat.orth.pairs.060208.csv";
#my $in_merged_clus_fl = "_rbph.and._coat.I31.clus";

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

#if($debug>5) { showHashTable(\%coat_ortho);}

  open (IN, "<$in_merged_clus_fl");
  
  my $count = 0;
  while ( my $line = <IN> ) {
     $count ++;
     chomp $line;
     my @els = split (/\t/, $line);
     
     #find out coat BG ids
     my @current_coatBGs = ();
     foreach my $el (@els) {
        if ( exists $coat{$el} ) {
           push (@current_coatBGs, $el );
           if ($debug>=3) { print "$count \@current_coatBGs are [@current_coatBGs]\n";}
        }
     }
     
     my $mydir = $current_coatBGs[0];
     if ( $#current_coatBGs >= 1) {
        for ( my $i=1; $i <= $#current_coatBGs; $i++) {
            $mydir .= '-'.$current_coatBGs[ $i ]; 
        }
     }
     if ($debug>1) { print "\$mydir = [$mydir]\n"}
     
     if ( ! -e $mydir ) { mkdir $mydir ; }
     chdir $mydir; 
     system ( "pwd" );

     open( OUT, ">ids.tab" );
     foreach my $id (sort @els ) {
        if ( $id !~ /^\s*$/ ) {  #remove empty ones
          print OUT "$id\n";
        }
     }

     chdir "..";
     system ( "pwd" );
  }

  close (IN);

exit;
#
# END
#

########################

