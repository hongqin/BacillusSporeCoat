use strict; use warnings; 
use lib '/Users/hongqin/lib/perl'; use Util; use FASTAshort; 
my $debug= 1 ;
use Bio::TreeIO;

open(IN,"<truncated.nwk"); my @lines = <IN>; close (IN);
open(OUT, ">truncated2.nwk");
foreach my $line (@lines) {
  $line =~ s/[\;|\s*]\s+\n/\;\n/go; #remove extra space, ensure semicolon
  my $rightnum = $line =~ tr/\)/\)/; 
  my $leftnum  = $line =~ tr/\(/\(/;
  if ($rightnum == $leftnum) { 
     print OUT $line; 
   } else {
     print "\nuneven tree::$line \n";
   }
}
close(OUT);

my $treeio = Bio::TreeIO->new(-format => 'newick', -file =>'truncated2.nwk');
#my $tree = $treeio->next_tree;

my $out = Bio::TreeIO->new(-file => '>_test', -format => 'tabtree');

open(OUT,">_tmp.out.txt");         
my $count = 0; 
while( my $tree = $treeio->next_tree ) {
   $out->write_tree($tree);
   my @nodes = $tree->get_nodes(-order => 'b|breadth' );
   my @leaves = $tree->get_leaf_nodes; 
   my $number_nodes = $tree->number_nodes; 
 
  foreach my $node ( @nodes) {
     #print OUT $node->id.', ';      #.$node->to_string()." ,"; 
   } 


   print OUT "Tree $count:";
   foreach my $ln (@leaves){
     my @lnids;
     push (@lnids, $ln->id); 
     # @lnids = sort @lnids; 
     @lnids = order_big2small_array( @lnids);
     print OUT "@lnids ";
   }
   print OUT "\n$lines[$count]----------\n"; 
   $count = $count + 1; 
}	 
close(OUT);

exit;

my @Htres = ( 
#"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,Ban,Bth)));", #H0
#"(Bha,(Bpu,(Bli,(Bam ,(Bmo #1,Bsu #1)))),(Bwe,(Bce,Ban,Bth)));", #H1a
#"(Bha,(Bpu,(Bli,(Bam #1,(Bmo #1,Bsu #1)#1))),(Bwe,(Bce,Ban,Bth)));", #H1b
#"(Bha,(Bpu,(Bli #1,(Bam #1,(Bmo #1,Bsu #1)#1)#1)),(Bwe,(Bce,Ban,Bth)));", #H1c
#"(Bha,(Bpu #1,(Bli #1,(Bam #1,(Bmo #1,Bsu #1)))),(Bwe,(Bce,Ban,Bth)));", #H1d
#"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,Ban #1,Bth)));", #H2a
#"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce #1,Ban #1,Bth #1)));", #H2b
#"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #1,(Bce #1,Ban #1,Bth #1)#1));", #H2c
#"(Bha,(Bpu,(Bli #2,(Bam #2,(Bmo #2,Bsu #2)#2)#2)),(Bwe #1,(Bce #1,Ban #1,Bth #1)#1));", #H3a
#"(Bha,(Bpu #2,(Bli #2,(Bam #2,(Bmo #2,Bsu #2)#2)#2)#2),(Bwe #1,(Bce #1,Ban #1,Bth #1)#1));" #H3b
 "(Bha,(Bpu ,(Bli,(Bam ,(Bmo ,Bsu )) ) ),(Bwe #2,(Bce #2,Ban #1,Bth #2)#2) );", #H3c
 "(Bha,(Bpu ,(Bli#3, (Bam #3,(Bmo #3 ,Bsu #3)#3) #3) ),(Bwe #2,(Bce #2,Ban #1,Bth #2)#2) );" #H4a
);

