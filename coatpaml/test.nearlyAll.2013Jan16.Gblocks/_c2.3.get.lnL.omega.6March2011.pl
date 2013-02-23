use strict; use warnings; 
BEGIN {    push(@INC,"/home/hqin/lib/perl/", "/Users/hongqin/lib/perl");};
#use lib $ENV{'QINPLIB'};  #use lib "/Users/hongqin/lib/perl"; 
use Util; 
use Bio::TreeIO;
use Bio::Tools::Phylo::PAML::Result;
use Bio::Tools::Phylo::PAML;

my $HOME = $ENV{'HOME'}; 
my $debug= 9 ;

my @HXs = qw(H0 H2C1S1); 
my @codemlmodels = qw(0 2); 
my $omegaReportFile = "_coat.omega.6March2011.tab"; 
my $lnLReportFile = "_coat.lnL.6March2011.tab"; 

#my $pwd = "$HOME/coat.protein07/ortholog.analysis/nonCE.codeml.21Jan2011";
my $pwd = "$HOME/coat.protein07/ortholog.analysis/coatpaml/test.nearlyAll.Jan18,2011";

my $tmpfile = "/tmp/___tmp.coat.codeml.txt". get_short_time_stamp_US(); 
system("ls -d wkdir/BG* > $tmpfile");
open (IN,"<$tmpfile"); my @lines = <IN>; close (IN); 
my @bgs = ();
foreach my $line (@lines) {  chomp $line;  my ($tmp, $bg, @res) = split( /\//, $line);  push (@bgs, $bg);  }

if($debug>9) { @bgs = splice(@bgs, 1, 50); }
#if($debug>9) { @bgs = splice(@bgs, 0, 500); }

open (LNL, ">$lnLReportFile"); 
open (OMEGA, ">$omegaReportFile"); 

my $count = 0;
foreach my $bg ( @bgs ) {
  $count ++; print "Cluster:$count I am now working on [$bg]::\n"; 
 
  foreach my $i ( 0 .. $#HXs ) {
   my $h = $HXs[$i];  chdir $pwd."/wkdir/$bg/$h"; 

   my $tree = _get_tree("bacillus.tree$h.nwk"); 
   print LNL "#Group\t$bg######$tree\n";
   print OMEGA "#Group\t$bg######\n";
   my $paml_parser = Bio::Tools::Phylo::PAML->new(-file => "results.$h.txt", -dir => "./");

  if( my $result = $paml_parser->next_result() ) {
   print LNL "$bg\t$h\t";
   if ($h eq "H0") {
     print LNL "df\t0\t"; 
   } else {
     my $df = _determine_df("bacillus.tree$h.nwk"); 
     print LNL "df\t$df\t"; 
   }
   my $lnL = _get_lnL("results.$h.txt");    
   print LNL "lnL\t", $lnL, "\n";

   while ( my $tree = $result->next_tree ) {
    for my $node ( $tree->get_nodes ) {
     my $id;
     # first we do some work to figure out what the ID should be. for a leaf or tip node this is just the taxon label
     if( $node->is_Leaf() ) {
        $id = $node->id;
     } else {
     # for the internal nodes it is just the name of all the sub-nodes put together, much like how Sanderson represents internal nodes in r8s
        $id = "(".join(",", map { $_->id } grep { $_->is_Leaf }
                                  $node->get_all_Descendents) .")";
     }
     if( ! $node->ancestor || ! $node->has_tag('t') ) {
       # skip when no values have been associated with this node (like the root node)
       next;
     }
     printf OMEGA "$bg\t$h\t%s\tt\t%.3f\tS\t%.1f\tN\t%.1f\tdN/dS\t%.4f\tdN\t%.4f\t".
             "dS\t%.4f\tS*dS\t%.1f\tN*dN\t%.1f\n",
      $id,map { ($node->get_tag_values($_))[0] }
      qw(t S N dN/dS dN dS), 'S*dS', 'N*dN';
    }#for $node loop
   }#while loop
  }#if loop
 }#HXs loop
}#for $bg loop

close (OMEGA); close (LNL); 

exit; 

##
## END
##

####################
#  my $lnL = _get_lnL("results.$h.txt");    
 sub _get_lnL {
   my ($fl) = @_;
   my $debug = 0; 
   open (IN2, "<$fl"); my @lines2=<IN2>; close (IN2);   
   my @tmplines = grep (/^lnL/, @lines2);
   my ($tmp, $second, $third, $fourth) = split( /:/, $tmplines[0] ); 
   my @els = split (/\s+/, $fourth); 
   if ($debug) {
     print $tmplines[0];
     for (my $i=0; $i<=$#els; $i++) {
       print "[$i:][$els[$i]]"; 
     }
     print "\n"; 
   }
   return $els[1]; 
 }

####################
#bacillus.treeH2C1S1.nwk
#     my $df = _determine_df("bacillus.tree$h.nwk"); 
 sub _determine_df {
   my ($treefl) = @_; 
   my $debug = 1; 
   open (IN, "<$treefl"); my @lines = <IN>; close (IN);
   my $tree = $lines[0]; 

   if (( $tree =~ /\#2/ )| ($tree =~ /\$2/)) { return 2; 
   } elsif (( $tree =~ /\#1/ ) | ($tree =~/\$1/) ){ return 1; 
   } else { return 0; }			  
 }


############
#     my $tree = _get_tree("bacillus.tree$h.nwk"); 
 sub _get_tree {
   my ($treefl) = @_; 
   my $debug = 0; 
   open (IN, "<$treefl"); my @lines = <IN>; close (IN);
   my $tree = $lines[0]; 
   chomp $tree; 
   return $tree; 
 }
