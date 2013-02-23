#BEGIN { unshift(@INC,"/opt/local/lib/perl5/site_perl/5.8.9"); }

use strict; use warnings; 
use lib "/Users/hongqin/lib/perl"; use Util; 
use Bio::TreeIO;
use Bio::Tools::Phylo::PAML::Result;
use Bio::Tools::Phylo::PAML;

my $HOME = '/Users/hongqin';
my $debug= 1 ;

my $modeldir = "H1C"; #This is the model, but omega is free changing 
my $codemlReportFile = "_free.omega.GBlocked.coat.paml.Jan17,2013.txt"; 

my $inprofilefile = '_coat.profile.for.truncated.codeml.test.Jan18,2011.tab';
#my $homedir = "$HOME/coat.protein07/ortholog.analysis/coatpaml/test.nearlyAll.Jan18,2011";
my $homedir = "$HOME/coat.protein07/ortholog.analysis/coatpaml/test.nearlyAll.2013Jan16.Gblocks";

open (IN,"<$inprofilefile"); my @lines = <IN>; close (IN); shift @lines;
my @bgs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($bg, $num, $paml, @res) = split( /\t/, $line);
  if ($paml =~ m/\s*YES\s*/ ) {
     push (@bgs, $bg);
  }
}

if($debug>9) { @bgs = splice(@bgs, 0, 1);  } #if($debug>9) { @bgs = ($bgs[5], $bgs[7], $bgs[11]); }

open (REPORT, ">$codemlReportFile"); 

my $count = 0;
foreach my $mydir ( @bgs ) {
  chdir $homedir; 
  chomp $mydir;
  chdir "wkdir/$mydir";
  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 

  # push (@outlines, '#'.$mydir);

  #foreach my $i ( 0 .. $#Hs) {
   my $h = $modeldir; 
   chdir $h;
   system("pwd; ls results*");
   open (RES, "<results.$h.txt"); my @results = <RES>; close (RES);

   my @sublines = splice( @results, $#results -3, 1 );
   if ($debug) { print "<@sublines>"; }
   open (omegaTREE, ">_$mydir.$h.out.nwk");
   print omegaTREE $sublines[0]; 
   close (omegaTREE); 
   
   #parse trees
   #my $treeio = Bio::TreeIO->new(-format => 'newick', -file =>'_$mydir.$h.out.nwk'); #now working   
  
   print REPORT "#Group\t$mydir######\n";
    my $paml_parser = Bio::Tools::Phylo::PAML->new(-file => "results.$h.txt", -dir => "./");
  
  if( my $result = $paml_parser->next_result() ) {
   while ( my $tree = $result->next_tree ) {
    for my $node ( $tree->get_nodes ) {
     my $id;
     # first we do some work to figure out what the ID should be.
     # for a leaf or tip node this is just the taxon label
     if( $node->is_Leaf() ) {
        $id = $node->id;
     } else {
       # for the internal nodes it is just the name of all the sub-nodes
       # put together, much like how Sanderson represents internal nodes
       # in r8s
        $id = "(".join(",", map { $_->id } grep { $_->is_Leaf }
                                  $node->get_all_Descendents) .")";
     }
     if( ! $node->ancestor || ! $node->has_tag('t') ) {
       # skip when no values have been associated with this node
       # (like the root node)
       next;
     }
     printf REPORT "$mydir\t%s\tt\t%.3f\tS\t%.1f\tN\t%.1f\tdN/dS\t%.4f\tdN\t%.4f\t".
             "dS\t%.4f\tS*dS\t%.1f\tN*dN\t%.1f\n",
      $id,map { ($node->get_tag_values($_))[0] }
      qw(t S N dN/dS dN dS), 'S*dS', 'N*dN';
    }
   }
  }
}
close (REPORT); 

exit;

#
# END
#

