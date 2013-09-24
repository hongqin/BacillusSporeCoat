#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

#summarize for mysql 

my $debug = 5;

my $genomedir = '/home/hqin/projects/coat.protein07/genomes';
my $rootdir = '/home/hqin/coat.protein07/orthog.analysis/phylo012708';
my $clusterdir = '/home/hqin/projects/coat.protein07/orthog.analysis/all.vs.all/simple.hits.mcl/single.cluster';

my $outfile = "_rbHits.coat.020108.csv";

#my %coats = ();
my %coat2name = ();

### get the coat ids
open (IN, "</home/hqin/projects/coat.protein07/key.data/BGcoat.csv");
while (my $line = <IN> ) {
   if ($line !~ /NCBI/ ) {
	chomp $line; 
	my($name, $bg, @els) = split (/\t/, $line); 
	#$coats{$bg} = '';
	$coat2name{$bg} = $name; 
   }
}
close(IN);
if ($debug) { showHashTable( \%coat2name); 
 my @els = sort (values %coat2name) ;
 print "Values are [@els]\n";
}

my %species = (); 
#my %id2species = (); 

chdir $rootdir;
open (OUT, ">mysql/$outfile");
foreach my $bg ( sort ( keys %coat2name) ) {
  #parse the hits for each coat prot, store in %coats; %species; 
  #chdir $rootdir;  chdir $bg; 
  open(IN, "<$rootdir/$bg/$bg-ortho.csv");   my @lines = <IN>;    close (IN);
  
  if (scalar @lines > 0 ) {
	foreach my $line (@lines) {
		chomp $line;
		my ($id, @els) =  split( /\s+/, $line );
		my $spec = join (' ', @els);
		#$coats{ $bg }   =  $coats{ $bg } . ' '. $id; 
		$species{ $spec } += 1; 
		#$id2species{ $id } = $spec; 
		print OUT "$coat2name{$bg}\t$bg\t$id\t$spec\n";
	}
  } else {
  	print OUT print OUT "$coat2name{$bg}\t$bg\tNA\tNA\n";
  } 
}
close (OUT);

chdir $rootdir;

if($debug) { #showHashTable(\%species);  #showHashTable(\%coats);   #showHashTable(\%id2species);
}

open(OUT2, ">/tmp/_hits.by.species.csv");
#my @sorted_species = ('BG168');
foreach my $spec (sort (keys %species))  {
   print OUT2 "$spec\t$species{$spec}\n";
};
close(OUT2);

 #if($debug>2) { print "\n**[@sorted_species]***\n";   #showHashTable(\%species);   # showHashTable(\%coats);



exit;
