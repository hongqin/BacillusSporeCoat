# summarize hits for each coat
#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

my $debug = 0;

my $genomedir = '/home/hqin/projects/coat.protein07/genomes';
my $rootdir = '/home/hqin/coat.protein07/orthog.analysis/phylo012708';
my $clusterdir = '/home/hqin/projects/coat.protein07/orthog.analysis/all.vs.all/simple.hits.mcl/single.cluster';

my $outfile = "_rbHits.coat.012708.csv";

### get the coat ids
system("ls -d BG* > /tmp/_ls-file.txt");
my %coats = ();
open (IN, "</tmp/_ls-file.txt");
while (my $line = <IN> ) {
   chomp $line; 
   my($bg, @els) = split (/-/, $line); 
   $coats{$bg} = $bg;
}
close(IN);
if ($debug) { showHashTable( \%coats); }

open (OUT, ">$outfile");
my $count = 0; 
foreach my $bg ( sort ( keys %coats ) ) {
  chdir $rootdir;
  chdir $bg; 
  if ($debug) {    system( "pwd" );  }

  #get the headers
  my $bgfna = "$bg.original.fna";
  system( "grep \">\" $bgfna > /tmp/_headers.txt");
  open(IN, "</tmp/_headers.txt");   my @lines = <IN>;  close(IN);

  #parse the headers
  my %hits = ();
  foreach my $line (@lines) {
      $line =~ s/^>\s*//o;
      chomp $line; 
      my( $id, @res ) = split (/\s+/, $line);
      my $left  = index( $line, '{', 1);
      my $right = index( $line, '}', 1); 
      my $species = substr( $line, ($left+1), ($right-$left-1) );
      if($debug>4){print "* $id $left $right $species\n";}
      if( $id =~ /^BG/ ) { $species = 'BG168' ;  
      } elsif ( $id =~ /^Bmo/ ) { $species = 'Bmo'; }
      $hits{ $id } = $species; 
  }
  if($debug>4) {showHashTable( \%hits); }

  write_hash2file(\%hits, "$bg-ortho.csv", "\t");
  system("ln -s $bg-ortho.csv ortho.csv");

  $count ++;
  print "**** count = $count $bg \n\n";
  if (($debug > 2)&&($count > 2)) {  exit; }
}

close(OUT);

exit;

