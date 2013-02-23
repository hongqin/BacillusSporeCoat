#090508 pick BG clusters for paml 
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; 

my $debug= 2 ;

my $inclusterfile = 'my.curated.essen.gene.csv.112608b';
#my $sourcedir = '/home/hqin/projects/coat.protein07/orthog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1/';
my $allfaa = "/home/hqin/projects/coat.protein07/genomes/_merged.all.062108.faa";
my $allfna = "/home/hqin/projects/coat.protein07/genomes/_merged.all.062108.fna";

open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); 
my @clusterdirs = ();
foreach my $line (@lines) {
  chomp $line;
  $line =~ s/^\s+//o; 
  my ($nu, $bg, @res) = split( /\s+/, $line);
  push (@clusterdirs, $bg);
}

if( $debug > 9 ) { @clusterdirs = splice(@clusterdirs, 1, 3); } 

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 

  my $cmd = "pick_fasta_records_by_ids.1b.pl -if $allfaa -of ids.faa -id _ids.tab";
  if ($debug) { print $cmd . "\n"; }
  system ( $cmd );
  my $cmd2 = "pick_fasta_records_by_ids.1b.pl -if $allfna -of ids.fna -id _ids.tab";
  if ($debug) { print $cmd2 . "\n"; }
  system ( $cmd2 );
  
  #check parsing consistency
  my $cmd3 = "wc -l _ids.tab "; print $cmd3."\n";
  system ( $cmd3 );
  my $cmd4 = "grep \"^>\" ids.faa | wc -l "; print $cmd4."\n";
  system ( $cmd4 );
  my $cmd5 = "grep \"^>\" ids.fna | wc -l "; print $cmd5."\n"; 
  system ( $cmd5 );

 open (IN, "<_ids.tab"); my @lines = <IN>; close (IN);
  my $num_of_ids = $#lines + 1;
  print "\$num_of_ids=[$num_of_ids]\n";
    
  open (IN, "<ids.faa"); @lines = <IN>; close (IN);
  my @lines2 = grep ( /^>/, @lines );
  my $num_of_faaLines = $#lines2 + 1;
  print "\$num_of_faaLines=[$num_of_faaLines]\n";
  
  open (IN, "<ids.fna"); @lines = <IN>; close (IN);
  my @lines3 = grep (/^>/, @lines );
  my $num_of_fnaLines = $#lines3 + 1;
  print "\$num_of_fnaLines=[$num_of_fnaLines]\n";  
  
  if ( ( ($num_of_ids - $num_of_faaLines) == 0 ) && ( ($num_of_ids - $num_of_fnaLines) == 0 ) ) {
     system( "touch __yes.faafna.entries.match.id.numbers.112608.txt" );
  } else {
     system( "touch __no.faafna.entries.NOT.match.id.numbers.112608.txt" );
  } 


  chdir "..";
}

exit;

#
# END
#


