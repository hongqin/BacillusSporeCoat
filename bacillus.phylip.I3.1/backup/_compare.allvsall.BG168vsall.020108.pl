#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

# 020108 compare orthologous hits from BG168.vs.all all.vs.all approaches

my $debug = 5;
 
my $datadir = '/home/hqin/projects/coat.protein07/key.data';
my $bghitfile  = "$datadir/_coat-phylo-profiles.122107.csv";
my $allhitfile = "$datadir/_rbHits.coat.012808.csv";

my $rootdir = '/home/hqin/coat.protein07/orthog.analysis/phylo012708/mysql';

#####

chdir $rootdir;

# the BG168 vs all hits
open(IN, "<$bghitfile"); my @bglines = <IN>; close(IN);
my $firstline = shift (@bglines); chomp $firstline; 
my @species = split( /\t/, $firstline);
if( $debug ) { print 'There are '.(scalar @species) . " species:[@species]\n"; }

#the allvsall hits
open(INALL, "<$allhitfile"); my @all_lines =<INALL>; close (INALL);


# now compare them
open( OUT1, ">_bghits.flaged.by.allhits.020108.csv");
open( OUT2, ">_allhits.flaged.by.bghits.020108.csv");
my $count = 0;
foreach my $bgline (@bglines ) { #loop over every entry for bg168vsall 
    ++$count; 
    chomp $bgline; 
    my ($bg1,@hits) = split( /\t+/,$bgline);
    
    #parse the allvsall hits for $bg1
    my @sub_alllines = ();
    my %sub_allhits = ();
    foreach my $allhit (@all_lines) {
       chomp $allhit; 
       my ($name2,$bg2,$hit2,$spcies2, @res2) = split( /\t/, $allhit );
       if( $bg1 eq $bg2 ) {
           push ( @sub_alllines, $allhit );
           $sub_allhits{$hit2} = $allhit;
       }
    }
    
    if($debug) { print "\n\n$count bgline=[$bgline]\n";
    	print "sub_allhits=[". ( keys %sub_allhits) . "]\n";
    }
    
    #now compare @hits to %sub_allhits
    my $subcount = 1;
    my %simplehits = ();
    if ( (scalar @hits) > 0 ) {
    foreach my $hit1 (@hits) {
        if (exists $sub_allhits{$hit1} ) {
            print OUT1 "$bg1\t$hit1\tMCLPassed\n";
        } else {
            print OUT1 "$bg1\t$hit1\tMCLfailed\n";
            print "warning simplehit $bg1 $hit1 is not in allhits, freq ". $subcount ++ ."\n";
        } 
        $simplehits{$hit1} = '';
    }
    } else { #no simple hit
       print "*************************";
       showHashTable(\%sub_allhits);
       print "*************************";
       if ((scalar @sub_alllines) == 1) {
          my ($name2,$bg2,$hit2,$spcies2, @res2) = split( /\t/, $sub_alllines[0] );
          if ($hit2 eq 'NA') {
            print OUT1 "$bg1\tNA\tMCLPassed\n";
          } else {
            print OUT1 "$bg1\tNA\tMCLFailed\n";
          } 
       } else {
          print OUT1 "$bg1\tNA\tMCLfailed\n";
       } 
    }

    #now compare %sub_allhits to %simplehits;
    my $num = 0;
    foreach my $hit2 ( keys %sub_allhits ) {
        if (exists $simplehits{$hit2}) {
            print OUT2 "$bg1\t$hit2\tSimpleHitPassed\n";
        }elsif ($hit2 ne 'NA') {
            print OUT2 "$bg1\t$hit2\tSimpleHitFailed\n";
            print "warning allhit $bg1 $hit2 is not in simplehits, freq ". ++ $num . "\n";
        } else {
            print OUT2 "$bg1\t$hit2\tSimpleHitPassed\n";
        } 
    }
}

close(OUT1); close(OUT2);

exit;
