use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; 

my $debug= 1 ;

my $inclusterfile = '_selected.BG.for.paml.121708.tab';
my $root = '/home/hqin/coat.protein07/orthog.analysis/coat.codeml121708/';
my $wkdir = "4sub4cer121708";  

open (IN,"<$inclusterfile"); my @clusterdirs = <IN>; close (IN); 

my @Hs = ('H0', 'H1a', 'H1bwe','H1bce', 'H2a','H2b','H2c','H2d','H2e');

my @Htres = ( 
#"( Bpu,(Bli,(Bam,Bsu)),(Bwe,(Bce,(Ban,Bth))));",  #H0
#"( Bpu,(Bli,(Bam,Bsu)),(Bwe #1,(Bce #1,(Ban #1,Bth #1)#1)#1));", #H1a
#"( Bpu,(Bli,(Bam,Bsu)),(Bwe #1,(Bce,(Ban,Bth))) );", #H1bwe
#"( Bpu,(Bli,(Bam,Bsu)),(Bwe #2,(Bce #1,(Ban #1,Bth #1)#1)#1));", #H2a
#"( Bpu,(Bli,(Bam,Bsu)),(Bwe #2,(Bce #2,(Ban #1,Bth #1)#1)#1));", #H2b
#"( Bpu,(Bli,(Bam,Bsu)),(Bwe #2,(Bce #2,(Ban #1,Bth #1)#1)#2));", #H2c
"( Bpu,(Bli,(Bam,Bsu)),(Bwe #1,(Bce #2,(Ban #1,Bth #1)#1)#1));", #H2d
);

if (($debug>9) ) { @clusterdirs = splice( @clusterdirs, 0, 5); }

my @outlines = ();

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $root;     chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 
  chdir $wkdir;

  push (@outlines, '#'.$mydir);

  foreach my $i ( 0 .. $#Hs) {
   my $h = $Hs[$i];
   # if ( ! (-e $h) ) { mkdir $h; }
   chdir $h;
   system("pwd; ls results*");
   open (RES, "<results.4sub4cereus.$h.txt") || die "Cannot open results\n";  ### Change here
   my @results = <RES>; 
   close (RES);

   #push (@outlines, "  $h");

   if ( $h eq 'H0' ) {
     my @sublines = splice( @results, $#results -2, 2 );
     my ($tmp, $dN) = split( /:\s+/, $sublines[0]);
     my ($tmp2, $dS) = split( /:\s+/, $sublines[1]);
     #print @sublines."\n"; print "\$dN=$dN  \$dS=$dS\n";
     push( @outlines, "$mydir\t$h\tAll\t".$dN/$dS."\t\t\t\t"); 
   } elsif ( $h =~ /H1/ ) {
    my @sublines = splice( @results, $#results -2, 2 );
     my @els = split( /Bsu[#\d|\s]+#/, $sublines[1]);
     if(defined $els[1]) {
        my ($bsuomega, @res) = split (/\s+/, $els[1]);
        my @els2 = split (/Bwe[#\d|\s]+#/, $sublines[1]);
        my ($bweomega, @res2) = split (/\s+/, $els2[1]);
        my @els3 = split (/Bce[#\d|\s]+#/, $sublines[1]);
        my ($bceomega, @res3) = split (/\s+/, $els3[1]);
        push (@outlines, "$mydir\t$h\tBsu\t$bsuomega\tBwe\t$bweomega\tBce\t$bceomega");
     } else {
        push (@outlines, "$mydir\t$h\tBsu\tNA\tBwe\tNA\tBce\tNA");
     }
   } elsif ( $h =~ /H2[abcde]/) { #H2s
#"( Bpu,(Bli,(Bam,Bsu)),(Bwe #2,(Bce #1,(Ban #1,Bth #1)#1)#1));", #H2a
#"( Bpu,(Bli,(Bam,Bsu)),(Bwe #2,(Bce #2,(Ban #1,Bth #1)#1)#1));", #H2b
#"( Bpu,(Bli,(Bam,Bsu)),(Bwe #2,(Bce #2,(Ban #1,Bth #1)#1)#2));", #H2c
#"( Bpu,(Bli,(Bam,Bsu)),(Bwe #1,(Bce #2,(Ban #1,Bth #1)#1)#1));", #H2d
     my @sublines = splice( @results, $#results -2, 2 );
     my @els = split( /Bsu #/, $sublines[1]);
     my ($bsuomega, @res) = split (/\s+/, $els[1]);
     my @els2 = split (/Bwe#. #/, $sublines[1]);
     my ($bweomega, @res2) = split (/\s+/, $els2[1]);
     my @els3 = split (/Bce#. #/, $sublines[1]);
     my ($bceomega, @res3) = split (/\s+/, $els3[1]);
     my @els4 = split (/Ban#. #/, $sublines[1]);
     my ($banomega, @res4) = split (/\s+/, $els4[1]);
     push (@outlines, "$mydir\t$h\tBsu\t$bsuomega\tBwe\t$bweomega\tBce\t$bceomega\tBan\t$banomega");
   } else { #error
     print "\nWARNING $mydir $h\n";
   } 

   chdir "..";
  }

}

#print "\n--------\n"; print @outlines; #"\n";

chdir $root; 
open( OUT, ">__omega.coat.4sub4cereus.010609.csv");
foreach my $line (@outlines) {
 print OUT $line."\n";
}
close (OUT);

exit;

#
# END
#

