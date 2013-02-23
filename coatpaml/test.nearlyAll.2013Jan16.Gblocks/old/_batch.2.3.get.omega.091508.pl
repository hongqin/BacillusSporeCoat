use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; 

my $debug= 1 ;

my $inclusterfile = '_073108.test.csv';

open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); 
  
shift @lines;

my @clusterdirs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($bg, @res) = split( /\t/, $line);
    push (@clusterdirs, $bg);
}

my $root = "/home/hqin/projects/coat.protein07/orthog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1/paml/test1.5sub4cer/";
my $wkdir = "Bha5sub4cer090108";

#my $long_choices = "Bha Bpu Bli Bam Bmo Bsu Bwe Bce Ban Bth";
#my @choices = split( /\s+/, $long_choices);

my @Hs = ('H0','H1a', 'H1b', 'H1c', 'H1d','H2a','H2b', 'H2c', 'H3a','H3b');
#my @Hs = ('H0');

my @Htres = ( "(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,Ban,Bth)));", #H0
"(Bha,(Bpu,(Bli,(Bam ,(Bmo #1,Bsu #1)))),(Bwe,(Bce,Ban,Bth)));", #H1a
"(Bha,(Bpu,(Bli,(Bam #1,(Bmo #1,Bsu #1)#1))),(Bwe,(Bce,Ban,Bth)));", #H1b
"(Bha,(Bpu,(Bli #1,(Bam #1,(Bmo #1,Bsu #1)#1)#1)),(Bwe,(Bce,Ban,Bth)));", #H1c
"(Bha,(Bpu #1,(Bli #1,(Bam #1,(Bmo #1,Bsu #1)))),(Bwe,(Bce,Ban,Bth)));", #H1d
"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,Ban #1,Bth)));", #H2a
"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce #1,Ban #1,Bth #1)));", #H2b
"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #1,(Bce #1,Ban #1,Bth #1)#1));", #H2c
"(Bha,(Bpu,(Bli #2,(Bam #2,(Bmo #2,Bsu #2)#2)#2)),(Bwe #1,(Bce #1,Ban #1,Bth #1)#1));", #H3a
"(Bha,(Bpu #2,(Bli #2,(Bam #2,(Bmo #2,Bsu #2)#2)#2)#2),(Bwe #1,(Bce #1,Ban #1,Bth #1)#1));" #H3b
);
#"(Bha,((Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #1,(Bce #1,(Ban #1,Bth #1)#1)#1)));",
#"(Bha,((Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,(Ban #1,Bth)))));"

#my @Hctls = ( "$root/model.H0.ctl", "$root/model.H1.ctl", "$root/model.H2.ctl", "$root/model.H3.ctl");
#my $Hxctl = "$root/model.HX.ctl";

my $count = 0;

my @outlines = ();

if (($debug>9) ) { 
@clusterdirs = splice( @clusterdirs, 0, 3);
}

foreach my $mydir ( @clusterdirs ) {
 
  chdir $root;   
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 
  chdir $wkdir;

  push (@outlines, '#'.$mydir);

  foreach my $i ( 0 .. $#Hs) {
   my $h = $Hs[$i];
   # if ( ! (-e $h) ) { mkdir $h; }
   chdir $h;
   #system("pwd; ls results*");
   open (RES, "<results.$h.txt"); my @results = <RES>; close (RES);

   #push (@outlines, "  $h");

   if ( $h eq 'H0' ) {
     my @sublines = splice( @results, $#results -2, 2 );
     my ($tmp, $dN) = split( /:\s+/, $sublines[0]);
     my ($tmp2, $dS) = split( /:\s+/, $sublines[1]);
     #print @sublines."\n"; print "\$dN=$dN  \$dS=$dS\n";
     push( @outlines, "$mydir\t$h\tAll\t".$dN/$dS."\t\t\t\t"); 
   } elsif ( $h =~ /H1/ ) { #two omega
     my @sublines = splice( @results, $#results -2, 2 );
     my @els = split( /Bsu#1 #/, $sublines[1]);
     #print @sublines; print "::::\n"; 
     #print $els[1]; 
     my ($bsuomega, @res) = split (/\s+/, $els[1]);
     my @els2 = split (/Bha #/, $sublines[1]);
     my ($bhaomega, @res2) = split (/\s+/, $els2[1]);
     push (@outlines, "$mydir\t$h\tAll\t$bhaomega\t\t\tBsu\t$bsuomega");
   } elsif ( $h =~ /H2/ ) { #two omega     
     my @sublines = splice( @results, $#results -2, 2 );
     my @els = split( /Ban#1 #/, $sublines[1]);
     #print @sublines; print "::::\n"; 
     #print $els[1]; 
     my ($banomega, @res) = split (/\s+/, $els[1]);
     my @els2 = split (/Bha #/, $sublines[1]);
     my ($bhaomega, @res2) = split (/\s+/, $els2[1]);
     push (@outlines, "$mydir\t$h\tAll\t$bhaomega\tBan\t$banomega\t\t");     
   } elsif ( $h =~ /H3/ ) { #3 omega
     my @sublines = splice( @results, $#results -2, 2 );
     my @els = split( /Ban#1 #/, $sublines[1]);
     print @sublines; print "::::\n"; 
     print $els[1]; 
     my ($banomega, @res) = split (/\s+/, $els[1]);
     my @els2 = split (/Bha #/, $sublines[1]);
     my ($bhaomega, @res2) = split (/\s+/, $els2[1]);
     my @els3 = split (/Bsu#2 #/, $sublines[1]);
     my ($bsuomega, @res3) = split( /\s+/, $els3[1]);
     push (@outlines, "$mydir\t$h\tAll\t$bhaomega\tBan\t$banomega\tBsu\t$bsuomega");     
   } else { #error
     print "\nWARNING $mydir $h\n";
   } 

   chdir "..";
  }

}

#print "\n--------\n"; print @outlines; #"\n";

chdir $root; 
open( OUT, ">__omega.coat.091008.csv");
foreach my $line (@outlines) {
 print OUT $line."\n";
}
close (OUT);

exit;

#
# END
#

