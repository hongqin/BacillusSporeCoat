 
use strict; use warnings; 
#use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 0 ;

my $inclusterfile = '_073108.test.csv';

open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); shift @lines;
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
my $Hxctl = "$root/model.HX.ctl";

my $count = 0;

foreach my $mydir ( @clusterdirs ) {
  if (($debug>9) &&($count>2) ) { exit; }
 
  chdir $root;   chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  chdir $wkdir;

  foreach my $i ( 0 .. $#Hs) {
   my $h = $Hs[$i];
   if ( ! (-e $h) ) { mkdir $h; }
   chdir $h;
   open( OUT, ">tree.nwk"); print OUT $Htres[$i]."\n"; close (OUT);

   #system( "cp $Hctls[$i] codeml.ctl" );
   open (IN, "<$Hxctl"); my @lines2 = <IN>; close (IN);

   $lines2[2] =~ s/HX/$Hs[$i]/o;  #outfile 

   #omega choice
   if ($i == 0) { $lines2[11] =  "       model = 0    \n"; #H0 one omega
   } else { $lines2[11] =  "       model = 2    \n"; }  #alternatives, 2 omega

   open (OUT, ">codeml.ctl"); 
   foreach my $line2 (@lines2) {
     print OUT $line2;
   }
   close(OUT);

   system( "ln -sf ../cds2.phylip cds2.phylip" );
   system( "ln -sf ../cds2.phylip cds.phylip" );

   system( "codeml codeml.ctl");

   chdir "..";
  }

}

exit;

#
# END
#

#my $H0tre = "(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,Ban,Bth)));";
#my $H1tre = "(Bha,((Bpu,(Bli #1,(Bam #1,(Bmo #1,Bsu #1)#1)#1)),(Bwe,(Bce,(Ban,Bth)))));";
#my $H2tre = "(Bha,((Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #1,(Bce #1,(Ban #1,Bth #1)#1)#1)));";
#my $H3tre = "(Bha,((Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,(Ban #1,Bth)))));";

my @___tmpHtres = ( "(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,Ban)));", #H0
"(Bha,(Bpu,(Bli,(Bam ,(Bmo #1,Bsu #1)))),(Bwe,(Bce,Ban)));", #H1a
"(Bha,(Bpu,(Bli,(Bam #1,(Bmo #1,Bsu #1)))),(Bwe,(Bce,Ban)));", #H1b
"(Bha,(Bpu,(Bli #1,(Bam #1,(Bmo #1,Bsu #1)))),(Bwe,(Bce,Ban)));", #H1c
"(Bha,(Bpu #1,(Bli #1,(Bam #1,(Bmo #1,Bsu #1)))),(Bwe,(Bce,Ban)));", #H1d
"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,Ban #1)));", #H2a
"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce #1,Ban #1)));", #H2b
"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #1,(Bce #1,Ban #1)));" #H2c
);
