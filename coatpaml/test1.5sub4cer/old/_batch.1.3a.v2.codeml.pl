#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

# generate sevearl hypothesis testing direcotories. 
# different codeml.ctl
# different trees

my $debug= 1 ;

my $inclusterfile = '_073108.test.csv';
open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); shift @lines;
my @clusterdirs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($bg, @res) = split( /\t/, $line);
  push (@clusterdirs, $bg);
}

my $root = "/home/hqin/projects/coat.protein07/orthog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1/paml/test1.5sub4cer";

my $wkdir = "fulltreeV2.080308";

#my $H0ctl = "$root/model.H0.ctl";
my $Hxctl = "$root/model.HX.ctl";

my @Hs = ('H0','H1a', 'H1b', 'H1c', 'H1d','H2a','H2b', 'H2c');
#my @Hs = ('H0');

#my $H0tre = 

my @Htres = ( "(Bcl,Bha,((Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,(Ban,Bth)))));", #H0
  "(Bcl,Bha,((Bpu,(Bli,(Bam,(Bmo #1,Bsu #1)))),(Bwe,(Bce,(Ban,Bth)))));", #H1a
  "(Bcl,Bha,((Bpu,(Bli,(Bam #1,(Bmo #1,Bsu #1)#1))),(Bwe,(Bce,(Ban,Bth)))));", #H1b
  "(Bcl,Bha,((Bpu,(Bli #1,(Bam #1,(Bmo #1,Bsu #1)#1)#1)),(Bwe,(Bce,(Ban,Bth)))));", #H1c
  "(Bcl,Bha,((Bpu #1,(Bli #1,(Bam #1,(Bmo #1,Bsu #1)#1)#1)#1),(Bwe,(Bce,(Ban,Bth)))));", #H1d
  "(Bcl,Bha,((Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,(Ban #1,Bth #1)))));", #H2a
  "(Bcl,Bha,((Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce #1,(Ban #1,Bth #1)#1))));", #H2b
  "(Bcl,Bha,((Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #1,(Bce #1,(Ban #1,Bth #1)#1)#1)));" #H2c
);

####################################################################

my $count = 0;
foreach my $mydir (@clusterdirs) {
  if (($debug>9) &&($count>2) ) { exit; }
  chdir $root;
  chdir $mydir; 

  print "\n******* count = $count $mydir \n"; #  system( "pwd" );

  if ( !(-e $wkdir) ) { mkdir $wkdir; }
  chdir $wkdir;

  foreach my $i ( 0 .. $#Hs) {
   my $h = $Hs[$i];
   if ( ! (-e $h) ) { mkdir $h; }
   chdir $h;
   open( OUT, ">tree.nwk"); print OUT $Htres[$i]."\n"; close (OUT);

   open (IN, "<$Hxctl"); my @lines2 = <IN>; close (IN);
   $lines2[2] =~ s/HX/$Hs[$i]/o;  #outfile 

   #omega choice #H0 one omega
   if ($i == 0) { $lines2[11] =  "       model = 0    \n"; } 

   open (OUT, ">codeml.ctl"); 
   foreach my $line2 (@lines2) {      print OUT $line2;   }
   close(OUT);

   system( "ln -sf ../../cds.phylip cds.phylip" );
   system( "codeml codeml.ctl");

   chdir "..";
  }

  $count ++;
}
close (IN);

exit;


##