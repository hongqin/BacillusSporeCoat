 
use strict; use warnings; 
#use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 1 ;

my $root = '/home/hqin/coat.protein07/orthog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1/paml/BhaBcl4sub4cer';
my $wkdir = 'BhaBcl4sub4cer';
my $inprofilefile = '_11genome.codeml.coat.csv';

open (IN,"<$inprofilefile"); my @lines = <IN>; close (IN);
my $header = shift @lines; 

if($debug>9) { @lines = splice (@lines, 0, 1); }

#my @Hs = ('H0', 'H1w','H1c','H2C2');
my @Hs = ('H2Cw', 'H2Cc');

my @Htres = ( 
#"((Bha,Bcl),(Bpu,(Bli,(Bam,Bsu))),(Bwe,(Bce,(Ban,Bth))));", #H0
#"((Bha,Bcl),(Bpu,(Bli,(Bam,Bsu))),(Bwe #1,(Bce,(Ban,Bth))));", #H1w
#"((Bha,Bcl),(Bpu,(Bli,(Bam,Bsu))),(Bwe,(Bce #1,(Ban,Bth))));", #H1c
#"((Bha,Bcl),(Bpu,(Bli,(Bam,Bsu))),(Bwe #1,(Bce #1,(Ban #2,Bth #2)#2)#1));", #H2C2
"((Bha,Bcl),(Bpu,(Bli,(Bam,Bsu))),(Bwe #1,(Bce #2,(Ban #2,Bth #2)#2)#2));", #H2Cw
"((Bha,Bcl),(Bpu,(Bli,(Bam,Bsu))),(Bwe #2,(Bce #1,(Ban #2,Bth #2)#2)#2));", #H2Cc
);

my $Hxctl = "$root/model.HX.ctl";


my $count = 0;
foreach my $line ( @lines ) { 
  chdir $root;
  my @els = split( /\t/, $line );
  my $mydir = $els[0];
  $mydir =~ s/\"//g;
  chdir $mydir;
 
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  chdir $wkdir;

  foreach my $i ( 0 .. $#Hs) {
   my $h = $Hs[$i];
   if ( ! (-e $h) ) { mkdir $h; }
   chdir $h;
   open( OUT, ">tree.$h.nwk"); print OUT $Htres[$i]."\n"; close (OUT);

   #system( "cp $Hctls[$i] codeml.ctl" );
   open (IN, "<$Hxctl"); my @lines2 = <IN>; close (IN);

   $lines2[1] =~ s/HX/$Hs[$i]/o;  #tree file
   $lines2[2] =~ s/HX/$Hs[$i]/o;  #report file 

   #omega choice
   if ($h eq 'H0' ) { $lines2[11] =  "       model = 0    \n"; #H0 one omega
   } else { $lines2[11] =  "       model = 2    \n"; }  #alternatives, 2 omega

   open (OUT, ">codeml.ctl"); 
   foreach my $line2 (@lines2) {
     print OUT $line2;
   }
   close(OUT);

   system( "ln -sf ../cds2.phylip cds2.phylip" );
#   system( "ln -sf ../cds2.phylip cds.phylip" );

   system( "codeml codeml.ctl");

   chdir "..";
  }

}

exit;

#
# END
#
