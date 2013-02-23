#site models 
use strict; use warnings; 
#use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 1 ;

my $inprofilefile = '_10g.coat.profile.011009.csv';

system("pwd > __mydir.txt"); open(IN,"<__mydir.txt"); my @tmplines=<IN>; close(IN);
my $root = shift @tmplines; chomp $root; 

my $wkdir = '4cer';

open (IN,"<$inprofilefile"); my @lines = <IN>; close (IN);
my $header = shift @lines; 

if($debug>9) { @lines = splice (@lines, 0, 1); }

#my @Hs = ('H0', 'H1w','H1c','H2Ccw','H2Cw', 'H2Cc');
my @Hs = ('M1a','M2a');

my @Htres = ( "(Bwe,(Bce,(Ban,Bth)));",
"(Bwe,(Bce,(Ban,Bth)));");

my $Hxctl = "$root/site.model.MX.ctl";

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
   open( OUT, ">tree.cereus.$h.nwk"); print OUT $Htres[$i]."\n"; close (OUT);

   open (IN, "<$Hxctl"); my @lines2 = <IN>; close (IN);

   $lines2[1] =~ s/MX/$Hs[$i]/o;  #tree file
   $lines2[2] =~ s/MX/$Hs[$i]/o;  #report file 

   #omega choice
   if ($h eq 'M1a' ) { 
     $lines2[11] =  "       model = 0    \n"; # M1a
     $lines2[15] =  "      NSsites = 1     * \n";
   } elsif ($h eq 'M2a') { #M2a 
      $lines2[11] =  "       model = 0    \n"; #M2a
      $lines2[15] =  "      NSsites = 2     * \n";
    }

   open (OUT, ">codeml.ctl"); 
   foreach my $line2 (@lines2) {
     print OUT $line2;
   }
   close(OUT);

   system( "ln -sf ../cds2.phylip cds2.phylip" );

   system( "codeml codeml.ctl");

   chdir "..";
  }

}

exit;

#
# END
#
