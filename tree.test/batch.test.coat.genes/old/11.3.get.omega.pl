use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; 

my $debug= 1;

my $root = '/home/hqin/coat.protein07/orthog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1/paml/BhaBcl';
my $wkdir = '11g';

my $ingrepfile = 'out.11g/grep.lnL.txt';
open(IN, "<$ingrepfile"); my @lines = <IN>; close (IN);

open(OUT,">out.11g/out.coat.11glnL.csv");
foreach my $line ( @lines ) {  
  my ($bg, $mydir, $h, $rest) = split( /\//, $line );
  print OUT "$bg\t$h\t$line";
}
close(OUT);

if($debug>9) { @lines = splice( @lines, 0, 10); }

#my @Hs = ('H0','H1a','H1c','H1t','H1w','H2C2a','H2C2t','H2C2c','H2C2w','H2C2cw');

#my @Htres = ( 
#"((Bha,Bcl),(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,(Ban,Bth))));", #H0
#"((Bha,Bcl),(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,(Ban #1,Bth))));", #H1a
#"((Bha,Bcl),(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce #1,(Ban,Bth))));", #H1c
#"((Bha,Bcl),(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,(Ban,Bth #1))));", #H1t
#"((Bha,Bcl),(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #1,(Bce,(Ban,Bth))));", #H1w
#"((Bha,Bcl),(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #1,(Bce #1,(Ban #2,Bth #1)#1)#1));", #H2C2a
#"((Bha,Bcl),(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #1,(Bce #1,(Ban #1,Bth #2)#1)#1));", #H2C2t
#"((Bha,Bcl),(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #1,(Bce #2,(Ban #1,Bth #1)#1)#1));", #H2C2c
#"((Bha,Bcl),(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #2,(Bce #1,(Ban #1,Bth #1)#1)#1));", #H2C2w
#"((Bha,Bcl),(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #2,(Bce #2,(Ban #1,Bth #1)#1)#2));", #H2C2cw
#);

my $count = 0;
my @outlines = ();
my %bgtrace = ();
foreach my $line ( @lines ) {
  chdir $root;
  
  my ($bg, $mydir, $h, $rest) = split( /\//, $line );
  #if (!exists($bgtrace{$bg}) ) { 
  #  $bgtrace{$mydir} ++; 
  #  push (@outlines, '#'.$mydir); 
  #}

  chdir "$bg/$mydir/$h";

   print "\t==>I am in $bg $mydir $h\n";
   system("pwd; ls *txt");
   my $light = 'green'; 
   my @results = (); 
   unless ( open (RES, "<results.$h.txt") ){ $light='red'};  ### Change here
   if ( $light eq 'green' ) {
      @results = <RES>; 
   close (RES);
   } else {
     system ("touch /tmp/error-$bg-$mydir-$h"); 
   } 

   if ( $light eq 'green' ) {
   if ( $h eq 'H0' ) {
     my @sublines = splice( @results, $#results -2, 2 );
     my ($tmp, $dN) = split( /:\s+/, $sublines[0]);
     my ($tmp2, $dS) = split( /:\s+/, $sublines[1]);
     #print @sublines."\n"; print "\$dN=$dN  \$dS=$dS\n";
     if (defined $dS) {
       push( @outlines, "$bg\t$h\tAll\t".$dN/$dS."\t\t\t\t"); 
     } else {
       push( @outlines, "$bg\t$h\tAll\tNA\t\t\t\t"); 
     } 
   } elsif ( $h =~ /H[12]/ ) { 
     my @specs = ('Bsu','Bwe','Bce','Bth','Ban');
     my @sublines = splice( @results, $#results -2, 2 );
     my $result = "$bg\t$h";
     foreach my $spec( @specs ) {
        my @els = split( /$spec[#\d|\s]+#/, $sublines[1]);
        if(defined $els[1]) {
          my ($specomega, @res) = split (/\s+/, $els[1]);
          $result = $result."\t$spec\t$specomega";
        } else {
         $result =  $result."\t$spec\tNA";
        }
     }
     push (@outlines, $result);      
   } else { #error
     print "\nWARNING $mydir $h\n";
   } 
   } 

}

#print "\n--------\n"; print @outlines; #"\n";

chdir $root; 
open( OUT, ">__omega.coat.11g.011009.csv");
foreach my $line (@outlines) {
 print OUT $line."\n";
}
close (OUT);

exit;

#
# END
#

