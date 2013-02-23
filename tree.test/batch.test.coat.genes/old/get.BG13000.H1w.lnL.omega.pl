use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; 

my $debug= 1;

system("pwd > __mydir.txt"); open(IN,"<__mydir.txt"); my @tmplines=<IN>; close(IN);
my $root = shift @tmplines; chomp $root; 

my $wkdir = 'Bha5sub4cer';

my $ingrepfile = '_lnL.BG13000.H1w.txt';
open(IN, "<$ingrepfile"); my @lines = <IN>; close (IN);

#open(OUT,">out.Bha5sub4cer/out.BG13000.H1w.Bha5sub4cer.lnL.csv");
#foreach my $line ( @lines ) {  
#  my ($bg, $mydir, $h, $rest) = split( /\//, $line );
#  print OUT "$bg\t$h\t$line";
#}
#close(OUT);

if($debug>9) { @lines = splice( @lines, 0, 10); }

my $count = 0;
my @outlines = ();
my %bgtrace = ();
foreach my $line ( @lines ) {
  chdir $root;
  
  my ($bg, $mydir, $h, $rest) = split( /\//, $line );
  if (!exists($bgtrace{$bg}) ) { 
    $bgtrace{$bg} ++; 
    push (@outlines, '#'.$bg); 
  }

  my (@tokens) = split( /\s+/, $line );
  my $lnL = $tokens[4];

  my (@tok2 ) = split( /\:/, $line );
  my (@tok3 ) = split( /\//, $tok2[0] );
  my $resultFile = $tok3[ $#tok3 ];


  chdir "$bg/$mydir/$h";

   print "\t==>I am in $bg $mydir $h\n";
   system("pwd; ls *txt");
   my $light = 'green'; 
   my @results = (); 
   #unless ( open (RES, "<results.$h.txt") ){ $light='red'};  ### Change here
   unless ( open (RES, "<$resultFile") ){ $light='red'};  ### Change here
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
       push( @outlines, "$bg\t$h\t$lnL\tAll\t".$dN/$dS."\t\t\t\t"); 
     } else {
       push( @outlines, "$bg\t$h\t$lnL\tAll\tNA\t\t\t\t"); 
     } 
   } elsif ( $h =~ /H[12]/ ) { 
     my @specs = ('Bsu','Bwe','Bce','Bth','Ban');
     my @sublines = splice( @results, $#results -2, 2 );
     my $result = "$bg\t$h\t$lnL";
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


chdir $root; 
open( OUT, ">__omega.BG13000.H1w.Bha5sub4cer.011109.csv");
foreach my $line (@outlines) {
 print OUT $line."\n";
}
close (OUT);

exit;

#
# END
#

