#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

#protein ML tree in paup.

my $debug = 0;

my $rootdir = '/home/hqin/coat.protein07/orthog.analysis/phylo012708';

### get the coat ids
system("ls -d BG* > /tmp/_ls-file.txt");
my %coats = ();
open (IN, "</tmp/_ls-file.txt");
while (my $line = <IN> ) {
   chomp $line; 
   my($bg, @els) = split (/-/, $line); 
   $coats{$bg} = $bg;
}
close(IN);
if ($debug) { showHashTable( \%coats); }

my $count = 0; 
foreach my $bg ( sort ( keys %coats ) ) {
  chdir $rootdir;
  chdir $bg; 
  if ($debug) {    system( "pwd" );  }

  open (IN2, "<prot.fasta");   my @tmps = <IN2>;   close (IN2);
  my @headers = grep (/^>/, @tmps );

  if ( (scalar @headers ) >= 3 ) { 

  ### generate nj 
  open (IN, "<prot-paup.nex"); my @lines = <IN>; close (IN);
  pop @lines; ### remove the end;
  pop @lines; ### remove charset
  pop @lines; ### remove begin paups;

  open (OUT, ">prot-paup.nexus"); 
  foreach my $line (@lines) { 
    if ( $line =~ "datatype = DNA" ) { 
	$line =~ s/DNA/protein/o;
    } 
    print OUT $line; 
  }

  my $pauplines = "
BEGIN PAUP;
   NJ;
   savetrees file=prot-paup-nj.tre brlens;
END;
  ";
  print OUT "$pauplines";
  close (OUT);
 
  system( "paup prot-paup.nexus -n " );

  ###########boot strap
  open (OUT, ">prot-paup2.nexus"); 
  foreach my $line (@lines) { 
    if ( $line =~ "datatype = DNA" ) { 
	$line =~ s/DNA/protein/o;
    } 
    print OUT $line; 
  }

  my $pauplines2 = "
BEGIN PAUP;
   NJ;
   bootstrap;
   savetrees file=prot-paup-nj-btstrp.tre brlens;
END;
  ";
  print OUT "$pauplines2";
  close (OUT);
  } else {
     print "I am skipping $bg because it has ".(scalar @headers)." hits\n"
  }
  #system( "paup prot-paup2.nexus -N" );

  $count ++;   print "**** count = $count $bg \n\n";
  
}


exit;


