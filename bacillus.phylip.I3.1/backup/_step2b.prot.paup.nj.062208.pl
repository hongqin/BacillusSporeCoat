# 062208 generate PAUP nexus file for protein nj
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; 

my $debug= 2 ;

my $tmpfl = "/tmp/_coatBGclus.062208.txt";
 #   system( "ls -d BG* | head -n 2 > $tmpfl " );
system( "ls -d BG*  > $tmpfl " );
open (IN, "<$tmpfl" ); my @clusterdirs = <IN>; close (IN);
system( "rm $tmpfl" );

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  open  (OUT, ">seqCat.in");   print OUT "prot-aligned.fasta\n";   close (OUT);  
  my $cmd2 = "seqCat.pl -dseqCat.in -if";
  system( $cmd2 );
  system( "mv seqCat_sequences.nex _tmp.prot-paup.nex "  );
  open (IN, "<_tmp.prot-paup.nex"); my @lines = <IN>; close (IN);
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
   savetrees file=prot-paup-nj-brlens.tre brlens=YES;
   savetrees file=prot-paup-nj.tre brlens=NO;
END;
  ";
  print OUT "$pauplines";
  close (OUT);
 
  system( "paup prot-paup.nexus -n " );

  if ( ! -e "prot-paup-nj.tre" ) {
    system("touch __FAILED.prot-nj.txt");
  }

  chdir "..";
}

exit;

#
# END
#

########################

