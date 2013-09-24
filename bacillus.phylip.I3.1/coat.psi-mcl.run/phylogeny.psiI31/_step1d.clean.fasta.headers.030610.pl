# 030610
# 062208 pick only the ids for fasta header, which is required for PAUP-nj
# 
use strict; use warnings; 
#use lib '/home/hqin/lib/perl'; use Util; 
use lib '/Users/hongqin/lib/perl'; use Util; 

my $debug= 2 ;

my $tmpfl = "/tmp/_coatBGclus.062208.txt";
  # system( "ls -d BG* | head -n 2 > $tmpfl " );
system( "ls -d BG*  > $tmpfl " );
open (IN, "<$tmpfl" ); my @clusterdirs = <IN>; close (IN);
system( "rm $tmpfl" );

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  system( "clean_fasta_headers.01.pl -if ids.faa -of prot.faa -pos 0 " );
  system( "cp prot.faa prot.fas " );

  system( "clean_fasta_headers.01.pl -if ids.fna -of cds.fna  -pos 0 " );
  system( "cp cds.fna cds.fas " );

  chdir "..";
}

exit;

#
# END
#

########################
  #system( "clean_fasta_headers.00.pl -if ids.faa -of /tmp/_tmp.faa -pos 0-1 " );
  #system( "remove_1stSpace_in_fasta_header.00.pl -if /tmp/_tmp.fna -of /tmp/_tmp2.fna " );

