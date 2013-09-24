# 070508 bootstrap nj
# 062208 protein nj in PAUP
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 2 ;

#my $idfl   = "/home/hqin/projects/coat.protein07/key.data/_id2genomes.csv";
#my $specfl = "/home/hqin/projects/coat.protein07/key.data/_genome2species.csv";
my $clusterfl = "infile-BGclusters.bootstrap.txt";

#my $tre ="((((((Bmo,Bsu),Bam),Bli),Bpu),((Bce,Ban,Bth),Bwe)),(Bha,Bcl))";
# OR in a more readble format: 
#	( 
#	(
#	 	( ( ( (Bmo,Bsu),Bam), Bli), Bpu ), 
#		( (Bce,Ban, Bth),Bwe)
#	), (Bha, Bcl) 
#	);

open (IN, "<$clusterfl" ); my @clusterdirs = <IN>; close (IN);

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  open  (OUT, ">seqCat.in");   print OUT "prot-bacillus-aligned.fasta\n";   close (OUT);  
  my $cmd1 = "seqCat.pl -dseqCat.in -if";
  system( $cmd1 );
  system( "mv seqCat_sequences.nex _tmp.prot-paup.nex "  );
  open (IN, "<_tmp.prot-paup.nex"); my @lines = <IN>; close (IN);
  pop @lines; ### remove the end;
  pop @lines; ### remove charset
  pop @lines; ### remove begin paups;

  open (OUT, ">prot-bacillus-paup.nexus"); 
  foreach my $line (@lines) { 
    if ( $line =~ "datatype = DNA" ) { 
        $line =~ s/DNA/protein/o;
    } 
    print OUT $line; 
  }


#   savetrees file=prot-bacillus-paup-njbt-brlens.tre SaveBootP=both from=1 to=1;
# root with the 1st taxon

  my $pauplines = "
BEGIN PAUP;
   outgroup 1;
   NJ;
   ROOTTREES ROOTMETHOD=OUTGROUP;
   savetrees file=prot-bacillus-paup-nj-brlens.tre brlens=YES;
   savetrees file=prot-bacillus-paup-nj.tre brlens=NO;

   Bootstrap nreps=500 treefile = bootstrap1.tre replace=no brlen=yes/ start=stepwise addseq=random;
   Contree all /majrule=yes strict=no treefile=finalbt.tre grpfreq=yes;
   savetrees file=prot-bacillus-paup-njbt.tre brlens=NO from=1 to=1 SaveBoot=NodeLabel;
END;
  ";
  print OUT "$pauplines";
  close (OUT);
 
  system( "paup prot-bacillus-paup.nexus -n " );

  if ( ! -e "prot-paup-nj.tre" ) {
    system("touch __FAILED.prot-bacillus-nj.txt");
  }

  chdir "..";
}

exit;

#
# END
#

