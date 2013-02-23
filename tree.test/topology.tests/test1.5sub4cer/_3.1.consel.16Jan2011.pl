use strict; use warnings; 
#use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 1 ;

my $inclusterfile = '_test.csv';

open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); shift @lines;
my @clusterdirs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($num, $bg, @res) = split( /\t/, $line);
  push (@clusterdirs, $bg);
}

if ($debug>9) { @clusterdirs = splice(@clusterdirs, 0, 5); }

my $home = "/Users/hongqin"; 
my $root = "$home/coat.protein07/ortholog.analysis/essetial.mcl.rbh.090208/bacillus.phylip.I3.1/topology.tests/test1.5sub4cer"; 

#"/home/hqin/projects/coat.protein07/orthog.analysis/essetial.mcl.rbh.090208/bacillus.phylip.I3.1/paml/test1.5sub4cer/";
#"/home/hqin/projects/coat.protein07/orthog.analysis/merge.mcl.rbh.061208/bacillus.phylip.I3.1/paml/test1.5sub4cer/";

#my $wkdir = "Bha5sub4cer090108";
#my $long_choices = "Bha Bpu Bli Bam Bmo Bsu Bwe Bce Ban Bth";
#my @choices = split( /\s+/, $long_choices);

#my @Hs = ('H0','H1a', 'H1b', 'H1c', 'H1d','H2a','H2b', 'H2c', 'H3a','H3b', 'H3c','H4a');
#my @Hs = ('H3c','H4a');

my @Htres = ( 
#"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,Ban,Bth)));", #H0
#"(Bha,(Bpu,(Bli,(Bam ,(Bmo #1,Bsu #1)))),(Bwe,(Bce,Ban,Bth)));", #H1a
#"(Bha,(Bpu,(Bli,(Bam #1,(Bmo #1,Bsu #1)#1))),(Bwe,(Bce,Ban,Bth)));", #H1b
#"(Bha,(Bpu,(Bli #1,(Bam #1,(Bmo #1,Bsu #1)#1)#1)),(Bwe,(Bce,Ban,Bth)));", #H1c
#"(Bha,(Bpu #1,(Bli #1,(Bam #1,(Bmo #1,Bsu #1)))),(Bwe,(Bce,Ban,Bth)));", #H1d
#"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,Ban #1,Bth)));", #H2a
#"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce #1,Ban #1,Bth #1)));", #H2b
#"(Bha,(Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #1,(Bce #1,Ban #1,Bth #1)#1));", #H2c
#"(Bha,(Bpu,(Bli #2,(Bam #2,(Bmo #2,Bsu #2)#2)#2)),(Bwe #1,(Bce #1,Ban #1,Bth #1)#1));", #H3a
#"(Bha,(Bpu #2,(Bli #2,(Bam #2,(Bmo #2,Bsu #2)#2)#2)#2),(Bwe #1,(Bce #1,Ban #1,Bth #1)#1));" #H3b
 "(Bha,(Bpu ,(Bli,(Bam ,(Bmo ,Bsu )) ) ),(Bwe #2,(Bce #2,Ban #1,Bth #2)#2) );", #H3c
 "(Bha,(Bpu ,(Bli#3, (Bam #3,(Bmo #3 ,Bsu #3)#3) #3) ),(Bwe #2,(Bce #2,Ban #1,Bth #2)#2) );" #H4a
);


#my @Hctls = ( "$root/model.H0.ctl", "$root/model.H1.ctl", "$root/model.H2.ctl", "$root/model.H3.ctl");
#my $Hxctl = "$root/model.HX.ctl";

my $count = 0;

foreach my $mydir ( @clusterdirs ) {
  if (($debug>9) &&($count>2) ) { exit; }
 
  chdir $root;  chomp$mydir;  chdir "wkdir/$mydir";
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

   system("pwd");
   system( "cp trees.nwk $mydir.tpl ");

   # modify lnf file for consel (Add the number of trees to the first line). 
   #system( "cp lnf $mydir.lnf ");
   open (IN, "<lnf"); my @tmplines = <IN>; close (IN); 
   open (OUT, ">$mydir.lnf"); 
   my $firstline = shift @tmplines; 
   chomp $firstline; 
   my @tokens = split( /\s+|\t+/, $firstline); 
   for my $i ( 0..$#tokens ) {
    #if ( $i==1 ) { $tokens[$i] = 9; } ## THERE ARE 9 trees. HARD code!!!!!
    if ( $i==1 ) { $tokens[$i] = 10; } ## THERE ARE 10 trees. HARD code!!!!!
    print OUT "    $tokens[$i]"; 
   }
   foreach my $l (@tmplines) { print OUT $l; }
   close (OUT); 

   system("makermt --paml $mydir.lnf");   
   system("consel $mydir");
   system("catpv $mydir > $mydir.consel.report.txt"); 
   system("ln -s $mydir.pv consel.pv"); 
  
   chdir "../..";

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
