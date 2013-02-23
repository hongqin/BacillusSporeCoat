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

my $H0ctl = "$root/model.H0.ctl";
my $H1ctl = "$root/model.H1.ctl";
my $H2ctl = "$root/model.H2.ctl";
#my $H0tre = "(NT01GK,(NT03BC,(BH,((NT03BL,(Bmo,BG)),(BC,(BA,NT02BT))))));";  # 
#my $H0tre = "(Gka,(Bcl,(Bha,((Bli,(Bmo,Bsu)),(Bce,(Ban,Bth))))));";
my $H0tre = "(Bcl,(Bha,((Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe,(Bce,(Ban,Bth))))));";
my $H1tre = "(Bcl,(Bha,((Bpu,(Bli #1,(Bam #1,(Bmo #1,Bsu #1)#1)#1)),(Bwe,(Bce,(Ban,Bth))))));";
my $H2tre = "(Bcl,(Bha,((Bpu,(Bli,(Bam,(Bmo,Bsu)))),(Bwe #1,(Bce #1,(Ban #1,Bth #1)#1)#1))));";

#my $H1tre = "(NT01GK,(NT03BC,(BH,((NT03BL #1,(Bmo #1,BG #1) #1) #1,(BC,(BA,NT02BT))))));"; 
#my $H1tre = "(Gka,(Bcl,(Bha,((Bli #1,(Bmo #1,Bsu #1) #1 )#1,(Bce,(Ban,Bth))))));";
#my $H2tre = "(Gka,(Bcl,(Bha,((Bli,(Bmo,Bsu)),(Bce #1,(Ban #1,Bth #1) #1)#1 ))));";

my $longline = "Bcl     NT03BC
Bha     BH
Bli     NT03BL
Bmo     Bm
Bpu     BP
Bam     RB
Bsu     BG
Bce     BC
Ban     BA
Bth     NT02BT
Bwe     Bc
";#I only letters that can sufficiently distinguish unique genomes
# Bwe has to follow Bce

my %prefix_hash = ();
 _set_prefix_hash (\%prefix_hash, $longline);
showHashTable( \%prefix_hash);

####################################################################

my $count = 0;
foreach my $mydir (@clusterdirs) {
  if (($debug>9) &&($count>1) ) { exit; }
  chdir $root;
  chdir $mydir; 

  print "\n******* count = $count $mydir \n"; #  system( "pwd" );

  my @orfs=();
  open(IN, "<_ids.tab"); 
  while( my $line = <IN>) {
    if ( $line !~ /^\s*$/ ) {
       chomp $line; push (@orfs, $line);
    }
  }
  close (IN);

  ###H0
  mkdir "H0";
  my $newh0tre = _replace_ids_in_tree2( $H0tre, \@orfs, \%prefix_hash);
  open( OUT, ">H0/tree.H0.tre"); print OUT $newh0tre."\n"; close (OUT);  
  system( "cp $H0ctl H0/codeml.ctl" );
  system( "ln -sf ../cds.phylip H0/." );

  ###H1
  mkdir "H1";
  my $newh1tre = _replace_ids_in_tree2( $H1tre, \@orfs, \%prefix_hash);
  open( OUT, ">H1/tree.H1.tre"); print OUT $newh1tre."\n"; close (OUT);
  system( "cp $H1ctl H1/codeml.ctl" );
  system( "ln -sf ../cds.phylip H1/." );

  ###H2
  mkdir "H2";
  my $newh2tre = _replace_ids_in_tree2( $H2tre, \@orfs, \%prefix_hash);
  open( OUT, ">H2/tree.H2.tre"); print OUT $newh2tre."\n"; close (OUT);
  system( "cp $H2ctl H2/codeml.ctl" );
  system( "ln -sf ../cds.phylip H2/." );

  ### 
  $count ++;
}
close (IN);

exit;

#######################################
#  my $newh2tre = _replace_ids_in_tree2( $H2tre, \@ids, \%idhash);
sub _replace_ids_in_tree2 {
 my $debug = 1;
 my ($tre, $reforfs, $refhash) = @_;
 my @orfs = @$reforfs;
 my %hash = %$refhash; 
 if ($debug>2) { print "[ @orfs ]\n"; showHashTable( \%hash ); }

 my %spec2orf = (); #species name to orf 
 my $num = 1;
 foreach my $orf ( @orfs ) {
  	my $prefix = substr($orf, 0, 2);
  	if ($prefix eq "NT") { $prefix = substr( $orf, 0, 6); }
  	$spec2orf{ $hash{$prefix} } = $orf; 
	if($debug) { print " $num\tprefix[$prefix]\torf[$orf]\tid[$hash{$prefix}]\n"; $num++;}
	
 }

 if($debug) { showHashTable(\%spec2orf ); }

 #replace species in the tree in specified order to avoid double-hits, such Bwe->Bce errors
 my $order = "Bcl Bha Bli Bmo  Bpu Bam Bsu Bce Ban Bth Bwe";
 my @orders = split( /\s+/, $order) ;

 foreach my $spec (@orders) {
   	$tre =~ s/$spec/$spec2orf{$spec}/g;
 	if($debug>1) { print " $spec tre [$tre]\n"; }
 }

 return $tre;
}

###############################33
sub _set_prefix_hash {
 my ($rf, $idline) = @_;
 my @lines = split( /\n/, $idline );
 foreach my $line (@lines) {
  my ($species, $prefix, @rest) = split (/\s+/, $line );
  $rf->{$prefix} = $species; 
 }
}
