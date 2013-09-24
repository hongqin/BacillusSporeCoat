#!/usr/bin/perl
use strict; use warnings;
use lib '/home/hqin/lib/perl'; 
use Util;

# generate sevearl hypothesis testing direcotories. 
# different codeml.ctl
# different trees

my $debug = 0;

my $root = "/home/hqin/coat.protein07/orthog.analysis/phylo122307";

my $coatfl = "$root/_simple.coat.orth.list.122307.versionD.csv";
if ($debug) { $coatfl = "_test.csv"; }

my $H0ctl = "$root/model.H0.ctl";
#my $H0tre = "(NT01GK,(NT03BC,(BH,((NT03BL,(Bmo,BG)),(BC,(BA,NT02BT))))));";  # 
my $H0tre = "(Gka,(Bcl,(Bha,((Bli,(Bmo,Bsu)),(Bce,(Ban,Bth))))));";

my $H1ctl = "$root/model.H1.ctl";
#my $H1tre = "(NT01GK,(NT03BC,(BH,((NT03BL #1,(Bmo #1,BG #1) #1) #1,(BC,(BA,NT02BT))))));"; 
my $H1tre = "(Gka,(Bcl,(Bha,((Bli #1,(Bmo #1,Bsu #1) #1 )#1,(Bce,(Ban,Bth))))));";
 # different omega in bsu clade

my $H2ctl = "$root/model.H2.ctl";
my $H2tre = "(Gka,(Bcl,(Bha,((Bli,(Bmo,Bsu)),(Bce #1,(Ban #1,Bth #1) #1)#1 ))));";
 # different omega in bce clase 


my $idline = "Gka     NT01GK
Bcl     NT03BC
Bha     BH
Bli     NT03BL
Bmo     Bm
Bsu     BG
Bce     BC
Ban     BA
Bth     NT02BT
";

my %idhash = ();
 _set_id_hash (\%idhash, $idline);
#showHashTable( \%idhash);

####################################################################

open (IN, "<$coatfl");
my $count = 0;
while (my $line = <IN>) {
  chomp $line; 
  my ( $bg, @ids) = split ( /\t/, $line );

  chdir $root;
  chdir $bg; 

   system( "pwd" );

  unshift( @ids, $bg );

  ###H0
  mkdir "H0";
  my $newh0tre = _replace_ids_in_tree2( $H0tre, \@ids, \%idhash);
  open( OUT, ">H0/tree.H0.tre"); print OUT $newh0tre."\n"; close (OUT);  
  system( "cp $H0ctl H0/codeml.ctl" );
  system( "ln -sf ../cds.phylip H0/." );

  ###H1
  mkdir "H1";
  my $newh1tre = _replace_ids_in_tree2( $H1tre, \@ids, \%idhash);
  open( OUT, ">H1/tree.H1.tre"); print OUT $newh1tre."\n"; close (OUT);
  system( "cp $H1ctl H1/codeml.ctl" );
  system( "ln -sf ../cds.phylip H1/." );

  ###H2
  mkdir "H2";
  my $newh2tre = _replace_ids_in_tree2( $H2tre, \@ids, \%idhash);
  open( OUT, ">H2/tree.H2.tre"); print OUT $newh2tre."\n"; close (OUT);
  system( "cp $H2ctl H2/codeml.ctl" );
  system( "ln -sf ../cds.phylip H2/." );

  ### 
  $count ++;
  print "**** count = $count \n\n";
}
close (IN);

exit;

#######################################
#  my $newh2tre = _replace_ids_in_tree2( $H2tre, \@ids, \%idhash);
sub _replace_ids_in_tree2 {
 my $debug = 0;
 my ($tre, $refids, $refhash) = @_;
 my @orfs = @$refids;
 my %hash = %$refhash; 
 if ($debug) { print "[ @orfs ]\n"; showHashTable( \%hash ); }

 foreach my $orf ( @orfs ) {
  	my $prefix = substr($orf, 0, 2);
  	if ($prefix eq "NT") { $prefix = substr( $orf, 0, 6); }
	if($debug) { print "\tprefix[$prefix]\torf[$orf]\tid[$hash{$prefix}]\n"; }
   	$tre =~ s/$hash{$prefix}/$orf/g;
 	if($debug) { print "\ttre [$tre]\n"; }
 }
 return $tre;
}

###############################33
#  _set_id_hash (\%idhash, $idlines); 

sub _set_id_hash {
 my ($rf, $idline) = @_;
 my @lines = split( /\n/, $idline );
 foreach my $line (@lines) {
  my ($id, $prefix, @rest) = split (/\s+/, $line );
  $rf->{$prefix} = $id; 
 }
}
