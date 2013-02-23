use lib '/Users/hongqin/lib/perl';use Util;
use strict; use warnings; 

my $debug = 0;
my $mapfile = "gene.id.summary.062210.tab";

my $user = "/Users/hongqin";
my $wkdir = "$user/coat.protein07/orthog.analysis/pairwise.matrix.omega.062310";
my $IDdir = "$user/coat.protein07/key.data.062210";

open(ERROR, ">_error.dump.txt");

my ( %strains) = (); #store file names for coat-id-file, essen-id-file, nonCE-id-file, and gene-id Abbreviations

#my %pvalues = (); ##THIS is the p-value MATRIX !!!
my %pvalues_w_coat = (); ##THIS is the p-value MATRIX !!!
my %pvalues_dN_coat = (); ##THIS is the p-value MATRIX !!!
my %pvalues_dS_coat = (); ##THIS is the p-value MATRIX !!!
my %pvalues_w_essen = (); ##THIS is the p-value MATRIX !!!
my %pvalues_dN_essen = (); ##THIS is the p-value MATRIX !!!
my %pvalues_dS_essen = (); ##THIS is the p-value MATRIX !!!

open(IN, "<$mapfile"); my @maplines =<IN>; close(IN);
shift @maplines;
foreach my $mapline (@maplines) {
  chomp $mapline; 
  my ($coatfile, $essenfile, $nonCEfile, $abbr, $genome, @res) = split( /\s+/, $mapline ); 
  my ($strain, @res2) = split( /\//, $genome );
  $strains{ $strain } = {
			 'COAT'=> $coatfile,
			 'ESSEN' => $essenfile,
			 'NONCE'  => $nonCEfile,
			 'ABBR' => $abbr
			}
}

if($debug>19) {
  my $i = 1; 
  foreach my $k (sort keys %strains) {
    print ERROR "\n\n$i coatfile for [$k] is [".$strains{$k}->{COAT}."]\n";
    #system("head -n 3 $IDdir/coat.genes/".$strains{$k}->{COAT}); 
    print ERROR "$i Essenfile for [$k] is [".$strains{$k}->{ESSEN}."]\n";
    #system("head -n 3 $IDdir/essential.genes/".$strains{$k}->{ESSEN}); 
    print ERROR "$i nonCEfile for [$k] is [".$strains{$k}->{NONCE}."]\n";
    system("head -n 3 $IDdir/noncoat.nonessen.genes/".$strains{$k}->{NONCE}); 
    $i ++; 
  }
}


#"-e" does not work, so I have to get the list of dNdS files
system( "ls dNdS.files/_dNdS*csv > /tmp/_062310.dNdS.list.txt");
my %dNdSFiles = ();
open(TMP, "</tmp/_062310.dNdS.list.txt");
while (my $dNdSline = <TMP>) {
  $dNdSline =~ s/dNdS.files\///;
  chomp $dNdSline;  
  $dNdSFiles{ $dNdSline } ++; 
}
close(TMP);
showHashTable( \%dNdSFiles );

#I will loop through the pairs,  find and partition omega results,  do t.test, 
# Finally, output to a matrix.

my @s = sort keys %strains;
if ($debug > 5) { @s = splice( @s, 0, 3); } 

for (my $i=1; $i<=$#s; $i++ ){
  for( my $j=0; $j< $i; $j++ ) {

    #Find the dNdS files
    my $dNdSFile1 = "_dNdS.rpbh.".$s[$i]."-".$s[$j].".csv";
    my $dNdSFile2 = "_dNdS.rpbh.".$s[$j]."-".$s[$i].".csv";
    my $dNdSFile = "NA";
    my $index = -1;
    # if ( -e "dNdS.files/$dNdSFile2") { $dNdSFile = $dNdSFile2; } #This is not working???
    # I have to call system to get the file list into a hash
    if ( exists $dNdSFiles{ $dNdSFile2 } ) { $dNdSFile = "dNdS.files/".$dNdSFile2; $index = $j}
    elsif (exists $dNdSFiles{ $dNdSFile1 } ) { $dNdSFile = "dNdS.files/".$dNdSFile1; $index =$i }

    if( $debug>20) { 
      system( "head -n 3 $dNdSFile" );  
    }
    print ERROR "\n\n****************\$i= $i, \$j=$j, $s[$i] $s[$j] -> [$dNdSFile]\n";		   

    #### read omega into omega
    my %omega = ();    my %dN = ();    my %dS = (); 
    open (OMEGA, "<$dNdSFile");
    while( my $oline = <OMEGA> ){
      chomp $oline;
      my @els = split( /\s+/, $oline );
      if ($debug> 15) {
	for( my $k=0; $k<= $#els; $k++ ) {
	  print "\n$k $els[$k]"; 
	}
      }

      if ($#els == 9) {
	my $id_pos = -1;
	if ($els[0] =~ /$strains{ $s[$index] }->{ABBR}/) { $id_pos=0 }
        elsif ( $els[1] =~ /$strains{ $s[$index] }->{ABBR}/ ) { $id_pos=1 } 
	my $current_id = "NA";
	if ($id_pos >= 0 ) { $current_id = $els[$id_pos]; }
	$dN{$current_id} = $els[5];
	$dS{$current_id} = $els[7];
	$omega{$current_id} = $els[9];
      } else { #print out the errors
	if ($debug > 15) { print ERROR "$oline\n"; }  
      }
    }
    close (OMEGA);
    # showHashTable( \%omega ); 

    ##### read coat ids, get their omegas
    my $coatidfile = $IDdir."/coat.genes/".$strains{ $s[$index] }->{COAT}; 
    print ERROR "\ncoatidfile is [$coatidfile]\r\n";
    my %coatW = (); my %coatdN =(); my %coatdS = (); 
    open (COAT, "<$coatidfile");
    while (my $cline = <COAT> ) {
      chomp $cline; 
      if ( $cline !~ /^\s*$/ ) {
       my @els = split (/\s+/, $cline);
       if ( $els[0] !~ /^\s*$/ ) {
	$coatW{ $els[0] } = 'NA'; 
	$coatdN{ $els[0] } = 'NA'; 
	$coatdS{ $els[0] } = 'NA'; 
	if (exists $omega{ $els[0] } ) {
	  $coatW{ $els[0] } = $omega{ $els[0] }  ; 
	  $coatdN{ $els[0] } = $dN{ $els[0] }  ; 
	  $coatdS{ $els[0] } = $dS{ $els[0] }  ; 
	}
       }
      }
    }
    close (COAT);
    #showHashTable( \%coatW );    #showHashTable( \%coatdN );    showHashTable( \%coatdS);

    #### read in essential gene ids, get their omegas
    my $essenidfile = $IDdir."/essential.genes/".$strains{ $s[$index] }->{ESSEN}; 
    print ERROR "\n Essenidfile is [$essenidfile]\r\n";
    my %essenW = ();     my %essendN = ();  my %essendS = (); 
    open (ESSEN, "<$essenidfile");
    while (my $eline = <ESSEN> ) {
      chomp $eline; 
      my @els = split (/\s+/, $eline);
      if ( $els[0] !~ /^\s*$/ ) {
	$essenW{ $els[0] } = 'NA';
	$essendS{ $els[0] } = 'NA';
	$essendN{ $els[0] } = 'NA';
	if (exists $omega{ $els[0] } ) {
	  $essenW{ $els[0] } = $omega{ $els[0] }  ; 
	  $essendS{ $els[0] } = $dN{ $els[0] }  ; 
	  $essendN{ $els[0] } = $dS{ $els[0] }  ; 
	}
      }
    }
    close (ESSEN);
    #showHashTable( \%essenW );

    #### read in nonCE gene id, get their omegas
    my $nonCEidfile = $IDdir."/noncoat.nonessen.genes/".$strains{ $s[$index] }->{NONCE}; 
    print ERROR "\n nonCEidfile is [$nonCEidfile]\r\n";
    my %nonCEW = ();my %nonCEdN = (); my %nonCEdS=();  
    open (NONCE, "<$nonCEidfile");
    while (my $eline = <NONCE> ) {
      chomp $eline; 
      my @els = split (/\s+/, $eline);
      if ( $els[0] !~ /^\s*$/ ) {
	$nonCEW{ $els[0] } = 'NA';
	$nonCEdN{ $els[0] } = 'NA';
	$nonCEdS{ $els[0] } = 'NA';
	if (exists $omega{ $els[0] } ) {
	  $nonCEW{ $els[0] }  = $omega{ $els[0] }  ; 
	  $nonCEdN{ $els[0] } = $dN{ $els[0] }  ; 
	  $nonCEdS{ $els[0] } = $dS{ $els[0] }  ; 
	}
      }
    }
    close (NONCE);
    #showHashTable( \%nonCEW );

    ####  t.test through R
    my ($bigger, $smaller) = order_big2small($s[$i] ,  $s[$j] ); 

    my @coatW = values %coatW;  print ERROR "\nOmega::". $coatW[0] . "*"; 
    my @essenW = values %essenW; print ERROR $essenW[0]. "*";  
    my @nonCEW = values %nonCEW; print ERROR $nonCEW[0]."*\n"; 
    $pvalues_w_coat { "$bigger, $smaller" } = _pairwise_test (\@coatW, \@nonCEW, *ERROR); 
    $pvalues_w_essen { "$bigger, $smaller" } = _pairwise_test (\@essenW, \@nonCEW, *ERROR); 

    my @coatdN = values %coatdN; print ERROR "\ndN:: ".$coatdN[0]."*";
    my @essendN = values %essendN; print ERROR $essendN[0]."*"; 
    my @nonCEdN= values %nonCEdN;  print ERROR $nonCEdN[0]."*\n"; 
    $pvalues_dN_coat { "$bigger, $smaller" } = _pairwise_test (\@coatdN, \@nonCEdN, *ERROR); 
    $pvalues_dN_essen { "$bigger, $smaller" } = _pairwise_test (\@essendN, \@nonCEdN, *ERROR); 

    my @coatdS = values %coatdS; print ERROR "\ndS::".$coatdS[0]."*";
    my @essendS = values %essendS; print ERROR $essendS[0]."*";
    my @nonCEdS= values %nonCEdS;  print ERROR $nonCEdS[0]."*\n";
    $pvalues_dS_coat { "$bigger, $smaller" } = _pairwise_test (\@coatdS, \@nonCEdS, *ERROR); 
    $pvalues_dS_essen { "$bigger, $smaller" } = _pairwise_test (\@essendS, \@nonCEdS, *ERROR); 
  }
}

close (ERROR);


print "\n\n***Omega***"; 
showHashTable( \%pvalues_w_coat ); 
_print_pvalues_dataframe_for_R( \"pmat.omega-coat-nonCE-wilcox-greater.tab",\%pvalues_w_coat ); 
showHashTable( \%pvalues_w_essen ); 
_print_pvalues_dataframe_for_R( \"pmat.omega-essen-nonCE-wilcox-greater.tab",\%pvalues_w_essen ); 

print "\n\n***dN***"; 
showHashTable( \%pvalues_dN_coat); 
_print_pvalues_dataframe_for_R( \"pmat.dN-coat-nonCE-wilcox-greater.tab",\%pvalues_dN_coat ); 
showHashTable( \%pvalues_dN_essen); 
_print_pvalues_dataframe_for_R( \"pmat.dN-essen-nonCE-wilcox-greater.tab",\%pvalues_dN_essen ); 

print "\n\n***dS***"; 
showHashTable( \%pvalues_dS_coat); 
_print_pvalues_dataframe_for_R( \"pmat.dS-coat-nonCE-wilcox-greater.tab",\%pvalues_dS_coat ); 
showHashTable( \%pvalues_dS_essen); 
_print_pvalues_dataframe_for_R( \"pmat.dS-essen-nonCE-wilcox-greater.tab",\%pvalues_dS_essen ); 

exit;
##############
# END #
#############


###############################
# Usage: $pvalue =  _pairwise.test (\@first, \@second, ERROR );
sub _pairwise_test {
  use strict; use warnings; my $debug = 0;
  my ($ref_array1, $ref_array2, $fh) = @_;
  my $pvalue = -1; 

    open (FIRST, ">/tmp/_first.070310.tab");
    foreach my $first (@$ref_array1) { print FIRST "$first\n"; }
    close (FIRST);

    open (SECOND, ">/tmp/_second.070310.tab");
    foreach my $second (@$ref_array2) { print SECOND "$second\n"; }
    close (SECOND);

    _print_R2(); 
    system("R --no-save < /tmp/_062410.2.R > /tmp/_output.062410.2.txt ");

    open (RTXT, "</tmp/_output.062410.2.txt"); my @Rlines = <RTXT>; close (RTXT); 
    # t = 1.1112, df = 40.16, p-value = 0.1365
    my @testlines = grep( /p-value =| p-value </, @Rlines);
    #my ( $ttmp, $dftmp, $ptmp ) = split( /\,/, $testlines[0] );
    #my ($txt, $pvalue) = split( /=/, $ptmp);
    if (! defined $testlines[0] ) { 
         die; 
    }
    my @tokens = split( /p-value =| p-value </, $testlines[0] );
    $pvalue = $tokens[1];  
    $pvalue =~ s/\s+|\n//g;

    print $fh  "$testlines[0]";
    print $fh "\n@Rlines\n"; 

  return $pvalue; 
}


############
sub _print_R2 { 
 my $lone_line = "
 firsttb = read.table('/tmp/_first.070310.tab');
 secondtb = read.table('/tmp/_second.070310.tab');
 summary(firsttb);
 summary(secondtb);
 wilcox.test( firsttb[,1], secondtb[,1], alter='gr'); 
 ";
 #one sided wilcox test

 open (OUT, ">/tmp/_062410.2.R");
 print OUT $lone_line; 
 close(OUT);
}

############
sub _print_R { 
 my $lone_line = "
 ctb = read.table('/tmp/coatW.tab');
 nonCEtb = read.table('/tmp/nonCEW.tab');
 summary(ctb);
 summary(nonCEtb);
 #t.test( ctb[,1], nonCEtb[,1], alter='gr'); 
 wilcox.test( ctb[,1], nonCEtb[,1], alter='gr'); 
 ";
 #one sided wilcox test

 open (OUT, ">/tmp/_062410.R");
 print OUT $lone_line; 
 close(OUT);
}

############################################################
# Usage: $flag =  _print_pvalues_dataframe_for_R(\$filename, \%pvalues);

sub _print_pvalues_dataframe_for_R {
  use strict; use warnings; my $debug = 0;
  my ($ref_fl, $ref_h) = @_;
  my %pvalues = %$ref_h; 

 # This is the order of genomes in the heatmap matrix
 my $_longline = "B.anthracis.Ames B.cereus.ATCC14579  B.cereus.E33L B.cereus.ATCC10987 B.thuringiensis  B.weihenstephanensis  B.subtilis168 B.amyloliquefaciens  B.licheniformis  B.pumilus  B.claussi  B.halodurans";
 # B.mojavensis    G.kaustophilus  O.iheyensisHTE831

 my @ordered_strains = split( /\s+/, $_longline ); 

 #open (PMAT, ">pmat-coat-nonCE-wilcox-greater.tab");
 open (PMAT, ">$$ref_fl");
  # print with row names for R
  #header
  print PMAT "\t"; foreach my $ss (@ordered_strains) { print PMAT $ss."\t"; }
  print PMAT "\n";

  #body
  foreach my $row (@ordered_strains) {
   print PMAT $row; #This is row name in R
   foreach my $column (@ordered_strains ) {
     my ($bigger, $smaller) = order_big2small( $row, $column ); 
     if( $bigger eq $smaller ) {
       print PMAT "\tNA"; 
     } elsif ( ! (defined $pvalues{"$bigger, $smaller"}) ) {
       print PMAT "\tNA";
     } else {
       print PMAT "\t". $pvalues{"$bigger, $smaller"} 
     }
   }
   print PMAT "\n";
 } 
 close(PMAT);
 return 1; # for success
}
