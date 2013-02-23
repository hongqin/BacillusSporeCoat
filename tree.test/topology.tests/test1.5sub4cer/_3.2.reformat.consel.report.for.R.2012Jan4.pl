#Jan 4,2013
use strict; use warnings; 
use lib '/Users/hongqin/lib/perl'; use Util; 

my $debug= 2 ;

my $inclusterfile = '_test.csv';
my $sourcedir = '/Users/hongqin/projects/coat.protein07/ortholog.analysis/essetial.mcl.rbh.090208/bacillus.phylip.I3.1/wkdir';

open (IN,"<$inclusterfile"); my @lines = <IN>; close (IN); shift @lines;
my @clusterdirs = ();
foreach my $line (@lines) {
  chomp $line;
  my ($num, $bg, @res) = split( /\t/, $line);
  print "[$num][$bg]------";
  push (@clusterdirs, $bg);
}

open ( OUT, ">_conselreport.34essenGene.2012Jan4.tab" ); 
print OUT "BG\trank\titem\tobs\tau\tnp\tbp\tpp\tkh\tsh\twkh\twsh\n"; 

my $count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  if ( !( -e "wkdir/$mydir/$mydir.consel.report.txt" ) ) { 
    print "Error: $mydir has no consel report"; 
  } else {
    $count ++; print "Cluster:$count I am now working on [$mydir]::\n"; 
  	open (IN, "<wkdir/$mydir/$mydir.consel.report.txt"); @lines = <IN>; close (IN);
  	shift @lines; shift @lines; shift @lines;  #remove the first three lines
  	pop @lines; #remove the last empty line
  	foreach my $line (@lines) {
  	   $line =~ s/^#\s+//o;
  	   $line =~ s/\|//g; 
  	   chomp $line;
  	   my @els = split( /\s+/, $line);
  	   print OUT $mydir;
  	   foreach my $el (@els) { print OUT "\t$el"; }
  	   print OUT "\n";     
  	}
    #chdir "wkdir/$mydir";
    # chdir "../..";
  }
}

close (OUT);

exit;

#
# END
#

########################