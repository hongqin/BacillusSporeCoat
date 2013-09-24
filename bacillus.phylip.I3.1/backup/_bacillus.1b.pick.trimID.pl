# 071208 pick entries in the Bacillus clade,trim IDs to 10 letters for PHYLIP
# generate id converstion tables, 
 
use strict; use warnings; 
use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

my $debug= 2 ;

my $idfl   = "/home/hqin/projects/coat.protein07/key.data/_id2genomes.csv";
my $specfl = "/home/hqin/projects/coat.protein07/key.data/_genome2species.csv";
my $clusterfl = "infile-BGclusters.txt";

#my $tre ="((((((Bmo,Bsu),Bam),Bli),Bpu),((Bce,Ban,Bth),Bwe)),(Bha,Bcl))";
# OR in a more readble format: 
#	( 
#	(
#	 	( ( ( (Bmo,Bsu),Bam), Bli), Bpu ), 
#		( (Bce,Ban, Bth),Bwe)
#	), (Bha, Bcl) 
#	);

### find out the cluster-directories
#my $tmpfl = "/tmp/_coatBGclus.062208.txt";
#   system( "ls -d BG* | head -n 2 > $tmpfl " );
#system( "ls -d BG*  > $tmpfl " );
open (IN, "<$clusterfl" ); my @clusterdirs = <IN>; close (IN);
#system( "rm $tmpfl" );

my %id2genomes = ();
open (IN, "<$idfl" );
while (my $line = <IN> ) {
  chomp $line;
  my ($id, $g, @res) = split (/\t/, $line);
  $id2genomes{$id} = $g;
}
close (IN);

my %genome2species = ();
open (IN, "<$specfl");
my $count = 1;
while (my $line = <IN> ) {
  if ( $count >= 2 ) {  #skip the header
  chomp $line;
  my ($genome, $spec, $flag, @res) = split (/\t/, $line);
  print "[$genome] [$spec] [$flag]\t";
  $genome2species{ $genome } = { 
   'SPEC' => $spec,
   'FLAG' =>  $flag };
  } else { $count ++;
  } 
}
close (IN);

print "\n*************************\n";
foreach my $key (keys %genome2species ) {
  print $genome2species{$key}->{SPEC} ."::".$genome2species{$key}->{FLAG} . "\t\t";
}
print "\n*************************\n";

$count = 0;
foreach my $mydir ( @clusterdirs ) {
  chomp $mydir;
  chdir $mydir;
  $count ++; print "\n-------\nCluster:$count I am now working on [$mydir]::\n"; 

  FASTAshort::paddle_fasta_file2( "prot.faa", "/tmp/__tmp.$mydir.prot.faa" );

  open (IN, "</tmp/__tmp.$mydir.prot.faa"); my @lines=<IN>; close (IN);

  open ( OUT, ">prot-bacillus.spec.faa" );
  for ( my $i=0; $i<= ($#lines -1); $i=$i+2) {
     my $cur_line = $lines[$i];
     chomp $cur_line; 
     $cur_line =~ s/^>\s*//o;
     my ($id, @res) = split (/\s+/, $cur_line);
     my $myflag = $genome2species{ $id2genomes{$id} }->{FLAG};
     if ($debug) { print "$i id=[$id] flag =[$myflag]\n"; }
     if ( $myflag eq 'Y' ) {
         print OUT ">$id".'_'.$genome2species{ $id2genomes{$id} }->{SPEC}."\n";
         print OUT $lines[ $i + 1 ];
     }
     
  }
  close (OUT);
  
  my %table = ();
  open( IN,  "<prot-bacillus.spec.faa" );
  open( OUT, ">prot-bacillus.10char.faa" );
  while( my $line = <IN> ) {
    if ( $line =~ /^>/ ) {
       chomp $line; 
       $line =~ s/_...$//o;
       $line =~ s/^>//o;
       my $len = length( $line );
       my $delta = $len - 10;
       my $newline = '';
       if ( $delta >= 1 ) { # id is longer than 10
          $newline = substr( $line,  $delta, 10 );
          if($debug) { print "\$delta=$delta \$len=$len [$line]->[$newline]\n"; }
       } elsif ( $delta <=-1 ) { # id is shorter than 10
          $newline = $line; 
          foreach my $i ( 1 .. (-$delta) ) {
             $newline .= 'x';
          }
       } else { #id is 10 chars
          $newline = $line; 
       }
       $table{ $newline } = $line; 
       $line = $newline; 
       print  "[$line]\n"; 
       print OUT ">$line\n"; 
    } else {
       print OUT $line; 
    } 
  }
  close(IN);
  close(OUT);

  open (OUT, ">_ids.bacillus.csv");
  foreach my $short ( keys %table ) {
     print OUT $short."\t".$table{$short}."\n";
  }
  close (OUT);

  chdir "..";
}

exit;

#
# END
#

