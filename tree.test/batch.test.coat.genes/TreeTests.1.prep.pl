#Tree tests
use strict; use warnings; 
#use lib '/home/hqin/lib/perl'; use Util; use FASTAshort; 

system("ls -d BG*/Bha*/H0 > /tmp/_bg.txt"); 

open (IN,"</tmp/_bg.txt"); my @lines = <IN>; close (IN);

foreach my $line (@lines) {
 chomp $line; 
 system("cp trees.nwk $line/.");
 system("cp codeml.ctl $line/.");
}
}


exit; 


#
# END
#
