
  library(ape);
  files = read.table( "__tree.txt")
  files = files[,1];

  pdf ("_prot.nj.plot.070508.pdf", width=10,height=8);
  par(mar=c(2,2,2,2)+0.1);

  for( i in 1:length(files) ) {
     file = as.character(files[i]);
     tr = read.nexus( file  );
     plot( tr, main=paste("BG10608",file), show.node.label=T );
     plot( tr, main=paste("BG10608",file), type="u" );
     plot( tr, main=paste("BG10608",file), type="c", show.node.label=T  );
     i;
  }

  tnj = read.nexus( as.character(files[1]) );
  tbs = read.nexus( as.character(files[2]) );

  d = dist.topo( tnj, tbs);
  if( d[1]== 0 ) { 
    write.table( d, "__nj.unchanged.in.bootstrap.070508.txt");
  } else {
    write.table( d, "__nj.changed.in.bootstrap.070508.txt");
  } 

  dev.off();
  