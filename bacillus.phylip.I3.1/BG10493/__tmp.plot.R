library(ape);
   postscript("_tree.071808.ps",width=8,height=8,horizontal=F);
   t=read.tree("protbaci.label.nwk");
   plot(t, main="BG10493 final");
   dev.off();
  