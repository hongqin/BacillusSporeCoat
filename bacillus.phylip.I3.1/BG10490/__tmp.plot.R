library(ape);
   postscript("_tree.073108.ps",width=8,height=8,horizontal=F);
   t=read.tree("protbaci.label.curated.nwk");
   plot(t, main="BG10490 final");
   dev.off();
  