library(ape);
   postscript("_tree.071808.ps",width=8,height=8,horizontal=F);
   t=read.tree("protbaci.label.nwk");
   plot(t, main="BG10499-BG10498 final");
   dev.off();
  