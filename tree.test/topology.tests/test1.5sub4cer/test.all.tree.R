list.files()

require(ape)
help(package=ape)

#read.tree("trees1.nwk")

trees = read.tree("treesall.nwk") #list of 10
trees[[1]]
trees[[2]]
dist.topo( trees[[1]], trees[[10]])

d.tree = matrix( nrow=10, ncol=10 )
for ( i in 1:10 ) {
  for (j in i:10) {
    print ( c(i, j))
    d.tree[i,j] = dist.topo( trees[[i]], trees[[j]]  )
  }
}

d.tree


