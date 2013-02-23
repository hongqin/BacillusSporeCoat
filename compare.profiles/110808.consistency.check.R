
 rm(list=ls());

 ### read the 110608 profile 
 ctb = read.table( "_coat.profile.110608.csv",sep=",",header=T);
 ctb$id = as.character( ctb$id );
 row.names( ctb ) = as.character( ctb$id );
 ctb = ctb[match(sort(ctb$id), ctb$id), ]


 ctb2 = ctb[, -c(2:9, 21)];
 names(ctb2) =  c( "NCBI", names(ctb2)[2:12] );
 bacillus.specs = names(ctb)[10:20];

 names( ctb2 ) [c(2,5,7,8,9,11,12)] #these are Henriques07 columns

 ctb4 = ctb2[, c(2,5,7,8,9,11,12)] 
 #row.names(ctb4) = as.character( ctb4$Bsu );
 ctb4 = ifelse( is.na(ctb4), 0, 1);
 ctb4 = data.frame(ctb4);
 head(ctb4); #passed

###### read the old profile
 
old = read.csv("/home/hqin/coat.protein07/key.data/_coat-phylo-profiles.122107.csv",
 sep="\t");
old$Bsu = as.character(old$Bsu);

old = old[, c(1,3,11,17)]
old2 = old;
old2b = old2; 

for( c in 2:4 ) {
 old2[,c] = as.numeric( old2[,c] );
}

assign = function( id, short ) {
  ret = 0;
  if ( substr(as.character(id), 1, 2) == short ) {
    ret = 1;
  } 
  ret; 
}
shorts = c('BA','BC','NT');
for(r in 1:73) {
  for( c in 2:4) {
    old2[r,c] = assign(old2b[r,c], shorts[c-1])
  }
}

old3 = old2[match(sort(old2$Bsu), old2$Bsu), ];

#######
#### compare new and old profiles
ctb5 = ctb4[,c(3,4,5)]

names(old3) = c( "Bsu", names(ctb5) );

row.names(ctb5) %in% old3$Bsu;

plot( match( row.names(ctb5), old3$Bsu ) );

d= numeric(length(ctb5[,1]))
names(d) = row.names(ctb5);


for( r in 1:length(ctb5[,1]) ){  
 d[r] = dist( rbind(ctb5[r,], old3[r,2:4]) ); #Euclidean distance
}
d = d^2;

changed = d[d>0]; #15 changed genes. 
 write.table(changed,"_coat.changed.cereus.profile.110808.csv",quote=F,col.name=F);

x = c();
y = c();
for( r in 1:length(ctb5[,1]) ){  
 if ( ctb5[r,1] < old3[r,2] ) { x = c(x, old3$Bsu[r]); }
 if ( ctb5[r,1] > old3[r,2] ) { y = c(y, old3$Bsu[r]); }
}

old[match(x,old$Bsu),]
ctb5[x,]
ctb[x,]




### 




quit('no');

