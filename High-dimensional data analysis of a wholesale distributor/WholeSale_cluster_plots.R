#PLOT DENDOGRAMS
hc.c<-hclust(distance)
plot(hc.c)

hc.s <- hclust(distance, method='single')
dends<-as.dendrogram(hc.s)
labels(dends)<-NULL
plot(hc.s)

hc.w<-hclust(distance, method='ward.D2') #ward method 
dend_colored <- color_branches(as.dendrogram(hc.w), k = 2, col = c("purple", "orange"))
labels(dend_colored) <- NULL
plot(dend_colored, main="Ward method, 2 clusters")
abline(h=28, col= 'black')

hc.ce<-hclust(distance, method='centroid') 
dendce<-as.dendrogram(hc.ce)
labels(dendce)<-NULL
plot(hc.ce)

#WARD METHOD SILHOUTTE
for (i in numeric_vector) {
  sil_w_i <- silhouette(cutree(hc.w,k=i),distance) #silhoutte method for model evaluation of the ward method (gives higher silhoutte score but one observation one cluster)
  rownames(sil_w_i) <- rownames(ssd) 
  plot(sil_w_i)
}

#NEAREST NEIGHBOUR
for (i in numeric_vector) {
  sil_s_i <- silhouette(cutree(hc.s,k=i),distance) #silhoutte method for model evaluation of the ward method (gives higher silhoutte score but one observation one cluster)
  rownames(sil_c_i) <- rownames(ssd) 
  plot(sil_s_i)
}

#FARTEST NEIGBOUR
numeric_vector <- c( 2, 3, 4, 5)
for (i in numeric_vector) {
  sil_c_i <- silhouette(cutree(hc.c,k=i),distance) #silhoutte method for model evaluation of the ward method (gives higher silhoutte score but one observation one cluster)
  rownames(sil_c_i) <- rownames(ssd) 
  plot(sil_c_i)
}


#plot k-means
fviz_cluster(kc.2, data = scale(ssd), palette = c("orange", "purple","#000000","grey"), 
             geom = "point", 
             ellipse.type = "convex", 
             ggtheme = theme_bw()
             
)
