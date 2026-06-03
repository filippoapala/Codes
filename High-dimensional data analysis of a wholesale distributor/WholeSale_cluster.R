
#Cluster analysis generic
ssd<-ld
distance=dist(ssd)
wss<- numeric(0)

#k-means cluster analysis, hirerechical clustering
wss<- numeric(0)
set.seed(2)
kc.2<-kmeans(ssd,2)
kc.3<-kmeans(ssd,3)
kc.4<-kmeans(ssd,4)
kc.5<-kmeans(ssd,5)
kc.6<-kmeans(ssd,6)
kc.7<-kmeans(ssd,7)

for (i in 2:20) {
  # Compute kmeans for each number of clusters and extract withinss
  wss[i] <- sum(kmeans(ssd, centers = i)$withinss)
}
# Plot the within-group SS against the number of clusters
plot(1:20, wss, type = "b", xlab = "Number of Clusters", ylab = "Within-group SS")


#plot k-means
fviz_cluster(kc.2, data = scale(ssd), palette = c("orange", "purple","#000000","grey"), 
             geom = "point", 
             ellipse.type = "convex", 
             ggtheme = theme_bw()
             
)


