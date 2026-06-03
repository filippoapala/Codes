library(readxl)
library(ggplot2)
library(scatterplot3d)
library(plotly)
library(datasets)
library(cluster)
library(factoextra)
library(GPArotation)
library(dendextend)
library(NbClust)
library(dplyr)
library(reshape2)

#load data
data= read_xlsx("Wholesale_customers_data.xlsx")

#delete categorical variables
d<-data[,!colnames(data)%in% c("Region","Channel")]

#take logs
ld<-log(d)

#print correlation matrix
corLd<-cor(ld)
melted_cor_matrix <- melt(corLd)

ggplot(data = melted_cor_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  coord_fixed() +
  labs(title = "Correlation Heatmap", x = "Variables", y = "Variables")


#PCA 
pca <- princomp(ld,cor = TRUE)
screeplot<-(pca$sdev)^2/(sum((pca$sdev)^2))
eigen<-(pca$sdev)^2
plot(eigen)
plot(screeplot)

#loadings
loadings<-pca$loadings
print(loadings)

#get scores from PCA
scores <- as.data.frame(pca$scores)


