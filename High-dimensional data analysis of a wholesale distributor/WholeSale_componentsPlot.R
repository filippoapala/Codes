
# 3D plot
#get custom colors
colors1 <- ifelse(data$Channel == "1", "purple",
                  ifelse(data$Channel == "2", "orange" , "red"))

#merge data
merged_data <- cbind(scores, data, colors1)

#2D plot  with different channels
Dplot1 <- plot_ly(data = merged_data, 
                  x = ~scores[,1], y = ~scores[,2],
                  color = ~Channel,
                  colors = colors1,
                  type = "scatter", mode = "markers")

# Display the  plots
print(Dplot1)

#get custom colors
merged_data2 <- cbind(scores, data)

#2D plot  with different regions
Dplot2 <- plot_ly(data = merged_data2, 
                  x = ~scores[,1], y = ~scores[,2],
                  color = ~Region,
                  type = "scatter", mode = "markers")

# Display the  plots
print(Dplot2)