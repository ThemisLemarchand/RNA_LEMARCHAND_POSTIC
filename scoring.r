#LEMARCHAND Thémis
#20222535 M2GENIOMHE

library(ggplot2)

pairs_list = c("AA","AU","AC","AG","UU","UC","UG","CC","CG","GG")

#import all the data from files into data frames
tables_list <- list()

for (pair in pairs_list) {
  # Construct the file path
  file_path <- file.path(getwd(), pair)  # Adjust the extension and path as needed
  
  # Read the table from the file
  table_data <- read.csv(file_path)  # Adjust read function based on file format
  
  # Assign the table to a list element with the file name as the key
  tables_list[[pair]] <- table_data
}

# Plot the data for each interaction profile

for (pair in pairs_list) {
  data <- data.frame(
  distance = c(1:19),
  score = tables_list[[pair]]$X0
  )
  name=paste(pair,".pdf",sep="")
  pdf(name)

  ggplot(data, aes(x = distance, y = score)) +  geom_line() +  geom_point() +
  labs(title = "Interaction Profile",
       x = "Distance",
       y = "Score")
  
  dev.off()
}


