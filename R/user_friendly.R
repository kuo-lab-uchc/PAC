#import the data and R file
data=read.csv("example_data.csv",row.names = 1)
source("proteomic_age.R")
proteomic_age(data)
