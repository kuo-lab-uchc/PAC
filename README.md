# Proteomic Aging clock 
## Setup
You need to download and install R but no R package is required for the calculation of PAC.
## Examples
Example data can be found in data folder including 5 individuals with 128 proteins and chronological age. You can get PAC by running

```
data=read.csv("example_data.csv",row.names = 1)
source("proteomic_age.R")
proteomic_age(data)
```
And proteomic_age.R can be seen in R foler. It shows the total steps to calculate the Proteomic Aging clock. 
