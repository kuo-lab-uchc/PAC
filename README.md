# Proteomic Aging clock 
## Setup
Users are required to download and install R but no R package is needed to calculate the PAC proteomic age.
## Examples
An example dataset can be found in the data folder, including 5 individuals with chronological age and 128 proteins in PAC.

```
data=read.csv("example_data.csv",row.names = 1)
source("proteomic_age.R")
proteomic_age(data)
```
And proteomic_age.R can be seen in R foler. It shows the total steps to calculate the Proteomic Aging clock. 