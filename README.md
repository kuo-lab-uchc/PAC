# Proteomic Aging clock 
## Setup
Users are required to download and install R but no R package is needed to calculate the PAC proteomic age.

## Examples
An example dataset can be found in the data folder, including five individuals with chronological age and 128 proteins in PAC. The R code below shows you how to load the function in "pac_proteomic_age.R" to calculate the PAC proteomic ages using the input data from the five individuals in "example.data.csv".

```
source("proteomic_age.R")
pac_input=read.csv("example_data.csv",row.names = 1)
pac_proteomic_age(pac_input)

```

