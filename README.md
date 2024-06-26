# Proteomic Aging Clock 
## Publication
Kuo CL, Chen Z, Liu P, Pilling LC, Atkins JL, Fortinsky RH, Kuchel GA, Diniz BS. Proteomic aging clock (PAC) predicts age-related outcomes in middle-aged and older adults. Aging Cell. 2024 May 15:e14195. doi: 10.1111/acel.14195. Epub ahead of print. PMID: 38747160.

## Setup
Users are required to download and install R but no R package is needed to calculate the PAC proteomic age.
## Examples
An example dataset (pac_example_data.csv") can be found in the data folder, including five individuals with chronological age and 128 proteins in PAC. The R code below shows you how to load the function in "pac_proteomic_age.R" to calculate the PAC proteomic ages using the input data from the five individuals.

```
source("pac_proteomic_age.R")
pac_input=read.csv("pac_example_data.csv")
pac_proteomic_age(pac_input)

```

