# Proteomic Aging Clock (PAC)
## Publication
Kuo CL, Chen Z, Liu P, Pilling LC, Atkins JL, Fortinsky RH, Kuchel GA, Diniz BS. Proteomic aging clock (PAC) predicts age-related outcomes in middle-aged and older adults. Aging Cell. 2024 May 15:e14195. doi: 10.1111/acel.14195. Epub ahead of print. PMID: 38747160.

## Setup
Users are required to download and install R but no R package is needed to calculate the PAC proteomic age.

## Input file
The input file, as shown in "pac_example_data.csv", should have the first column designated as the ID column, with any arbitrary column name. This should be followed by columns containing age and the required proteins. The input file may include additional proteins beyond those needed. The protein column names will be converted to lowercase by the R code "pac_proteomic_age.R".

## Example
The R code below shows you how to load the "pac_proteomic_age" function in "pac_proteomic_age.R" to calculate PAC proteomic ages for five subjects with input data in "pac_example_data.csv".

```
source("pac_proteomic_age.R")
pac_input=read.csv("pac_example_data.csv")
pac_proteomic_age(pac_input)

```

