# Proteomic Aging Clock (PAC)
## Publication
Kuo CL, Chen Z, Liu P, Pilling LC, Atkins JL, Fortinsky RH, Kuchel GA, Diniz BS. Proteomic aging clock (PAC) predicts age-related outcomes in middle-aged and older adults. Aging Cell. 2024 May 15:e14195. doi: 10.1111/acel.14195. Epub ahead of print. PMID: 38747160.

Abstract
Beyond mere prognostication, optimal biomarkers of aging provide insights into qualitative and quantitative features of biological aging and might, therefore, offer useful information for the testing and, ultimately, clinical use of gerotherapeutics. We aimed to develop a proteomic aging clock (PAC) for all-cause mortality risk as a proxy of biological age. Data were from the UK Biobank Pharma Proteomics Project, including 53,021 participants aged between 39 and 70 years and 2923 plasma proteins assessed using the Olink Explore 3072 assay®. 10.9% of the participants died during a mean follow-up of 13.3 years, with the mean age at death of 70.1 years. The Spearman correlation between PAC proteomic age and chronological age was 0.77. PAC showed robust age-adjusted associations and predictions for all-cause mortality and the onset of various diseases in general and disease-free participants. The proteins associated with PAC proteomic age deviation were enriched in several processes related to the hallmarks of biological aging. Our results expand previous findings by showing that biological age acceleration, based on PAC, strongly predicts all-cause mortality and several incident disease outcomes. Particularly, it facilitates the evaluation of risk for multiple conditions in a disease-free population, thereby, contributing to the prevention of initial diseases, which vary among individuals and may subsequently lead to additional comorbidities.

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

