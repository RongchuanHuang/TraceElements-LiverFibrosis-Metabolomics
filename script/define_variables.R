library(tidyverse)
library(bruceR)
covariates <- c("Age", "Gender", "BMI", "Smoke", "Drink", "Exercises", 
                "Education_new", "diet_score", "Hypertension_new", "Hyperlipidemia_new",
                "Diabetes_new", "eGFR")
exposures_name <- c("Fe", "Zn", "Cu", "Mn", "Se", "V", "Cr", "Ni", "Mo", "Co")
exposures_log_scaled <- paste0(exposures_name, "_log_scaled")


data4 <- import("data/row/四季血浆代谢组数据712人次_250414.sav")
metabolites_name <- names(data4)[-c(1:2)]

metabolites_log_scaled <- paste0(names(data4)[-c(1:2)],"_log_scaled")

