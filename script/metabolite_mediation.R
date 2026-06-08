library(tidyverse)
library(bruceR)

# 导入数据 ####
data7 <- import("data/processed/data7.rds")
MWAS_Zn_V2 <- import("result/draft/MWAS_Zn_V2.xlsx")
MWAS_LSM_V2 <- import("result/draft/MWAS_LSM_V2.xlsx")

MWAS_Zn_V2_sig <- MWAS_Zn_V2 %>% 
  filter(PFDR < 0.2) %>% 
  dplyr::select(all_of(c("Estimate", "Medi")))

MWAS_LSM_V2_sig <- MWAS_LSM_V2 %>% 
  filter(PFDR < 0.2) %>% 
  dplyr::select(all_of(c("Estimate", "Medi")))

overlap_metabolite <- inner_join(MWAS_Zn_V2_sig,MWAS_LSM_V2_sig,by = "Medi") %>% 
  mutate(Direction = Estimate.x * Estimate.y) %>% 
  filter(Direction < 0)

## CIT ####
library(cit)

CIT_table <- tibble()

for (medi in overlap_metabolite$Medi) {
  L_data = data7[["Zn_log_scaled"]]
  G_data = data7[[medi]]  # "M1" is the abbreviation of metabolite/lipid
  T_data = data7[["liverhard_log"]]
  C_data = data7[, covariates]
  
  model <- cit.cp(L_data, G_data, T_data, C_data, rseed = 2025)     
  df_table <- as.data.frame(t(as.data.frame(model))) %>% mutate(medi_name = medi)
  CIT_table <- rbind(CIT_table, df_table) 
}

CIT_table <- CIT_table %>% mutate(PFDR = p.adjust(p_cit, "BH"))
export(CIT_table, file = "result/draft/CIT_table.xlsx")



## 传统中介 ####
CIT_table <- import("result/draft/CIT_table.xlsx")
CIT_table_sig <- CIT_table %>% filter(PFDR < 0.05)
library(mediation)
set.seed(2024)
mediation_table_metabolite_CIT_sig <- tibble()
for (medi in CIT_table_sig$medi_name) {
  Metabolite_ID_new <- medi
  
  model1 <- lmer(data = data7, formula = get(Metabolite_ID_new) ~ Zn_log_scaled + Age + Gender + BMI + Smoke + Drink + Exercises + Education + diet_score + Hypertension_new + Hyperlipidemia_new + Diabetes_new + eGFR + (1 | ID))
  model2 <- lmer(data = data7, formula = liverhard_log ~ get(Metabolite_ID_new) + Zn_log_scaled + Age + Gender + BMI + Smoke + Drink + Exercises + Education + diet_score + Hypertension_new + Hyperlipidemia_new + Diabetes_new + eGFR + (1 | ID))
  medi_model <- mediate(model1, model2, treat = "Zn_log_scaled", mediator = "get(Metabolite_ID_new)", sims = 2000)
  df <- summary(medi_model)
  df_table <- tibble(
    name = c("ACME", "ADE", "Total Effect", "Prop. Mediated"),
    Estimate = c(df$d0, df$z0, df$tau.coef, df$n0),
    `95% CI Lower` = c(df$d0.ci[1], df$z0.ci[1], df$tau.ci[1], df$n0.ci[1]),
    `95% CI Upper` = c(df$d0.ci[2], df$z0.ci[2], df$tau.ci[2], df$n0.ci[2]),
    `p-value` = c(df$d0.p, df$z0.p, df$tau.p, df$n0.p)
  ) %>% mutate(
    exposure_name = "Zn",
    medi_name = medi,
    outcome_name = "LSM"
    
  )
  mediation_table_metabolite_CIT_sig <- rbind(mediation_table_metabolite_CIT_sig, df_table)
}

mediation_table_metabolite_CIT_sig_V2 <- mediation_table_metabolite_CIT_sig %>% process_mediation_results()
export(mediation_table_metabolite_CIT_sig_V2,"result/draft/mediation_table_metabolite_CIT_sig_V2.xlsx")

## 合并CIT 和 传统中介 ####
Table_3 <- CIT_table %>% 
  filter(PFDR<0.05) %>% 
  left_join(mediation_table_metabolite_CIT_sig_V2, by = "medi_name") %>% 
  mutate(Metabolite_ID = str_split_i(medi_name, "_log", 1)) %>% 
  left_join(code_book, by = "Metabolite_ID") %>% 
  dplyr::select(all_of(c("exposure_name", "Metabolite_Name", "outcome_name", "p_cit", "PFDR", "ACME (β, 95%CI)",
                         "p-value_ACME","Prop. Mediated","p-value_Prop. Mediated")))

export(Table_3, file = "result/Table/Table 3.xlsx")


