library(bruceR)
library(tidyverse)

data5 <- import("data/processed/data5.rds")

source("func/functions.R")
source("script/define_variables.R")

# 协变量转换 ####
data6 <- data5 |>
  mutate(
    Age = as.numeric(age),
    Gender = Gender,
    BMI = weight / (height / 100)^2,
    Smoke = Smoke,
    Drink = Drink,
    Exercises = exercises,
    Education_new = case_when(
      Education == 1 ~ 1,
      Education %in% c(2, 3) ~ 2,
      TRUE ~ 3
    ),
    
    ## 构建健康饮食评分
    Fruits = diet_function(Orange) + diet_function(Applepear) + diet_function(Banana) + diet_function(Berry) + diet_function(Mango) + diet_function(Pitaya) + diet_function(Watermelon) + diet_function(Hamimelon) + diet_function(Litchi) + diet_function(Jujube) + diet_function(Peach),
    Vegetables = diet_function(Darkvegetable) + diet_function(Lightvegetable) + diet_function(Rootvegetables) + diet_function(Radish) + diet_function(Carrot) + diet_function(Freshbean) + diet_function(Eggplant) + diet_function(Tomato) + diet_function(Chili) + diet_function(Mushroom) + diet_function(Agaric) + diet_function(Algae),
    Meat = diet_function(Meat) + diet_function(Beef) + diet_function(Mutton),
    Bean = diet_function(Bean) + diet_function(Mixedbeans) + diet_function(Driedslag) + diet_function(Slagsoup),
    Fish = diet_function(Freshwater) + diet_function(Saltwater) + diet_function(Crab) + diet_function(Squid),
    Fruits = ifelse(Fruits >= 7, 1, 0),
    Vegetables = ifelse(Vegetables >= 7, 1, 0),
    Meat = ifelse(Meat %in% c(1, 2, 3, 4, 5, 6), 1, 0), ## 每周1-6
    Bean = ifelse(Bean >= 4, 1, 0), ## 每周大于4
    Fish = ifelse(Fish >= 1, 1, 0), ## 每周大于1
    
    diet_score = Fruits + Vegetables + Meat + Bean + Fish,
    
    ## 疾病定义
    
    Hypertension_new = case_when(
      Hypertension == 1 | SBP >= 140 | DBP >= 90 | drug_hypertension == 1 ~ 1,
      TRUE ~ 0
    ),
    Hyperlipidemia_new = case_when(
      Hyperlipidemia == 1 | CHOL > 5.72 | TG > 1.7 | drug_lipid == 1 ~ 1,
      TRUE ~ 0
    ),
    Diabetes_new = case_when(
      Diabetes == 1 | GLU >= 7 | HbA1c >= 6.5 | drug_jiangtang == 1 | drug_insulin == 1 ~ 1,
      TRUE ~ 0
    ),
    Abdominal_obesity = case_when(
      (gender == 1 & (waist >= 90 | WHR >= 0.9)) | (gender == 0 & (waist >= 85 | WHR >= 0.85)) ~ 1,
      TRUE ~ 0
    ),
    
    # CKD-EPI核心计算
    CREA_mg_dl = CREA / 88.4,
    eGFR = case_when(
      # 男性分支
      Gender == "1" & CREA_mg_dl <= 0.9 ~ 141 * (CREA_mg_dl / 0.9)^(-0.411) * 0.993^Age,
      Gender == "1" & CREA_mg_dl > 0.9 ~ 141 * (CREA_mg_dl / 0.9)^(-1.209) * 0.993^Age,
      
      # 女性分支
      Gender == "2" & CREA_mg_dl <= 0.7 ~ 144 * (CREA_mg_dl / 0.7)^(-0.329) * 0.993^Age,
      Gender == "2" & CREA_mg_dl > 0.7 ~ 144 * (CREA_mg_dl / 0.7)^(-1.209) * 0.993^Age,
      
      # 意外值处理
      TRUE ~ NA_real_
    ),
    
    ## 结局定义
    liverhard_degree = case_when(
      liverhard >= 7.3 ~ 1,
      is.na(liverhard) ~ NA, # 将NA单独编码为NA
      TRUE ~ 0 # 其他情况（包括 <7.3 和未定义值）编码为0
    ),
    fat_decay_degree = case_when(
      fat_decay > 240 ~ 1,
      TRUE ~ 0
    )
  )





# 变量因子化 ####
data6 <- data6 |>
  mutate(across(
    .cols = all_of(c("Age", "BMI", "diet_score", "liverhard", "fat_decay", "eGFR")), # 选择列
    .fns = ~ as.numeric(.x) # 转换函数
  ))
data6 <- data6 |>
  mutate(across(
    .cols = all_of(c("Gender", "Smoke", "Drink", "Exercises", "Education_new", "Hypertension_new", "Hyperlipidemia_new", "Diabetes_new", "liverhard_degree", "fat_decay_degree", "ID")), # 选择列
    .fns = ~ as.factor(.x) # 转换函数
  ))

# 多重插补 ####
library(mice)
source("func/functions.R")
df <- as.data.frame(data6[, covariates])
imputed_data <- mice(df, m = 10, seed = 123)


completed_data_list <- map(1:10, ~ complete(imputed_data, .))

Drink_table1 <- as.data.frame(completed_data_list[[1]]) |> dplyr::select(Drink)
Drink_table2 <- as.data.frame(completed_data_list[[2]]) |> dplyr::select(Drink)
Drink_table3 <- as.data.frame(completed_data_list[[3]]) |> dplyr::select(Drink)
Drink_table4 <- as.data.frame(completed_data_list[[4]]) |> dplyr::select(Drink)
Drink_table5 <- as.data.frame(completed_data_list[[5]]) |> dplyr::select(Drink)
Drink_table6 <- as.data.frame(completed_data_list[[6]]) |> dplyr::select(Drink)
Drink_table7 <- as.data.frame(completed_data_list[[7]]) |> dplyr::select(Drink)
Drink_table8 <- as.data.frame(completed_data_list[[8]]) |> dplyr::select(Drink)
Drink_table9 <- as.data.frame(completed_data_list[[9]]) |> dplyr::select(Drink)
Drink_table10 <- as.data.frame(completed_data_list[[10]]) |> dplyr::select(Drink)

drink_combined <- bind_cols(Drink_table1, Drink_table2, Drink_table3, Drink_table4, Drink_table5, Drink_table6, Drink_table7, Drink_table8, Drink_table9, Drink_table10)


# 添加众数列
drink_combined <- drink_combined %>%
  rowwise() %>%
  mutate(
    Drink = get_mode(c_across(Drink...1:Drink...10))
  ) %>%
  ungroup()

Exercises_table1 <- as.data.frame(completed_data_list[[1]]) |> dplyr::select(Exercises)
Exercises_table2 <- as.data.frame(completed_data_list[[2]]) |> dplyr::select(Exercises)
Exercises_table3 <- as.data.frame(completed_data_list[[3]]) |> dplyr::select(Exercises)
Exercises_table4 <- as.data.frame(completed_data_list[[4]]) |> dplyr::select(Exercises)
Exercises_table5 <- as.data.frame(completed_data_list[[5]]) |> dplyr::select(Exercises)
Exercises_table6 <- as.data.frame(completed_data_list[[6]]) |> dplyr::select(Exercises)
Exercises_table7 <- as.data.frame(completed_data_list[[7]]) |> dplyr::select(Exercises)
Exercises_table8 <- as.data.frame(completed_data_list[[8]]) |> dplyr::select(Exercises)
Exercises_table9 <- as.data.frame(completed_data_list[[9]]) |> dplyr::select(Exercises)
Exercises_table10 <- as.data.frame(completed_data_list[[10]]) |> dplyr::select(Exercises)

Exercises_combined <- bind_cols(Exercises_table1, Exercises_table2, Exercises_table3, Exercises_table4, Exercises_table5, Exercises_table6, Exercises_table7, Exercises_table8, Exercises_table9, Exercises_table10)

# 添加众数列
Exercises_combined <- Exercises_combined %>%
  rowwise() %>%
  mutate(
    Exercises = get_mode(c_across(Exercises...1:Exercises...10))
  ) %>%
  ungroup()

eGFR_table1 <- as.data.frame(completed_data_list[[1]]) |> dplyr::select(eGFR)
eGFR_table2 <- as.data.frame(completed_data_list[[2]]) |> dplyr::select(eGFR)
eGFR_table3 <- as.data.frame(completed_data_list[[3]]) |> dplyr::select(eGFR)
eGFR_table4 <- as.data.frame(completed_data_list[[4]]) |> dplyr::select(eGFR)
eGFR_table5 <- as.data.frame(completed_data_list[[5]]) |> dplyr::select(eGFR)
eGFR_table6 <- as.data.frame(completed_data_list[[6]]) |> dplyr::select(eGFR)
eGFR_table7 <- as.data.frame(completed_data_list[[7]]) |> dplyr::select(eGFR)
eGFR_table8 <- as.data.frame(completed_data_list[[8]]) |> dplyr::select(eGFR)
eGFR_table9 <- as.data.frame(completed_data_list[[9]]) |> dplyr::select(eGFR)
eGFR_table10 <- as.data.frame(completed_data_list[[10]]) |> dplyr::select(eGFR)

eGFR_combined <- bind_cols(eGFR_table1, eGFR_table2, eGFR_table3, eGFR_table4, eGFR_table5, eGFR_table6, eGFR_table7, eGFR_table8, eGFR_table9, eGFR_table10)
eGFR_combined <- eGFR_combined %>%
  mutate(
    eGFR_mean = rowMeans(across(starts_with("eGFR...")), na.rm = TRUE)
  )

data6$Drink <- drink_combined$Drink
data6$Exercises <- Exercises_combined$Exercises
data6$eGFR <- eGFR_combined$eGFR_mean


# 排除缺失 ####
nrow(data6)
data6 <- data6 |> filter(!is.na(liverhard))
nrow(data6)
summary(as.factor(data6$season.x))
data6 <- data6 |> filter(!is.na(Mg))
nrow(data6)
data6 <- data6 |> filter(!is.na(M110T27))
nrow(data6)


# Fe,Cu 用1/2的最小值替代 ####
Fe_min <- min(data6$Fe, na.rm = T)
Cu_min <- min(data6$Cu, na.rm = T)

data6$Fe <- ifelse(is.na(data6$Fe), Fe_min / 2, data6$Fe)
data6$Cu <- ifelse(is.na(data6$Cu), Cu_min / 2, data6$Cu)

sum(is.na(data6$Cu))




# 输出 ####
export(data6, file = "data/processed/data6.rds")
