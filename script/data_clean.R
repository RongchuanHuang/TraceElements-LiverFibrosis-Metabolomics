library(bruceR)
library(tidyverse)

data5 <- import("data/processed/data5.rds")

source("func/functions.R")

# 协变量转换 ####
data6 <- data5 |>
  mutate(
    Age = as.numeric(age),
    Gender = Gender,
    BMI = weight / (height / 100)^2,
    Smoke = Smoke,
    Drink = Drink,
    Exercises = exercises,
    Education = case_when(
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

# 

