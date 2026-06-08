library(tidyverse)
library(bruceR)

data6 <- import("data/processed/data6.rds")

# 正态化 - 暴露 ####
data6 <- data6 %>%
  # 步骤1：对变量进行log转换，Cd使用log10(x+1)，其他变量使用log10(x)
  mutate(across(
    .cols = all_of(c(exposures_name, "liverhard", "fat_decay")), # 排除Cd的其他变量
    .fns = ~ log10(.x),
    .names = "{.col}_log"
  )) %>%
  # 步骤2：仅对exposures对应的log变量进行scale标准化
  mutate(across(
    .cols = all_of(paste0(exposures_name, "_log")), # 动态构造变量名
    .fns = ~ as.numeric(scale(.x)), # 确保输出为数值向量
    .names = "{.col}_scaled" # 新变量名：原变量名_log_scaled
  ))

# 正态化 - 代谢 ####
data7 <- data6 |>
  mutate(across(
    .cols = all_of(metabolites_name),
    .fns = ~ log10(.x),
    .names = "{.col}_log"
  )) |>
  mutate(across(
    .cols = all_of(paste0(metabolites_name, "_log")),
    .fns = ~ as.numeric(scale(.x)),
    .names = "{.col}_scaled"
  ))

# 输出 ####
export(data7, file = "data/processed/data7.rds")