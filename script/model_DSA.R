library(tidyverse)
library(bruceR)

data5 <- import("data/processed/data5.rds")
source("script/varibles_func.R")

summary(data5[, c(outcome_log, exposure_log_scale, covariate)])

# 加载必要的包
library(DSA)

set.seed(2026) # 设置种子以便于结果复现
# 生成一个没有重复的随机数列，范围在 1 到 1000 之间，选择 100 个数
random_numbers <- sample(1:100000, 1000, replace = FALSE)
export(as.tibble(random_numbers), file = "result/draft/random_numbers.txt")


DSA_model <- list()

for (i in random_numbers) {
  ID <- as.integer(as.factor(data5[, "ID"]))
  data_test <- data5[, c(exposure_log_scale, "AST_log", covariate)]
  formula <- as.formula(paste("AST_log", "~ 1 + ", paste0(covariate, collapse = "+")))
  DSA_model[[i]] <- DSA(
    formula = formula, data = data_test, maxsize = 21,
    maxorderint = 1, maxsumofpow = 1, vfold = 10,
    family = "gaussian", userseed = i, id = ID
  )
}

save(DSA_model, file = "result/model/DSA_model.Rdata")
# 计算列表中不是 NULL 的元素个数
load("result/model/DSA_model.Rdata")
random_numbers <- import("result/draft/random_numbers.txt")
all_vars <- c()
for (i in random_numbers$value) {
  model_df <- DSA_model[i]

  model_df <- unlist(model_df[[1]]$model.selected) |> as.character()
  vars <- gsub("I\\(|\\^1\\)", "", strsplit(model_df[3], " \\+ ")[[1]])
  all_vars <- c(all_vars, vars)
}
freq_table <- as.tibble(table(all_vars)) %>%
  arrange(n) %>%
  filter(n < 1000) %>% 
  mutate(Exposure = str_split_i(all_vars,"_log",i = 1)) 

export(freq_table, file = "result/draft/DSA_table.xlsx")
