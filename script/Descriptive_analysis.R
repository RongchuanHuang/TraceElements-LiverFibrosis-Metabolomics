library(bruceR)
library(tidyverse)


data6 <- import("data/processed/data6.rds")
source("script/define_variables.R")
# 基线表制作 ####
library(crosstable)
crosstable(data6, cols = "liverhard", by = "Gender")

table_base <- crosstable(data = data6, cols = c(covariates, "liverhard", "fat_decay"), by = "season.y", percent_digits = 2, num_digits = 2, test = T)
S1 <- crosstable(data = data6 |> filter(season.x == 1), cols = c(covariates, "liverhard", "fat_decay"), percent_digits = 2, num_digits = 2)
S2 <- crosstable(data = data6 |> filter(season.x == 2), cols = c(covariates, "liverhard", "fat_decay"), percent_digits = 2, num_digits = 2)
S3 <- crosstable(data = data6 |> filter(season.x == 3), cols = c(covariates, "liverhard", "fat_decay"), percent_digits = 2, num_digits = 2)
S4 <- crosstable(data = data6 |> filter(season.x == 4), cols = c(covariates, "liverhard", "fat_decay"), percent_digits = 2, num_digits = 2)

data_list <- list(S1, S2, S3, S4, table_base)
base_merge <- purrr::reduce(
  data_list,
  function(df1, df2) {
    left_join(df1, df2, by = c(".id", "label", "variable"))
  }
)

Table1_1 <- base_merge %>% 
  filter(label %in% c("Age","BMI","diet_score","eGFR")) %>% 
  filter(variable == "Mean (std)") %>% 
  select(all_of(c("label","variable", "value.x","value.y","value.x.x","value.y.y")))

Table1_2 <- base_merge %>% 
  filter(label %in% c("Gender","Smoke","Drink","Exercises","Education_new","Hypertension_new","Hyperlipidemia_new","Diabetes_new")) %>% 
  select(all_of(c("label","variable", "value.x","value.y","value.x.x","value.y.y")))

Table1_3 <- base_merge %>% 
  filter(label %in% c("liverhard", "fat_decay")) %>% 
  filter(variable == "Med [IQR]") %>% 
  select(all_of(c("label","variable", "value.x","value.y","value.x.x","value.y.y")))

Table1 <- bind_rows(Table1_1, Table1_2) %>% 
  mutate(label = factor(label, levels = c(covariates))) %>% 
  arrange(label) %>% 
  bind_rows(Table1_3)
  

export(Table1, file = "result/Table/Table1.xlsx")



# 暴露的分布 ####

## 浓度-使用crosstable 
exposure_table <- crosstable(data = data6 , cols = c(exposures_name), by = "season.x", num_digits = 2)

## 建模函数-LME
# analyze_longitudinal_consistency_lme <- function(data, id_var, time_var, variables) {
  
  # 确保 time 是 factor
#   data[[time_var]] <- factor(data[[time_var]])
#   
#   results <- tibble(
#     Exposure = character(),
#     p = numeric()
#   )
#   
#   for (var in variables) {
#     cat("\n分析变量:", var, "\n")
#     
#     # 构建 LME 公式
#     lme_formula <- as.formula(
#       paste(var, "~", time_var, "+ (1 |", id_var, ")")
#     )
#     
#     # 拟合模型
#     fit <- lmer(lme_formula, data = data, REML = FALSE)
#     
#     # 提取 season 的整体 p 值（Type III）
#     anova_res <- anova(fit)
#     
#     p_value <- anova_res[time_var, "Pr(>F)"]
#     
#     results <- results %>%
#       add_row(Exposure = var, p = round(p_value, 3))
#   }
#   
#   return(results)
# }

analyze_longitudinal_consistency <- function(data, id_var, time_var, variables) {
  results <- tibble(Exposure = character(), p = numeric())
  
  for (var in variables) {
    cat("\n分析变量:", var, "\n")
    
    # ANOVA分析
    anova_formula <- as.formula(paste(var, "~", time_var, "+ Error(", id_var, "/", time_var, ")"))
    anova_result <- summary(aov(anova_formula, data = data))
    
    p_value <- anova_result$`Error: ID:season`[[1]][1, "Pr(>F)"]
    
    results <- results %>% add_row(Exposure = var, p = round(p_value, 3))
  }
  
  return(results)
}

data6$season.x <- as.factor(data6$season.x)
## 拟合
# final_results <- analyze_longitudinal_consistency_lme(
#   data = data6,
#   id_var = "ID",
#   time_var = "season.x",
#   variables = exposures_name
# )

final_results <- analyze_longitudinal_consistency(data6, "ID", "season.x", exposures_name)

Table_S1 <- exposure_table %>% 
  filter(variable == "Med [IQR]") %>% 
  rename(Exposure = label) %>%
  left_join(final_results, by = "Exposure")


export(Table_S1,file = "result/Table/Table_S1.xlsx")



# ---------- 按季节绘制暴露相关性（四合一图） ----------
library(dplyr)
library(corrplot)
library(stringr)

# 设置字体（非Windows系统可去掉或使用 extrafont）
windowsFonts(Times = windowsFont("Times New Roman"))

# 打开图形设备，准备 2×2 布局
png(filename = "result/Figure/Figure S2.png",
    width = 10, height = 10, units = "in", res = 600,
    family = "serif")
par(mfrow = c(2, 2))

# 循环四个季节
for (s in 1:4) {
  
  # 筛选当前季节并提取暴露变量
  df <- data6 %>%
    filter(season.x == s) %>%
    select(all_of(exposures_name)) 
  
  # 计算 Spearman 相关及 p 值
  cor_matrix <- cor(df, method = "spearman", use = "pairwise.complete.obs")
  p_mat <- cor.mtest(df, conf.level = 0.95)$p
  
  # （可选）如需多重比较校正，请在此处处理，例如：
  # p_mat_adj <- matrix(p.adjust(p_mat, method = "holm"), ncol = ncol(p_mat))
  # 然后在下方的 p.mat 中使用 p_mat_adj
  
  # 第一步：上三角颜色 + 显著性标记
  corrplot(cor_matrix,
           method = "color",
           type = "upper",
           tl.col = "black",
           tl.srt = 45,
           tl.pos = "d",
           mar = c(0, 0, 2, 0),
           p.mat = p_mat,               # 用原始 p 值；若校正则改为 p_mat_adj
           sig.level = c(0.001, 0.01, 0.05),
           insig = "label_sig",
           pch.cex = 0.9,
           pch.col = "black")
  
  # 第二步：下三角纯黑色数字（叠加）
  corrplot(cor_matrix,
           method = "number",
           type = "lower",
           add = TRUE,
           tl.pos = "n",
           cl.pos = "n",
           diag = FALSE,
           number.cex = 0.8,
           col = "black")
  
  # 添加季节标题
  title(main = paste0("(",LETTERS[s],")"), line = 0.5, adj = 0)
}

dev.off()


