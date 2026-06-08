library(tidyverse)
library(bruceR)

source("func/PDM.R")
source("func/functions.R")

data7 <- import("data/processed/data7.rds")

LSM_pathway <- import("result\\draft\\LSM_pos\\pathway_results.csv")
Zn_pathway <- import("result\\draft\\Zn_neg\\pathway_results.csv") 
code_book <- import("data/row/代谢物ID_名称.xlsx")

LSM_pathway_sig <- LSM_pathway %>% filter(`Holm adjust` < 0.05)
Zn_pathway_sig <- Zn_pathway %>% filter(`Holm adjust` < 0.05)

pathway_merge <- inner_join(LSM_pathway_sig, Zn_pathway_sig, by = "V1")

## 读取KEGG 数据库 ####
library(KEGGREST)
pathway_list <- keggList("pathway", "hsa")

MWAS_Zn_V2 <- import("result/draft/MWAS_Zn_V2.xlsx")
MWAS_LSM_V2 <- import("result/draft/MWAS_LSM_V2.xlsx")

MWAS_Zn_V2_sig <- MWAS_Zn_V2 %>% 
  filter(PFDR < 0.2) %>% 
  filter(Estimate < 0) %>% 
  pull(KEGG)
MWAS_LSM_V2_sig <- MWAS_LSM_V2 %>% 
  filter(PFDR < 0.2) %>% 
  filter(Estimate > 0) %>% 
  pull(KEGG)

merge_metabolite <- unique(c(MWAS_Zn_V2_sig,MWAS_LSM_V2_sig))

path_meta_table <- tibble()
for (path in c(pathway_merge$V1)) {
  pathway <- grep(path, pathway_list, value = TRUE,fixed = TRUE)
  metabolites <- keggGet(names(pathway))
  metabolites <- metabolites[[1]]
  metabolites <- metabolites$COMPOUND
  metabolites <- intersect(names(metabolites), merge_metabolite)
  df_table <- tibble(Pathway_name = path, KEGG = metabolites)
  path_meta_table <- bind_rows(path_meta_table, df_table)
}

path_meta_table <- path_meta_table %>% left_join(code_book, by = "KEGG")
export(path_meta_table, file = "result/draft/path_meta_table.xlsx")


# 检查是否已加载 lmerTest 包
if ("package:lmerTest" %in% search()) {
  # 如果已加载，则卸载
  tryCatch(
    {
      detach("package:lmerTest", unload = TRUE, character.only = TRUE)
      cat("成功卸载 lmerTest 包\n")
    },
    error = function(e) {
      cat("卸载 lmerTest 包时出错:", e$message, "\n")
    }
  )
} else {
  # 如果未加载，则提示
  cat("lmerTest 包未加载，无需卸载\n")
}

## 通路中介 ####
library(mediation)
set.seed(2024)
mediation_table_pathway <- tibble()
for (path in unique(c(path_meta_table$Pathway_name))) {
  KEGGs <- path_meta_table |>
    filter(Pathway_name == path) |>
    mutate(Metabolite_ID = paste0(Metabolite_ID,"_log_scaled")) %>% 
    pull(Metabolite_ID)
  
  A <- data7[, "Zn_log_scaled"] |> as.matrix()
  Y <- data7[, "liverhard_log"] |> as.matrix()
  m1 <- data7[, KEGGs] |> as.matrix()
  m1 <- na.omit(m1)
  
  output <- PDM_1(
    x = A, y = Y, m = m1,
    imax = 1000, # 最大迭代次数
    tol = 1e-5, # 收敛阈值
    theta = rep(1, 5), # 参数初始值
    w1 = rep(1, ncol(m1)), # 初始权重（等权重）
    interval = 1e6, # 优化参数
    step = 1e4 # 优化步长
  )
  
  est_w <- output$w1
  # 将dm1添加为数据框的新列（确保与数据框行数一致）
  dm1 <- m1 %*% as.matrix(est_w, ncol = 1)
  
  data7$dm1 <- as.vector(dm1) # 转换矩阵为向量并添加到数据框
  
  # 修改后的模型公式（直接使用dm1变量名）
  model1 <- lmer(
    data = data7,
    formula = dm1 ~ Zn_log_scaled + Age + Gender + BMI +
      Smoke + Drink + Exercises + Education +
      diet_score + Hypertension_new +
      Hyperlipidemia_new + Diabetes_new + eGFR + (1 | ID)
  )
  
  model2 <- lmer(
    data = data7,
    formula = liverhard_log ~ dm1 + Zn_log_scaled + Age + Gender +
      BMI + Smoke + Drink + Exercises + Education +
      diet_score + Hypertension_new +
      Hyperlipidemia_new + Diabetes_new + eGFR + (1 | ID)
  )
  
  # 调整mediate参数（使用dm1作为中介变量）
  medi_model <- mediate(model1, model2,
                        treat = "Zn_log_scaled", # 直接暴露变量名
                        mediator = "dm1", # 固定为dm1
                        sims = 2000
  )
  
  df <- summary(medi_model)
  df_table <- tibble(
    name = c("ACME", "ADE", "Total Effect", "Prop. Mediated"),
    Estimate = c(df$d0, df$z0, df$tau.coef, df$n0),
    `95% CI Lower` = c(df$d0.ci[1], df$z0.ci[1], df$tau.ci[1], df$n0.ci[1]),
    `95% CI Upper` = c(df$d0.ci[2], df$z0.ci[2], df$tau.ci[2], df$n0.ci[2]),
    `p-value` = c(df$d0.p, df$z0.p, df$tau.p, df$n0.p)
  ) %>% mutate(
    exposure_name = "Zn",
    medi_name = path,
    outcome_name = "LSM"
  )
  mediation_table_pathway <- rbind(mediation_table_pathway, df_table)
}
export(mediation_table_pathway,"result/draft/mediation_table_pathway_row.xlsx")

mediation_table_pathway_1 <- mediation_table_pathway %>%
  filter(name != "Prop. Mediated") %>%
  mutate(
    Estimate = (10^Estimate - 1) * 100,
    `95% CI Lower` = (10^`95% CI Lower` - 1) * 100,
    `95% CI Upper` = (10^`95% CI Upper` - 1) * 100
  )

mediation_table_pathway_2 <- mediation_table_pathway %>%
  filter(name == "Prop. Mediated")

mediation_table_pathway_V2 <- bind_rows(mediation_table_pathway_1, mediation_table_pathway_2)

mediation_table_pathway_V2 <- mediation_table_pathway_V2 %>% process_mediation_results()
export(mediation_table_pathway_V2, file = "result/draft/pathway_mediation.xlsx")
mediation_table_pathway_V2 <- import("result/draft/pathway_mediation.xlsx")
mediation_table_pathway_V3 <- mediation_table_pathway_V2 %>% 
  mutate(PFDR_pathway = p.adjust(`p-value_Prop. Mediated`,"BH"))
mediation_table_pathway_V3 %>% filter(PFDR_pathway < 0.05)



## 代谢物中介 ####
library(mediation)
set.seed(2024)
mediation_table_metabolite <- tibble()
for (medi in c(unique(path_meta_table$Metabolite_ID))) {
  Metabolite_ID_new <- paste0(medi,"_log_scaled")
  
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
  mediation_table_metabolite <- rbind(mediation_table_metabolite, df_table)
}

export(mediation_table_metabolite,file = "result/draft/mediation_table_metabolite_row.xlsx")

mediation_table_metabolite_1 <- mediation_table_metabolite %>%
  filter(name != "Prop. Mediated") %>%
  mutate(
    Estimate = (10^Estimate - 1) * 100,
    `95% CI Lower` = (10^`95% CI Lower` - 1) * 100,
    `95% CI Upper` = (10^`95% CI Upper` - 1) * 100
  )

mediation_table_metabolite_2 <- mediation_table_metabolite %>%
  filter(name == "Prop. Mediated")

mediation_table_metabolite_V2 <- bind_rows(mediation_table_metabolite_1, mediation_table_metabolite_2)

mediation_table_metabolite_V2 <- mediation_table_metabolite_V2 %>% process_mediation_results()

export(mediation_table_metabolite_V2,file = "result/draft/mediation_table_metabolite_V2.xlsx")

## 合并通路和 代谢物 ####
mediation_table_metabolite_V2 <- import("result/draft/mediation_table_metabolite_V2.xlsx")

Table_2 <- mediation_table_pathway_V3 %>% 
  rename(Pathway_name = medi_name) %>% 
  left_join(path_meta_table, by = "Pathway_name") %>% 
  mutate(medi_name = Metabolite_ID ) %>% 
  left_join(mediation_table_metabolite_V2, by = "medi_name") %>% 
  dplyr::select(all_of(c("exposure_name.x", "Pathway_name", "ACME (β, 95%CI).x", "p-value_ACME.x", "Prop. Mediated.x",
                  "PFDR_pathway","p-value_Prop. Mediated.x","Metabolite_Name","ACME (β, 95%CI).y", "p-value_ACME.y", 
                  "Prop. Mediated.y","p-value_Prop. Mediated.y")))

Table_2 <- Table_2 %>% 
  group_by(Pathway_name) %>% 
  mutate(PFDR_metabolite = p.adjust(`p-value_Prop. Mediated.y`, "BH"))
  
Table_2 %>% filter(PFDR_metabolite < 0.05 & PFDR_pathway <0.05)

  
export(Table_2, file = "result/Table/Table 2.xlsx")


# 绘制桑基图 ####
plot_data <- import("result/draft/桑基图数据.xlsx")

library(ggsankey)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggsci)

# 1. 长标签自动换行
plot_data <- plot_data %>%
  mutate(
    pathway = str_wrap(pathway, width = 25),
    metabolite = str_wrap(metabolite, width = 20)
  )

# 2. 转换成长格式
df_long <- plot_data %>%
  make_long(Exposure, pathway, metabolite, LSM)

# 3. 绘图（直接用 geom_sankey 控制节点与流线）
P <- ggplot(df_long, aes(
  x = x,
  next_x = next_x,
  node = node,
  next_node = next_node,
  fill = factor(node),   # 流线和节点都按节点名称填充
  label = node           # 节点文字
)) +
  geom_sankey(
    flow.alpha = 0.2,    # 流线透明度保持 50%
    alpha = 0.2,         # 新增：节点透明度设为 70%
    node.color = "white",# 节点白色边框
    size = 0.3           # 流线边框粗细
  ) +            
  geom_sankey_text(
    size = 3,
    color = "black",
    family = "serif",
    fontface = "plain"
  ) +
  ggsci::scale_fill_nejm() +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(family = "serif", size = 14, 
                              hjust = 0.5, margin = margin(b = 8))
  )

ggsave("result/Figure/Figure sankey.pdf", plot = P, device = cairo_pdf,height = 3.5, width = 6
       ,dpi = 300)
