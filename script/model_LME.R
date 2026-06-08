library(tidyverse)
library(bruceR)

data7 <- import("data/processed/data7.rds")
source("script/define_variables.R")

# 建模前-查看分布 ####
summary(data7[,covariates])
summary(data7[,exposures_log_scaled])

# 建模-LME ####
library(lmerTest)
LME_table1 <- tibble()
for (out in c("liverhard_log","fat_decay_log")) {
  for (expo in exposures_log_scaled) {
    formular <- paste0(out, "~", expo, "+ Age + Gender + BMI + Smoke + Drink + Exercises + 
                       Education_new + diet_score + Hypertension_new + Hyperlipidemia_new + 
                       Diabetes_new + eGFR + (1|ID)")
    model <- lmer(formula = as.formula(formular), data = data7)
    df <- summary(model)
    df <- as.data.frame(df$coefficients)
    df <- df |> mutate(Exposure = expo, Out = out)
    LME_table1 <- bind_rows(LME_table1, df[2, ])
  }
}

LME_table2 <- LME_table1 |> 
  group_by(Out) %>% 
  mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_table1_sig <- LME_table2 |> filter(PFDR < 0.05)

# 百分比转化
source("func/functions.R")
LME_table3 <- LME_table2 %>% LME_percent_format()
export(LME_table3, file = "result/Table/Table S2.xlsx")
# 可视化 ####
plot1 <- LME_table3 |>
  # 1. 将所需列转化为数值型，特别是 PFDR
  mutate(across(c(Estimate, LOW_CI, UP_CI, PFDR), as.numeric)) |>
  
  # 2. 生成显著性星号标签
  mutate(P_lable = case_when(
    PFDR < 0.001 ~ "***",
    PFDR < 0.01  ~ "**",
    PFDR < 0.05  ~ "*",
    TRUE         ~ ""
  )) |>
  
  # 3. 重新编码 Out 变量，直接变成我们想要的分面标题
  mutate(Out = case_when(
    Out == "liverhard_log" ~ "(A) LSM",
    Out == "fat_decay_log" ~ "(B) CAP score"
  )) |>
  
  # 4. 固定因子顺序，确保出图的左右顺序和分面的上下顺序正确
  mutate(Exposure = factor(Exposure, levels = exposures_name)) |>
  mutate(Out = factor(Out, levels = c("(A) LSM", "(B) CAP score"))) |>
  
  # 开始绘图
  ggplot(aes(x = Exposure, y = Estimate)) +
  
  # 使用 facet_wrap 并允许 Y 轴自由缩放 (free_y)，按列垂直排列 (ncol = 1)
  facet_wrap(~ Out, scales = "free_y", ncol = 1) +
  
  # 添加 Y=0 的参考线
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) + 
  
  # 添加误差条
  geom_errorbar(aes(ymin = LOW_CI, ymax = UP_CI), width = 0.2, linewidth = 0.8) + 
  
  # 添加点
  geom_point(color = "black", stroke = 1, size = 2) + 
  
  # 添加 P 值星号（利用 vjust = -0.5 自动贴合在误差条上方，而不是固定加 0.25）
  geom_text(aes(x = Exposure, y = UP_CI, label = P_lable), 
            family = "Times New Roman", size = 6, fontface = "bold", vjust = -0.1) +
  
  # 调整坐标轴标签 (移除总标题，因为 A 和 B 已经作为分面标题了)
  labs(
    x = "",
    y = "Percent Change (95% CI) of LSM / CAP score"
  ) +
  
  # Y 轴数值格式化
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  
  # 主题美化
  theme_classic(base_size = 12) + 
  theme(
    text = element_text(family = "Times New Roman"), 
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    
    # 【核心配置】：将分面标签 (A) 和 (B) 加粗、靠左对齐 (hjust = 0)，并取消背景框
    strip.text = element_text(face = "bold", size = 12, hjust = 0),
    strip.background = element_blank() 
  )

ggsave(filename = "result/Figure/Figure 1.png", plot = plot1, device = "png",
       width = 7.5, height = 5, dpi = 300)
ggsave(filename = "result/Figure/Figure 1.pdf", plot = plot1, device = cairo_pdf,
       width = 7.5, height = 5, dpi = 300)

