library(tidyverse)
library(bruceR)

data7 <- import("data/processed/data7.rds")
source("script/define_variables.R")
source("func/functions.R")
# 去掉高危险性的 LSM ####
data8 <- data7 |> filter(HepatitisB == 0)
nrow(data8)
length(unique(data7$ID))
data8 <- data7 |> filter(Chronichepatitis == 0)
nrow(data8)

library(lmerTest)
LME_table5 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates, collapse = "+"), "+ (1|ID)")
  model <- lmer(formula = as.formula(formular), data = data8)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_table5 <- bind_rows(LME_table5, df[2, ])
}

LME_table5 <- LME_table5 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_table5_V2 <- LME_table5 |> LME_percent_format()



plot5 <- LME_table5_V2 |>
  mutate(P_lable = ifelse(PFDR < 0.05, "*", ifelse(PFDR < 0.01, "**", ifelse(PFDR < 0.001, "***", "")))) |>
  mutate(across(
    c(Estimate, LOW_CI, UP_CI),
    ~ as.numeric(.x)
  )) |>
  mutate(Exposure = factor(Exposure, levels = exposures_name)) |>
  ggplot(aes(x = Exposure, y = Estimate)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.8
  ) + # 优化参考线
  geom_errorbar(aes(ymin = LOW_CI, ymax = UP_CI),
                width = 0.2,
                linewidth = 0.8
  ) + # 误差条样式
  geom_point(
    stroke = 1
  ) + # 点样式（填充+边框）
  geom_text(aes(x = Exposure, y = UP_CI + 0.25, label = P_lable), family = "Times New Roman", size = 7, fontface = "bold") +
  labs(
    x = "",
    y = "%Change (95% CI) of LSM",
    title = NULL # 科研图通常不加标题
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + # 新增此行
  theme_classic(base_size = 12) + # 经典无网格主题
  theme(
    text = element_text(family = "Times New Roman"), # 指定字体
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
  )


ggsave(
  filename = "result/Figure/Figure S5.png", plot = plot5, device = "png",
  width = 7.5, height = 3.5, dpi = 300
)


# 去掉egfr ####
LME_table5 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-12], collapse = "+"), "+ (1|ID)")
  model <- lmer(formula = as.formula(formular), data = data7)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_table5 <- bind_rows(LME_table5, df[2, ])
}

LME_table5 <- LME_table5 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_table5_V2 <- LME_table5 |> LME_percent_format()



plot5 <- LME_table5_V2 |>
  mutate(P_lable = ifelse(PFDR < 0.05, "*", ifelse(PFDR < 0.01, "**", ifelse(PFDR < 0.001, "***", "")))) |>
  mutate(across(
    c(Estimate, LOW_CI, UP_CI),
    ~ as.numeric(.x)
  )) |>
  mutate(Exposure = factor(Exposure, levels = exposures_name)) |>
  ggplot(aes(x = Exposure, y = Estimate)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.8
  ) + # 优化参考线
  geom_errorbar(aes(ymin = LOW_CI, ymax = UP_CI),
                width = 0.2,
                linewidth = 0.8
  ) + # 误差条样式
  geom_point(
    stroke = 1
  ) + # 点样式（填充+边框）
  geom_text(aes(x = Exposure, y = UP_CI + 0.25, label = P_lable), family = "Times New Roman", size = 7, fontface = "bold") +
  labs(
    x = "",
    y = "%Change (95% CI) of LSM",
    title = NULL # 科研图通常不加标题
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + # 新增此行
  theme_classic(base_size = 12) + # 经典无网格主题
  theme(
    text = element_text(family = "Times New Roman"), # 指定字体
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
  )

ggsave(
  filename = "result/Figure/Figure S6.png", plot = plot5, device = "png",
  width = 7.5, height = 3.5, dpi = 300
)

# 分层 ####
### 年龄 ####
summary(data7[,c(covariates)])
data7_0 <- data7 |> filter(age < 60)
LME_0 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-1], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_0)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_0 <- bind_rows(LME_0, df[2, ])
}
LME_0 <- LME_0 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_0_V2 <- LME_0 |> LME_percent_format()
LME_0_V2 |> filter(PFDR < 0.05)



data7_1 <- data7 |> filter(age >= 60)
LME_1 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-1], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_1)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_1 <- bind_rows(LME_1, df[2, ])
}
LME_1 <- LME_1 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_1_V2 <- LME_1 |> LME_percent_format()
LME_1_V2 |> filter(PFDR < 0.05)


data7_inter <- data7 |> mutate(age_class = ifelse(age >= 60, 1, 0), age_class = as.factor(age_class))
LME_inter <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "*", "age_class", "+", paste0(covariates[-1], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_inter)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_inter <- bind_rows(LME_inter, tail(df, 1))
}
LME_inter <- LME_inter |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_inter_V2 <- LME_inter |> LME_percent_format()
LME_inter_V2 |> filter(`Pr(>|t|)` < 0.1)


LME_0_V2 <- LME_0_V2 |> mutate(group = "< 60 years")
LME_1_V2 <- LME_1_V2 |> mutate(group = "≥ 60 years")

data_plot <- bind_rows(LME_0_V2, LME_1_V2)
plot5 <- list()
plot5[["age"]] <- data_plot |>
  mutate(
    P_lable = case_when(
      PFDR < 0.001 ~ "***",
      PFDR < 0.01 ~ "**",
      PFDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    across(c(Estimate, LOW_CI, UP_CI), as.numeric),
    # 创建组合标签用于x轴
    Exposure = factor(Exposure, levels = unique(data_plot$Exposure))
  ) |>
  ggplot(aes(x = Exposure, y = Estimate, shape = group)) + # 颜色分组
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    shape = "grey40",
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(ymin = LOW_CI, ymax = UP_CI),
    width = 0.2,
    position = position_dodge(width = 0.5), # 防止重叠
    linewidth = 0.8
  ) +
  geom_point(
    stroke = 1,
    position = position_dodge(width = 0.5) # 与误差条对齐
  ) +
  geom_text(
    aes(y = UP_CI + 0.25, label = P_lable),
    position = position_dodge(width = 0.5),
    family = "Times New Roman",
    size = 5, # 减小星号尺寸
    fontface = "bold"
  ) +
  labs(x = "", y = "%Change (95% CI) ", shape = "Age") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Times New Roman"), # 指定字体
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(hjust = 0, face = "bold"), # 标题左对齐
    legend.text = element_text(hjust = 0),                # 标签左对齐
    legend.justification = "left",                        # 图例在区域内左对齐
    legend.box.just = "left",                             # 图例框左对齐
    legend.margin = margin(l = 0, r = 0),                 # 彻底消除图例左右两侧的间距
    legend.box.margin = margin(l = 0, r = 0),              # 彻底消除外部容器间距
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )

### BMI ####
data7_0 <- data7 |> filter(BMI < 24)
LME_0 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-3], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_0)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_0 <- bind_rows(LME_0, df[2, ])
}
LME_0 <- LME_0 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_0_V2 <- LME_0 |> LME_percent_format()
LME_0_V2 |> filter(PFDR < 0.05)



data7_1 <- data7 |> filter(BMI >= 24)
LME_1 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-3], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_1)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_1 <- bind_rows(LME_1, df[2, ])
}
LME_1 <- LME_1 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_1_V2 <- LME_1 |> LME_percent_format()
LME_1_V2 |> filter(PFDR < 0.05)


data7_inter <- data7 |> mutate(BMI_class = ifelse(BMI < 24, 1, 0), BMI_class = as.factor(BMI_class))
LME_inter <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "*", "BMI_class", "+", paste0(covariates[-3], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_inter)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_inter <- bind_rows(LME_inter, tail(df, 1))
}
LME_inter <- LME_inter |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_inter_V2 <- LME_inter |> LME_percent_format()
LME_inter_V2 |> filter(`Pr(>|t|)` < 0.1)


LME_0_V2 <- LME_0_V2 |> mutate(group = "< 24")
LME_1_V2 <- LME_1_V2 |> mutate(group = "≥ 24")

data_plot <- bind_rows(LME_0_V2, LME_1_V2)

plot5[["BMI"]] <- data_plot |>
  mutate(
    P_lable = case_when(
      PFDR < 0.001 ~ "***",
      PFDR < 0.01 ~ "**",
      PFDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    across(c(Estimate, LOW_CI, UP_CI), as.numeric),
    # 创建组合标签用于x轴
    Exposure = factor(Exposure, levels = unique(data_plot$Exposure))
  ) |>
  ggplot(aes(x = Exposure, y = Estimate, shape = group)) + # 颜色分组
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    shape = "grey40",
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(ymin = LOW_CI, ymax = UP_CI, shape = group),
    width = 0.2,
    position = position_dodge(width = 0.5), # 防止重叠
    linewidth = 0.8
  ) +
  geom_point(
    stroke = 1,
    position = position_dodge(width = 0.5) # 与误差条对齐
  ) +
  geom_text(
    aes(y = UP_CI + 0.25, label = P_lable),
    position = position_dodge(width = 0.5),
    family = "Times New Roman",
    size = 5, # 减小星号尺寸
    fontface = "bold"
  ) +
  labs(x = "", y = "%Change (95% CI) ", shape = expression(BMI ~ (kg / m^2))) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Times New Roman"), # 指定字体
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(hjust = 0, face = "bold"), # 标题左对齐
    legend.text = element_text(hjust = 0),                # 标签左对齐
    legend.justification = "left",                        # 图例在区域内左对齐
    legend.box.just = "left",                             # 图例框左对齐
    legend.margin = margin(l = 0, r = 0),                 # 彻底消除图例左右两侧的间距
    legend.box.margin = margin(l = 0, r = 0),              # 彻底消除外部容器间距
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )

### Hyperlipidemia_new ####
data7_0 <- data7 |> filter(Hyperlipidemia_new == 0)
LME_0 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-10], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_0)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_0 <- bind_rows(LME_0, df[2, ])
}
LME_0 <- LME_0 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_0_V2 <- LME_0 |> LME_percent_format()
LME_0_V2 |> filter(PFDR < 0.05)



data7_1 <- data7 |> filter(Hyperlipidemia_new == 1)
LME_1 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-10], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_1)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_1 <- bind_rows(LME_1, df[2, ])
}
LME_1 <- LME_1 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_1_V2 <- LME_1 |> LME_percent_format()
LME_1_V2 |> filter(PFDR < 0.05)


LME_inter <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "*", "Hyperlipidemia_new", "+", paste0(covariates[-10], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_inter)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_inter <- bind_rows(LME_inter, tail(df, 1))
}
LME_inter <- LME_inter |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_inter_V2 <- LME_inter |> LME_percent_format()
LME_inter_V2 |> filter(`Pr(>|t|)` < 0.1)



LME_0_V2 <- LME_0_V2 |> mutate(group = "No")
LME_1_V2 <- LME_1_V2 |> mutate(group = "Yes")

data_plot <- bind_rows(LME_0_V2, LME_1_V2)

plot5[["Hyperlipidemia"]] <- data_plot |>
  mutate(
    P_lable = case_when(
      PFDR < 0.001 ~ "***",
      PFDR < 0.01 ~ "**",
      PFDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    across(c(Estimate, LOW_CI, UP_CI), as.numeric),
    # 创建组合标签用于x轴
    Exposure = factor(Exposure, levels = unique(data_plot$Exposure))
  ) |>
  ggplot(aes(x = Exposure, y = Estimate, shape = group)) + # 颜色分组
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    shape = "grey40",
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(ymin = LOW_CI, ymax = UP_CI, shape = group),
    width = 0.2,
    position = position_dodge(width = 0.5), # 防止重叠
    linewidth = 0.8
  ) +
  geom_point(
    stroke = 1,
    position = position_dodge(width = 0.5) # 与误差条对齐
  ) +
  geom_text(
    aes(y = UP_CI + 0.25, label = P_lable),
    position = position_dodge(width = 0.5),
    family = "Times New Roman",
    size = 5, # 减小星号尺寸
    fontface = "bold"
  ) +
  labs(x = "", y = "%Change (95% CI) ", shape = "Hyperlipidemia") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Times New Roman"), # 指定字体
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(hjust = 0, face = "bold"), # 标题左对齐
    legend.text = element_text(hjust = 0),                # 标签左对齐
    legend.justification = "left",                        # 图例在区域内左对齐
    legend.box.just = "left",                             # 图例框左对齐
    legend.margin = margin(l = 0, r = 0),                 # 彻底消除图例左右两侧的间距
    legend.box.margin = margin(l = 0, r = 0),              # 彻底消除外部容器间距
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )


### Gender ####
data7_0 <- data7 |> filter(Gender == 1)
LME_0 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-2], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_0)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_0 <- bind_rows(LME_0, df[2, ])
}
LME_0 <- LME_0 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_0_V2 <- LME_0 |> LME_percent_format()
LME_0_V2 |> filter(PFDR < 0.05)



data7_1 <- data7 |> filter(Gender == 2)
LME_1 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-2], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_1)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_1 <- bind_rows(LME_1, df[2, ])
}
LME_1 <- LME_1 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_1_V2 <- LME_1 |> LME_percent_format()
LME_1_V2 |> filter(PFDR < 0.05)



LME_inter <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "*", "Gender", "+", paste0(covariates[-2], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_inter)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_inter <- bind_rows(LME_inter, tail(df, 1))
}
LME_inter <- LME_inter |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_inter_V2 <- LME_inter |> LME_percent_format()
LME_inter_V2 |> filter(`Pr(>|t|)` < 0.1)




LME_0_V2 <- LME_0_V2 |> mutate(group = "Male")
LME_1_V2 <- LME_1_V2 |> mutate(group = "Female")

data_plot <- bind_rows(LME_0_V2, LME_1_V2)

plot5[["Gender"]] <- data_plot |>
  mutate(
    P_lable = case_when(
      PFDR < 0.001 ~ "***",
      PFDR < 0.01 ~ "**",
      PFDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    across(c(Estimate, LOW_CI, UP_CI), as.numeric),
    # 创建组合标签用于x轴
    Exposure = factor(Exposure, levels = unique(data_plot$Exposure))
  ) |>
  ggplot(aes(x = Exposure, y = Estimate, shape = group)) + # 颜色分组
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    shape = "grey40",
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(ymin = LOW_CI, ymax = UP_CI, shape = group),
    width = 0.2,
    position = position_dodge(width = 0.5), # 防止重叠
    linewidth = 0.8
  ) +
  geom_point(
    stroke = 1,
    position = position_dodge(width = 0.5) # 与误差条对齐
  ) +
  geom_text(
    aes(y = UP_CI + 0.25, label = P_lable),
    position = position_dodge(width = 0.5),
    family = "Times New Roman",
    size = 5, # 减小星号尺寸
    fontface = "bold"
  ) +
  labs(x = "", y = "%Change (95% CI) ", shape = "Sex") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Times New Roman"), # 指定字体
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(hjust = 0, face = "bold"), # 标题左对齐
    legend.text = element_text(hjust = 0),                # 标签左对齐
    legend.justification = "left",                        # 图例在区域内左对齐
    legend.box.just = "left",                             # 图例框左对齐
    legend.margin = margin(l = 0, r = 0),                 # 彻底消除图例左右两侧的间距
    legend.box.margin = margin(l = 0, r = 0),              # 彻底消除外部容器间距
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )

### Smoke ####

data7_0 <- data7 |> filter(Smoke == 1)
LME_0 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-4], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_0)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_0 <- bind_rows(LME_0, df[2, ])
}
LME_0 <- LME_0 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_0_V2 <- LME_0 |> LME_percent_format()
LME_0_V2 |> filter(PFDR < 0.05)



data7_1 <- data7 |> filter(Smoke != 1)
LME_1 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-4], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_1)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_1 <- bind_rows(LME_1, df[2, ])
}
LME_1 <- LME_1 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_1_V2 <- LME_1 |> LME_percent_format()
LME_1_V2 |> filter(PFDR < 0.05)


data7_inter <- data7 %>% mutate(Smoke_class = ifelse(Smoke == 1, 1, 0), Smoke_class = as.factor(Smoke_class))
LME_inter <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "*", "Smoke_class", "+", paste0(covariates[-4], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_inter)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_inter <- bind_rows(LME_inter, tail(df, 1))
}
LME_inter <- LME_inter |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_inter_V2 <- LME_inter |> LME_percent_format()
LME_inter_V2 |> filter(`Pr(>|t|)` < 0.05)




LME_0_V2 <- LME_0_V2 |> mutate(group = "Current")
LME_1_V2 <- LME_1_V2 |> mutate(group = "Never/Former")

data_plot <- bind_rows(LME_0_V2, LME_1_V2)

plot5[["Smoke"]] <- data_plot |>
  mutate(
    P_lable = case_when(
      PFDR < 0.001 ~ "***",
      PFDR < 0.01 ~ "**",
      PFDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    across(c(Estimate, LOW_CI, UP_CI), as.numeric),
    # 创建组合标签用于x轴
    Exposure = factor(Exposure, levels = unique(data_plot$Exposure))
  ) |>
  ggplot(aes(x = Exposure, y = Estimate, shape = group)) + # 颜色分组
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    shape = "grey40",
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(ymin = LOW_CI, ymax = UP_CI, shape = group),
    width = 0.2,
    position = position_dodge(width = 0.5), # 防止重叠
    linewidth = 0.8
  ) +
  geom_point(
    stroke = 1,
    position = position_dodge(width = 0.5) # 与误差条对齐
  ) +
  geom_text(
    aes(y = UP_CI + 0.25, label = P_lable),
    position = position_dodge(width = 0.5),
    family = "Times New Roman",
    size = 5, # 减小星号尺寸
    fontface = "bold"
  ) +
  labs(x = "", y = "%Change (95% CI) ", shape = "Smoking status") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Times New Roman"), # 指定字体
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(hjust = 0, face = "bold"), # 标题左对齐
    legend.text = element_text(hjust = 0),                # 标签左对齐
    legend.justification = "left",                        # 图例在区域内左对齐
    legend.box.just = "left",                             # 图例框左对齐
    legend.margin = margin(l = 0, r = 0),                 # 彻底消除图例左右两侧的间距
    legend.box.margin = margin(l = 0, r = 0),              # 彻底消除外部容器间距
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )


### Drink ####

data7_0 <- data7 |> filter(Drink == 1)
LME_0 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-5], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_0)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_0 <- bind_rows(LME_0, df[2, ])
}
LME_0 <- LME_0 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_0_V2 <- LME_0 |> LME_percent_format()
LME_0_V2 |> filter(PFDR < 0.05)



data7_1 <- data7 |> filter(Drink != 1)
LME_1 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-5], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_1)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_1 <- bind_rows(LME_1, df[2, ])
}
LME_1 <- LME_1 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_1_V2 <- LME_1 |> LME_percent_format()
LME_1_V2 |> filter(PFDR < 0.05)


data7_inter <- data7 %>% mutate(Drink_class = ifelse(Drink == 1, 1, 0), Drink_class = as.factor(Drink_class))
LME_inter <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "*", "Drink_class", "+", paste0(covariates[-4], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_inter)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_inter <- bind_rows(LME_inter, tail(df, 1))
}
LME_inter <- LME_inter |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_inter_V2 <- LME_inter |> LME_percent_format()
LME_inter_V2 |> filter(`Pr(>|t|)` < 0.1)




LME_0_V2 <- LME_0_V2 |> mutate(group = "Current")
LME_1_V2 <- LME_1_V2 |> mutate(group = "Never/Former")

data_plot <- bind_rows(LME_0_V2, LME_1_V2)

plot5[["Drink"]] <- data_plot |>
  mutate(
    P_lable = case_when(
      PFDR < 0.001 ~ "***",
      PFDR < 0.01 ~ "**",
      PFDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    across(c(Estimate, LOW_CI, UP_CI), as.numeric),
    # 创建组合标签用于x轴
    Exposure = factor(Exposure, levels = unique(data_plot$Exposure))
  ) |>
  ggplot(aes(x = Exposure, y = Estimate, shape = group)) + # 颜色分组
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    shape = "grey40",
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(ymin = LOW_CI, ymax = UP_CI, shape = group),
    width = 0.2,
    position = position_dodge(width = 0.5), # 防止重叠
    linewidth = 0.8
  ) +
  geom_point(
    stroke = 1,
    position = position_dodge(width = 0.5) # 与误差条对齐
  ) +
  geom_text(
    aes(y = UP_CI + 0.25, label = P_lable),
    position = position_dodge(width = 0.5),
    family = "Times New Roman",
    size = 5, # 减小星号尺寸
    fontface = "bold"
  ) +
  labs(x = "", y = "%Change (95% CI) ", shape = "Drinking status") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Times New Roman"), # 指定字体
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(hjust = 0, face = "bold"), # 标题左对齐
    legend.text = element_text(hjust = 0),                # 标签左对齐
    legend.justification = "left",                        # 图例在区域内左对齐
    legend.box.just = "left",                             # 图例框左对齐
    legend.margin = margin(l = 0, r = 0),                 # 彻底消除图例左右两侧的间距
    legend.box.margin = margin(l = 0, r = 0),              # 彻底消除外部容器间距
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )


### Exercises ####
data7_0 <- data7 |> filter(Exercises == 0)
LME_0 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-6], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_0)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_0 <- bind_rows(LME_0, df[2, ])
}
LME_0 <- LME_0 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_0_V2 <- LME_0 |> LME_percent_format()
LME_0_V2 |> filter(PFDR < 0.05)



data7_1 <- data7 |> filter(Exercises == 1)
LME_1 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates[-6], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_1)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_1 <- bind_rows(LME_1, df[2, ])
}
LME_1 <- LME_1 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_1_V2 <- LME_1 |> LME_percent_format()
LME_1_V2 |> filter(PFDR < 0.05)



LME_inter <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "*", "Exercises", "+", paste0(covariates[-6], collapse = "+"), "+ (1|ID) ")
  model <- lmer(formula = as.formula(formular), data = data7_inter)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_inter <- bind_rows(LME_inter, tail(df, 1))
}
LME_inter <- LME_inter |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_inter_V2 <- LME_inter |> LME_percent_format()
LME_inter_V2 |> filter(`Pr(>|t|)` < 0.1)




LME_0_V2 <- LME_0_V2 |> mutate(group = "No")
LME_1_V2 <- LME_1_V2 |> mutate(group = "Yes")

data_plot <- bind_rows(LME_0_V2, LME_1_V2)

plot5[["Exercises"]] <- data_plot |>
  mutate(
    P_lable = case_when(
      PFDR < 0.001 ~ "***",
      PFDR < 0.01 ~ "**",
      PFDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    across(c(Estimate, LOW_CI, UP_CI), as.numeric),
    # 创建组合标签用于x轴
    Exposure = factor(Exposure, levels = unique(data_plot$Exposure))
  ) |>
  ggplot(aes(x = Exposure, y = Estimate, shape = group)) + # 颜色分组
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    shape = "grey40",
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(ymin = LOW_CI, ymax = UP_CI, shape = group),
    width = 0.2,
    position = position_dodge(width = 0.5), # 防止重叠
    linewidth = 0.8
  ) +
  geom_point(
    stroke = 1,
    position = position_dodge(width = 0.5) # 与误差条对齐
  ) +
  geom_text(
    aes(y = UP_CI + 0.25, label = P_lable),
    position = position_dodge(width = 0.5),
    family = "Times New Roman",
    size = 5, # 减小星号尺寸
    fontface = "bold"
  ) +
  labs(x = "", y = "%Change (95% CI) ", shape = "Exercises") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Times New Roman"), # 指定字体
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    
    legend.title = element_text(hjust = 0, face = "bold"), # 标题左对齐
    legend.text = element_text(hjust = 0),                # 标签左对齐
    legend.justification = "left",                        # 图例在区域内左对齐
    legend.box.just = "left",                             # 图例框左对齐
    legend.margin = margin(l = 0, r = 0),                 # 彻底消除图例左右两侧的间距
    legend.box.margin = margin(l = 0, r = 0)              # 彻底消除外部容器间距
  )


library(gridExtra)
library(cowplot)
plots_legend <- gridExtra::grid.arrange(plot5$age, plot5$Gender, plot5$BMI, plot5$Smoke, plot5$Drink, plot5$Exercises, plot5$Hyperlipidemia, nrow = 7)

legend_names <- c("age", "Gender", "BMI", "Smoke", "Drink", "Exercises", "Hyperlipidemia")

library(showtext)
showtext_auto()

# 开启设备
png(tempfile()) 

legend_list <- lapply(legend_names, function(x) {
  get_legend(
    plot5[[x]] + 
      theme(
        # 使用通用别名 serif，它在大多数系统中自动对应 Times New Roman
        # 且不会触发 PDF CID 字体映射错误
        text = element_text(family = "serif"), 
        legend.text = element_text(family = "serif"),
        legend.title = element_text(family = "serif"),
        legend.justification = "left", 
        legend.box.just = "left",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,0,0)
      )
  )
})

dev.off()
# 此时运行 legend_list 应该就不会报错了


library(showtext)
showtext_auto(FALSE)
library(cowplot)
library(ggplot2)

# Open a standard windows/X11 device instead of PDF/AGG for the calculation
# If on Windows:
if(.Platform$OS.type == "windows") windows() else x11()

legend_list <- lapply(legend_names, function(x) {
  # Add a safe theme to the plot BEFORE extracting the legend
  p_safe <- plot5[[x]] + 
    theme(
      text = element_text(family = "sans"), # Force a standard system font
      legend.text = element_text(family = "sans"),
      legend.title = element_text(family = "sans"),
      legend.justification = "left", 
      legend.box.just = "left",
      legend.margin = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0)
    )
  
  # Extract the legend
  return(get_legend(p_safe))
})

# Close the temporary window
dev.off()

legend_plots <- lapply(legend_list, function(lgd) {
  ggdraw() + draw_grob(lgd, x = 0, hjust = 0) # 强制在坐标轴 x=0 的位置开始画
})

combined_legends <- plot_grid(
  plotlist = legend_plots,
  ncol = 1,
  align = "v",
  axis = "l", 
  rel_heights = rep(1, length(legend_plots)) # 确保每一行高度一致
)

my_labels <- paste0("(", LETTERS[1:length(legend_names)], ")")

main_list <- lapply(1:length(legend_names), function(i) {
  name_label <- legend_names[i]
  letter_label <- my_labels[i]
  
  # 注意：这里每一项之间都必须用 + 连接
  plot5[[name_label]] + 
    labs(title = letter_label) + 
    theme(
      legend.position = "none",
      # 这里的 hjust = 0 确保标题在最左侧
      plot.title = element_text(
        family = "Times New Roman", 
        face = "bold", 
        size = 14, 
        hjust = 0, 
        margin = margin(b = 5) # 给标题和图之间留点空隙
      )
    )
})

# 合并主图列
combined_mian <- plot_grid(
  plotlist = main_list,
  ncol = 1,
  align = "v" # 确保多张图的坐标轴对齐
)

# 最终合并
final_figure <- plot_grid(
  combined_mian, 
  combined_legends, 
  ncol = 2, 
  rel_widths = c(0.9, 0.1),
  align = "h", # 关键：确保左边的图和右边的图例水平对齐
  axis = "t"   # 顶部对齐
)


ggsave(
  filename = "result/Figure/Figure S7.png", plot = final_figure, device = "png",
  width = 12.5, height = 15, dpi = 1200
)

# 反向因果-通路 ####
code_book <- import("data/row/代谢物ID_名称.xlsx")
path_meta_table <- import("result/draft/path_meta_table.xlsx")
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
mediation_sens1 <- tibble()
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
    formula = liverhard_log ~ Zn_log_scaled + Age + Gender + BMI +
      Smoke + Drink + Exercises + Education +
      diet_score + Hypertension_new +
      Hyperlipidemia_new + Diabetes_new + eGFR + (1 | ID)
  )
  
  model2 <- lmer(
    data = data7,
    formula = dm1 ~ liverhard_log + Zn_log_scaled + Age + Gender +
      BMI + Smoke + Drink + Exercises + Education +
      diet_score + Hypertension_new +
      Hyperlipidemia_new + Diabetes_new + eGFR + (1 | ID)
  )
  
  # 调整mediate参数（使用dm1作为中介变量）
  medi_model <- mediate(model1, model2,
                        treat = "Zn_log_scaled", # 直接暴露变量名
                        mediator = "liverhard_log", # 固定为dm1
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
    medi_name = "LSM",
    outcome_name = path
  )
  mediation_sens1 <- rbind(mediation_sens1, df_table)
}

export(mediation_sens1,"result/draft/mediation_sens1_row.xlsx")

mediation_sens1_1 <- mediation_sens1 %>%
  filter(name != "Prop. Mediated") %>%
  mutate(
    Estimate = (10^Estimate - 1) * 100,
    `95% CI Lower` = (10^`95% CI Lower` - 1) * 100,
    `95% CI Upper` = (10^`95% CI Upper` - 1) * 100
  )

mediation_sens1_2 <- mediation_sens1 %>%
  filter(name == "Prop. Mediated")

mediation_sens1_V2 <- bind_rows(mediation_sens1_1, mediation_sens1_2)

mediation_sens1_V2 <- mediation_sens1_V2 %>% process_mediation_results()
export(mediation_sens1_V2, file = "result/draft/mediation_sens1_V2.xlsx")
mediation_sens1_V2 <- import("result/draft/mediation_sens1_V2.xlsx")
mediation_sens1_V3 <- mediation_sens1_V2 %>% 
  mutate(PFDR_pathway = p.adjust(`p-value_Prop. Mediated`,"BH"))
mediation_sens1_V3 %>% filter(PFDR_pathway < 0.05)



## 代谢物中介 ####
library(mediation)
set.seed(2024)
mediation_sens2 <- tibble()
for (medi in c(unique(path_meta_table$Metabolite_ID))) {
  Metabolite_ID_new <- paste0(medi,"_log_scaled")
  
  model1 <- lmer(data = data7, formula = liverhard_log ~ Zn_log_scaled + Age + Gender + BMI + Smoke + Drink + Exercises + Education + diet_score + Hypertension_new + Hyperlipidemia_new + Diabetes_new + eGFR + (1 | ID))
  model2 <- lmer(data = data7, formula = get(Metabolite_ID_new) ~ liverhard_log + Zn_log_scaled + Age + Gender + BMI + Smoke + Drink + Exercises + Education + diet_score + Hypertension_new + Hyperlipidemia_new + Diabetes_new + eGFR + (1 | ID))
  medi_model <- mediate(model1, model2, treat = "Zn_log_scaled", mediator = "liverhard_log", sims = 2000)
  df <- summary(medi_model)
  df_table <- tibble(
    name = c("ACME", "ADE", "Total Effect", "Prop. Mediated"),
    Estimate = c(df$d0, df$z0, df$tau.coef, df$n0),
    `95% CI Lower` = c(df$d0.ci[1], df$z0.ci[1], df$tau.ci[1], df$n0.ci[1]),
    `95% CI Upper` = c(df$d0.ci[2], df$z0.ci[2], df$tau.ci[2], df$n0.ci[2]),
    `p-value` = c(df$d0.p, df$z0.p, df$tau.p, df$n0.p)
  ) %>% mutate(
    exposure_name = "Zn",
    medi_name = "LSM",
    outcome_name = medi
    
  )
  mediation_sens2 <- rbind(mediation_sens2, df_table)
}

export(mediation_sens2,file = "result/draft/mediation_sens2_row.xlsx")

mediation_sens2_1 <- mediation_sens2 %>%
  filter(name != "Prop. Mediated") %>%
  mutate(
    Estimate = (10^Estimate - 1) * 100,
    `95% CI Lower` = (10^`95% CI Lower` - 1) * 100,
    `95% CI Upper` = (10^`95% CI Upper` - 1) * 100
  )

mediation_sens2_2 <- mediation_sens2 %>%
  filter(name == "Prop. Mediated")

mediation_sens2_V2 <- bind_rows(mediation_sens2_1, mediation_sens2_2)
mediation_sens2_V2 <- mediation_sens2_V2 %>% process_mediation_results()

export(mediation_sens2_V2,file = "result/draft/mediation_sens2_V2.xlsx")

## 合并通路和 代谢物 ####
mediation_sens2_V2 <- import("result/draft/mediation_sens2_V2.xlsx")

Table_sens2 <- mediation_sens1_V3 %>% 
  rename(Pathway_name = outcome_name) %>% 
  left_join(path_meta_table, by = "Pathway_name") %>% 
  mutate(outcome_name = Metabolite_ID) %>% 
  left_join(mediation_sens2_V2, by = "outcome_name") %>% 
  dplyr::select(all_of(c("exposure_name.x", "Pathway_name", "ACME (β, 95%CI).x", "p-value_ACME.x", "Prop. Mediated.x",
                         "PFDR_pathway","p-value_Prop. Mediated.x","Metabolite_Name","ACME (β, 95%CI).y", "p-value_ACME.y", 
                         "Prop. Mediated.y","p-value_Prop. Mediated.y")))

Table_sens2 <- Table_sens2 %>% 
  group_by(Pathway_name) %>% 
  mutate(PFDR_metabolite = p.adjust(`p-value_Prop. Mediated.y`, "BH"))

export(Table_sens2, file = "result/Table/Table S4.xlsx")




# 进一步矫正白蛋白 ####
data7$ALB_log <- log(data7$ALB)

LME_table5 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates, collapse = "+"), "+ ALB + (1|ID)")
  model <- lmer(formula = as.formula(formular), data = data7)
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_table5 <- bind_rows(LME_table5, df[2, ])
}

LME_table5 <- LME_table5 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_table5_V2 <- LME_table5 |> LME_percent_format()


plot5 <- LME_table5_V2 |>
  mutate(P_lable = ifelse(PFDR < 0.05, "*", ifelse(PFDR < 0.01, "**", ifelse(PFDR < 0.001, "***", "")))) |>
  mutate(across(
    c(Estimate, LOW_CI, UP_CI),
    ~ as.numeric(.x)
  )) |>
  mutate(Exposure = factor(Exposure, levels = exposures_name)) |>
  ggplot(aes(x = Exposure, y = Estimate)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.8
  ) + # 优化参考线
  geom_errorbar(aes(ymin = LOW_CI, ymax = UP_CI),
                width = 0.2,
                linewidth = 0.8
  ) + # 误差条样式
  geom_point(
    stroke = 1
  ) + # 点样式（填充+边框）
  geom_text(aes(x = Exposure, y = UP_CI + 0.25, label = P_lable), family = "Times New Roman", size = 7, fontface = "bold") +
  labs(
    x = "",
    y = "%Change (95% CI) of LSM",
    title = NULL # 科研图通常不加标题
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + # 新增此行
  theme_classic(base_size = 12) + # 经典无网格主题
  theme(
    text = element_text(family = "Times New Roman"), # 指定字体
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
  )

ggsave(
  filename = "result/Figure/Figure S6.png", plot = plot5, device = "png",
  width = 7.5, height = 3.5, dpi = 300
)



# 加上季节 ####

data7$season.y <- as.numeric(data7$season.y)
summary(data7[,c(covariates,"season.y")])

LME_table5 <- tibble()
for (expo in exposures_log_scaled) {
  formular <- paste0("liverhard_log", "~", expo, "+", paste0(covariates, collapse = "+"), " + season.y + (1|ID)")
  model <- lmer(formula = as.formula(formular), data = data7)
  
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Exposure = expo)
  LME_table5 <- bind_rows(LME_table5, df[2, ])
}

LME_table5 <- LME_table5 |> mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
LME_table5_V2 <- LME_table5 |> LME_percent_format()



plot5 <- LME_table5_V2 |>
  mutate(P_lable = ifelse(PFDR < 0.05, "*", ifelse(PFDR < 0.01, "**", ifelse(PFDR < 0.001, "***", "")))) |>
  mutate(across(
    c(Estimate, LOW_CI, UP_CI),
    ~ as.numeric(.x)
  )) |>
  mutate(Exposure = factor(Exposure, levels = exposures_name)) |>
  ggplot(aes(x = Exposure, y = Estimate)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.8
  ) + # 优化参考线
  geom_errorbar(aes(ymin = LOW_CI, ymax = UP_CI),
                width = 0.2,
                linewidth = 0.8
  ) + # 误差条样式
  geom_point(
    stroke = 1
  ) + # 点样式（填充+边框）
  geom_text(aes(x = Exposure, y = UP_CI + 0.25, label = P_lable), family = "Times New Roman", size = 7, fontface = "bold") +
  labs(
    x = "",
    y = "%Change (95% CI) of LSM",
    title = NULL # 科研图通常不加标题
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + # 新增此行
  theme_classic(base_size = 12) + # 经典无网格主题
  theme(
    text = element_text(family = "Times New Roman"), # 指定字体
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
  )

ggsave(
  filename = "result/Figure/Figure S6.png", plot = plot5, device = "png",
  width = 7.5, height = 3.5, dpi = 300
)
