library(bruceR)
library(tidyverse)

data7 <- import("data/processed/data7.rds")
# GAMM ####
library(mgcv) # 用于GAM模

# 1. 设定当前分析的暴露变量
expo <- "Zn_log_scaled"
expo_x <- str_split(expo, "_", simplify = TRUE)[1, 1] # 提取前缀，例如 "Zn"

# 判断单位
units <- case_when(
  expo_x %in% c("Zn", "Rb") ~ "mg/L",
  expo_x %in% c("V", "Cr") ~ "ug/L",
  TRUE ~ "ug/L"
)

# 2. 数据清洗与截断 
required_vars <- c("liverhard_log", expo, covariates, "ID")

# 计算 2.5% 和 97.5% 分位数，去除异常值
q_low <- quantile(data7[[expo]], 0.025, na.rm = TRUE)
q_high <- quantile(data7[[expo]], 0.975, na.rm = TRUE)

data_subset <- data7 |>
  select(all_of(required_vars)) |>
  filter(!!sym(expo) >= q_low & !!sym(expo) <= q_high)

# 3. 构建公式并拟合 GAM 模型
formula_terms <- c(
  paste0("s(", expo, ", bs='cr', k=3)"),
  covariates,
  "s(ID, bs='re')" # 随机效应
)
formula_gam <- reformulate(formula_terms, response = "liverhard_log")

model <- gam(formula_gam, data = data_subset, family = gaussian(), method = "REML")

# 4. 提取非线性 P 值和 EDF（估计自由度）
summary_model <- summary(model)
smooth_p <- summary_model$s.pv[1]
smooth_edf <- summary_model$edf[1]

# 格式化 P 值和 EDF
p_text <- ifelse(smooth_p < 0.001, "< 0.001", paste0("== ", sprintf("%.3f", smooth_p)))
edf_text <- sprintf("%.2f", smooth_edf)

# 组合标签表达式 (支持 parse = TRUE，渲染成公式)
# 效果类似：EDF == 1.25, P < 0.001
anno_label <- sprintf('bolditalic(edf) == %s ~ "," ~ bolditalic(P) %s', edf_text, p_text)

# 5. 提取预测值
pred <- predict.gam(model, type = "terms", se.fit = TRUE)
expected_col <- paste0("s(", expo, ")")

if (!expected_col %in% colnames(pred$fit)) {
  stop(paste("错误：预测项未找到 -", expected_col))
}

mfit <- pred$fit[, expected_col]
sfit <- pred$se.fit[, expected_col]

# 6. 调整拟合值居中 (表示相对变化)
mfit_centered <- mfit + median(model$fitted.values) - median(mfit)

# 7. 构建绘图数据框并计算 95% 置信区间
dat <- tibble(
  expo_var = data_subset[[expo]],
  mfit = mfit_centered,
  sfit = sfit
) |>
  mutate(
    y.low = mfit - 1.96 * sfit,
    y.upp = mfit + 1.96 * sfit
  )

# 重命名第一列以匹配绘图代码的变量名
names(dat)[1] <- expo

# 8. 绘制图形
p_gam <- ggplot(dat, aes(x = !!sym(expo), y = mfit)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = y.low, ymax = y.upp), alpha = 0.2, fill = "#6BAED6") +
  geom_rug(aes(x = !!sym(expo)), sides = "b", color = "black", alpha = 0.5) +
  labs(
    x = bquote(paste(.(expo_x), " (", .(units), ", log"[10] * " and z-scores transformed)")),
    y = bquote(paste("Changes of LSM (kPa, log"[10] * " transformed)")),
    title = ""
  ) +
  theme_classic() +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  # 添加 EDF 和 P 值的文本
  annotate("text",
           x = -Inf, y = Inf,
           label = anno_label,
           hjust = -0.1, vjust = 1.5, size = 5,
           parse = TRUE, family = "Times New Roman"
  ) +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.tag = element_text(face = "bold", size = 14), # (A) 标签的格式
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(color = "black", size = 12),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5)
  )

# 查看生成的图像
print(p_gam)

# 9. 导出高分辨率图像
ggsave(
  filename = "result/Figure/Figure S3.png", 
  plot = p_gam, 
  device = "png", 
  width = 5, height = 4, dpi = 300
)
