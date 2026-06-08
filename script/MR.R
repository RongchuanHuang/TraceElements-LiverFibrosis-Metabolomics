library(TwoSampleMR)
library(dplyr)
library(ggplot2)

# --- 第一步：提取数据 ---
# 以血清锌 (ieu-a-1010) 和肝纤维化 (finn-b-K11_FIBRO) 为例
exp_dat <- extract_instruments(outcomes = "ieu-a-1010")
out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = "finn-b-K11_FIBRO")
dat <- harmonise_data(exp_dat, out_dat)

# --- 第二步：计算 F 统计量 (科研必备) ---
# 公式: F = (beta^2) / (se^2)
dat$F_stat <- (dat$beta.exposure^2) / (dat$se.exposure^2)
mean_F <- mean(dat$F_stat)
cat("Mean F-statistic is:", mean_F, "\n(通常建议 > 10)")

# --- 第三步：多方法 MR 分析 ---
res <- mr(dat)

# --- 第四步：计算 Odds Ratio (OR) 及其 95% CI ---
# 结局是二分类变量（如肝纤维化）时，必须报告 OR
res_or <- generate_odds_ratios(res)
write.csv(res_or, "MR_results_OR.csv") # 保存结果表


# 生成基础图对象
p_scatter <- mr_scatter_plot(res, dat)

# 使用 ggplot2 覆盖外观设置
final_p1 <- p_scatter[[1]] + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Causal Effect of Serum Zinc on Liver Fibrosis",
    x = "SNP effect on Serum Zinc",
    y = "SNP effect on Liver Fibrosis"
  )

# 保存为高分辨率 TIFF
ggsave("Scatter_Plot_Zinc.tiff", final_p1, width = 7, height = 6, dpi = 300)


# 构造适合绘图的数据框
res_or$method <- factor(res_or$method, levels = res_or$method)

final_p2 <- ggplot(res_or, aes(x = method, y = or, ymin = or_lci95, ymax = or_uci95)) +
  geom_pointrange(aes(color = method), size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") + # OR=1 参考线
  coord_flip() + # 横轴纵轴翻转
  theme_classic() +
  labs(x = "Analysis Method", y = "Odds Ratio (95% CI)", title = "Forest Plot of MR Estimates") +
  theme(legend.position = "none")

ggsave("Forest_OR_Zinc.pdf", final_p2, width = 8, height = 4)


# 1. 异质性检验 (Heterogeneity)
het <- mr_heterogeneity(dat)

# 2. 多效性检验 (Pleiotropy - MR-Egger Intercept)
pleio <- mr_pleiotropy_test(dat)

# 3. 汇总成标准表格
sensitivity_table <- data.frame(
  Test = c("IVW Q-statistic", "MR-Egger Intercept", "Intercept P-value"),
  Value = c(het$Q[het$method=="Inverse variance weighted"], pleio$egger_intercept, pleio$pval)
)
print(sensitivity_table)
