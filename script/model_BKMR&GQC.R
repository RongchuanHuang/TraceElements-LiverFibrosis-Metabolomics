library(tidyverse)
library(bruceR)

data7 <- import("data/processed/data7.rds")

# BKMR ####
# 建模  
library(bkmr)

bkmr_expo <- data7[, exposures_log_scaled] %>% as.matrix()
bkmr_cova <- as.matrix(model.matrix(~ . - 1, data = data7[, covariates]))
bkmr_out <- data7[, "liverhard_log"] %>% as.matrix()
ID <- data7[, "ID"] %>% as.matrix()

set.seed(2026)
bkmr_model <- kmbayes(
  y = bkmr_out, Z = bkmr_expo, X = bkmr_cova, iter = 20000, id = ID,
  varsel = T
)

save(bkmr_model, file =  "result/model/bkmr_model.Rdata")

# 绘图
load("result/model/bkmr_model.Rdata")

## PIP提取 
PIP_table <- ExtractPIPs(bkmr_model)
PIP_table <- PIP_table %>% 
  mutate(variable = str_split_i(variable, "_log", 1),
         PIP = round(PIP, 3))
export(PIP_table, file = "result/draft/PIP table.xlsx")


library(ggsci)
pip_plot <- ggplot(PIP_table, aes(x = reorder(variable, PIP), y = PIP)) +
  geom_col(fill = pal_lancet()(1), width = 0.7) +  # 使用Lancet期刊风格的单色
  geom_text(aes(label = sprintf("%.3f", PIP)), 
            hjust = -0.2, size = 3, family = "Times New Roman") +  # 添加数值标签
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # 扩展y轴范围以容纳标签
  coord_flip() +  # 横向柱形图，便于阅读标签
  labs(title = "(A)",
       x = "Essential trace elements", 
       y = "Posterior Inclusion Probabilities (PIPs)") +
  theme_classic(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = "white")
  )


## overall effect
risks.overall <- OverallRiskSummaries(
  fit = bkmr_model,
  qs = seq(0.25, 0.75, by = 0.05),
  q.fixed = 0.5
)

risks.overall <- risks.overall %>% 
  mutate(lower = est - 1.96 * sd, upper = est + 1.96 * sd) %>% 
  mutate(est = (10^est - 1) * 100,
         lower = (10^lower - 1) * 100,
         upper = (10^upper - 1) * 100,)

bkmr_plot1 <- ggplot(risks.overall, aes(quantile, est, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Percent change of LSM", x = "Quantile of essential trace elements mixtures", title = "(B)") +
  theme_classic() +
  theme(text = element_text(family = "Times New Roman"),
        plot.title = element_text(face = "bold"))

plots <- gridExtra::grid.arrange(pip_plot, bkmr_plot1, nrow = 1)



# GQC ####
library(qgcomp)
formular <- formula(paste0("liverhard_log ~ ", paste0(exposures_log_scaled, collapse = "+"), "+", paste0(covariates, collapse = "+")))
print(formular)
QG_model <- qgcomp.glm.boot(f = formular, data = data7, expnms = exposures_log_scaled, q = 4, B = 2000, seed = 1)
save(QG_model, file = "result/model/QG_model.Rdata")

load("result/model/QG_model.Rdata")

plot(QG_model)
QG_effect <- summary(QG_model)
QG_effect <- QG_effect$coefficients %>% as.data.frame()
QG_effect <- QG_effect %>% mutate(Estimate = (10 ^ Estimate - 1) * 100,
                                  `Lower CI` = (10 ^ `Lower CI` - 1) * 100,
                                  `Upper CI` = (10 ^ `Upper CI` - 1) * 100)

export(QG_effect, file = "result/draft/QG_table.xlsx")


## 提取权重 
coefficients <- GQ_model$fit$coefficients
coefficients <- coefficients[2:11]
names(coefficients) <- str_split_i(names(coefficients), pattern = "_log", i = 1)
QG_weight <- tibble(exposure = names(coefficients), est = coefficients)

QG_weight <- QG_weight %>%
  mutate(
    weight_type = ifelse(est > 0, "positive", "negative")
  ) %>%
  group_by(weight_type) %>%
  mutate(
    total_abs = sum(est),
    relative_weight = est / total_abs,
    relative_weight = ifelse(weight_type == "positive", relative_weight, -relative_weight)
  ) %>% 
  arrange(relative_weight)


stats_text <- paste0(
  "Percent Change: ", round(2.779, 2), "% ",
  "(95% CI: ", round(-1.138, 2), "%, ", round(6.852, 2), "%), ",
  "<i>P</i> = ", round(0.167, 3)
)

library(ggtext)
QG_plot <- ggplot(QG_weight, aes(y = reorder(exposure, relative_weight), x = relative_weight, fill = weight_type)) +
  ggsci::scale_fill_lancet() +
  labs(
    title = "(C)",
    x = "Relative Weight", 
    y = "", 
    fill = "Direction",
    caption = stats_text
  ) +
  # 添加自定义网格线
  geom_segment(
    aes(x = -Inf, xend = 0, y = reorder(exposure, relative_weight), yend = reorder(exposure, relative_weight)),
    color = alpha("gray", 0.5),
    data = QG_weight
  ) +
  geom_col() + 
  # 确保X轴从0开始
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    text = element_text(family = "Times New Roman", colour = "black"),
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    
    # 关键修改 1: X轴标题居中
    axis.title.x = element_text(hjust = 0.7),
    plot.caption = element_markdown(hjust = 0.1, size = 9, 
                                     color = "grey30", face = "bold")
  )


# 拼图
plots <- gridExtra::grid.arrange(
  pip_plot,
  bkmr_plot1,
  QG_plot,
  nrow = 2
)


ggsave(filename = "result/Figure/Figure S4.png", plot = plots, device = "png",
       width = 8, height = 8, dpi = 300)

