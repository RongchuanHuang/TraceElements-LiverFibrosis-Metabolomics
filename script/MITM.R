library(bruceR)
library(tidyverse)

# 导入数据 ####
data7 <- import("data/processed/data7.rds")
MWAS_Zn_V2 <- import("result/draft/MWAS_Zn_V2.xlsx")
MWAS_LSM_V2 <- import("result/draft/MWAS_LSM_V2.xlsx")

# 导出显著的代谢物列表 ####
MWAS_Zn_V2_sig_neg <- MWAS_Zn_V2 %>% filter(PFDR < 0.2) %>% filter(Estimate < 0) %>% filter(!is.na(KEGG))
MWAS_Zn_V2_sig_pos <- MWAS_Zn_V2 %>% filter(PFDR < 0.2) %>% filter(Estimate > 0) %>% filter(!is.na(KEGG))
 
MWAS_LSM_V2_sig_neg <- MWAS_LSM_V2 %>% filter(PFDR < 0.2) %>% filter(Estimate < 0) %>% filter(!is.na(KEGG))
MWAS_LSM_V2_sig_pos <- MWAS_LSM_V2 %>% filter(PFDR < 0.2) %>% filter(Estimate > 0) %>% filter(!is.na(KEGG))


data_MWAS_Zn_sig_neg <- data7 %>% select(sampleID,Zn_log_scaled, MWAS_Zn_V2_sig_neg$Medi) 
colnames(data_MWAS_Zn_sig_neg)[-c(1:2)] <- MWAS_Zn_V2_sig_neg$KEGG


data_MWAS_Zn_sig_pos <- data7 %>% select(sampleID,Zn_log_scaled, MWAS_Zn_V2_sig_pos$Medi)
colnames(data_MWAS_Zn_sig_pos)[-c(1:2)] <- MWAS_Zn_V2_sig_pos$KEGG


data_LSM_V2_sig_neg <- data7 %>% select(sampleID,liverhard_log, MWAS_LSM_V2_sig_neg$Medi)
colnames(data_LSM_V2_sig_neg)[-c(1:2)] <- MWAS_LSM_V2_sig_neg$KEGG


data_LSM_V2_sig_pos <- data7 %>% select(sampleID,liverhard_log, MWAS_LSM_V2_sig_pos$Medi)
colnames(data_LSM_V2_sig_pos)[-c(1:2)] <- MWAS_LSM_V2_sig_pos$KEGG



export(data_MWAS_Zn_sig_neg,file =  "result/draft/data_MWAS_Zn_sig_neg.csv")
export(data_MWAS_Zn_sig_pos,file =  "result/draft/data_MWAS_Zn_sig_pos.csv")
export(data_LSM_V2_sig_neg,file =  "result/draft/data_LSM_V2_sig_neg.csv")
export(data_LSM_V2_sig_pos,file =  "result/draft/data_LSM_V2_sig_pos.csv")


# 读取通路的结果 ####
LSM_pathway <- import("result\\draft\\LSM_pos\\pathway_results.csv") %>% 
  mutate(Derection = "Positive", VAR = "LSM")
Zn_pathway <- import("result\\draft\\Zn_neg\\pathway_results.csv") %>% 
  mutate(Derection = "Negative", VAR = "Zn")

plot_data <- bind_rows(Zn_pathway, LSM_pathway)

plot_data <-
  plot_data %>%
  mutate(
    point_color = if_else(
      VAR == "Zn", "Minerals-MetPath", "LSM-MetPath"
    ),
    VAR = factor(VAR, levels = c("Zn", "LSM")),
    FDR_lable = case_when(
      `Holm adjust` < 0.05 ~ "<0.05", 
      TRUE ~ "≥0.05" 
    )
  ) |>
  mutate(
    V1 = str_replace(V1, "alpha", "α"),
    V1 = str_replace(V1, "beta", "β")
  )
V1_sort <- rev(sort(unique(plot_data$V1)))

plot_data <- plot_data |> mutate(V1 = factor(V1, levels = V1_sort))


# 1. 斑马纹数据框 (保持不变)
bg_data <- data.frame(
  V1 = levels(plot_data$V1),
  # 偶数行灰色，奇数行白色
  bg = rep(c("grey95", "white"), length.out = length(levels(plot_data$V1))) 
)

# 2. 修正后的绘图代码
p3 <- ggplot() +
  # 【修改点在这里】：将 x = 1.5 改为 x = "Zn"，使其与你的 VAR 数据类型匹配
  geom_tile(data = bg_data, aes(x = "Zn", y = V1, fill = bg), 
            width = Inf, height = 1, show.legend = FALSE) +
  scale_fill_identity() +
  # 叠加实际的数据点
  geom_point(data = plot_data, aes(x = VAR, y = V1, size = Impact, color = Derection, shape = FDR_lable), 
             stroke = 0.8) +
  theme_minimal(base_size = 13) +
  scale_shape_manual(
    values = c("<0.05" = 19, "≥0.05" = 1),
    name = "Statistical significance",
    labels = c(
      "<0.05" = expression(italic(P)["Holm adjust"] < 0.05),
      "≥0.05" = expression(italic(P)["Holm adjust"] >= 0.05)
    )
  )  +
  scale_color_manual(
    name = "Direction",
    values = c("Negative" = "#0072B5FF", "Positive" = "#BC3C29FF") # 经典的Lancet红蓝
  ) +
  scale_size_continuous(range = c(2, 6)) +
  labs(x = "", y = "") +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"),
    # 去除自带的网格线，依赖斑马纹
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "right",
    legend.box.background = element_rect(color = "black", linewidth = 0.5)

  )

ggsave(filename = "result/Figure/Figure 3.pdf", plot = p3, device = cairo_pdf,width = 8, height = 9,
       dpi = 300)




