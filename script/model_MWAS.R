library(tidyverse)
library(bruceR)

data7 <- import("data/processed/data7.rds")
code_book <- import("D:\\科研\\项目_ing\\恩施\\TraceElements-LiverFibrosis-Metabolomics\\data\\row\\代谢物ID_名称.xlsx")

# 拟合模型 ####
# Zn-代谢组 
summary(data7[,covariates])

MWAS_Zn <- tibble()
for (medi in metabolites_log_scaled) {
    formula_terms <- c("Zn_log_scaled", covariates, "(1|ID)")
    formular <- reformulate(formula_terms, response = medi)
    print(formular)
    model <- lmer(formula = formular, data = data7)
    # 2. 提取
    df <- summary(model)
    df <- as.data.frame(df$coefficients)
    df <- df |> mutate(Exposure = "Zn", Medi = medi)
    MWAS_Zn <- bind_rows(MWAS_Zn, df[2, ])
}

MWAS_Zn_V2 <- MWAS_Zn |>
  mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH")) %>% 
  mutate(Metabolite_ID = str_split_i(Medi, pattern = "_", 1)) %>% 
  left_join(code_book, by = "Metabolite_ID") %>% 
  mutate(LOW_CI = Estimate - 1.96*`Std. Error`,
         UP_CI = Estimate + 1.96*`Std. Error`) %>% 
  mutate(Estimate_per = (10 ^ Estimate - 1) * 100,
         LOW_CI_per = (10 ^ LOW_CI - 1) * 100,
         UP_CI_per = (10 ^ UP_CI - 1) * 100)
export(MWAS_Zn_V2, file = "result/draft/MWAS_Zn_V2.xlsx")

# 代谢组-肝脏弹性 
MWAS_LSM <- tibble()
for (medi in metabolites_log_scaled) {
  # 1. 构建公式
  formula_terms <- c(medi, covariates, "(1|ID)")
  formular <- reformulate(formula_terms, response = "liverhard_log")
  print(formular)
  model <- lmer(formula = formular, data = data7)
  # 2. 提取
  df <- summary(model)
  df <- as.data.frame(df$coefficients)
  df <- df |> mutate(Medi = medi, Out = "LSM")
  MWAS_LSM <- bind_rows(MWAS_LSM, df[2, ])
}

MWAS_LSM_V2 <- MWAS_LSM |> 
  mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH")) %>% 
  mutate(Metabolite_ID = str_split_i(Medi, pattern = "_", 1)) %>% 
  left_join(code_book, by = "Metabolite_ID") %>% 
  mutate(LOW_CI = Estimate - 1.96*`Std. Error`,
         UP_CI = Estimate + 1.96*`Std. Error`) %>% 
  mutate(Estimate_per = (10 ^ Estimate - 1) * 100,
         LOW_CI_per = (10 ^ LOW_CI - 1) * 100,
         UP_CI_per = (10 ^ UP_CI - 1) * 100)

export(MWAS_LSM_V2, file = "result/draft/MWAS_LSM_V2.xlsx")


# 火山图 ####
MWAS_Zn_V2 <- import("result/draft/MWAS_Zn_V2.xlsx")
MWAS_LSM_V2 <- import("result/draft/MWAS_LSM_V2.xlsx")

dataplot <- MWAS_Zn_V2 %>%
  mutate(
    point_color = factor(
      ifelse(PFDR > 0.2, "Non-sig",
             ifelse(Estimate < 0, "Negative", "Positive")
      ),
      levels = c("Negative", "Positive", "Non-sig") # 确保图例顺序
    )
  )

Positive_num <- dataplot |>
  filter(point_color == "Positive") |>
  nrow()
Negative_num <- dataplot |>
  filter(point_color == "Negative") |>
  nrow()

plot_MWAS1 <- dataplot %>%
  ggplot(aes(x = Estimate, y = -log10(`Pr(>|t|)`), colour = point_color)) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Significance",
    values = c(
      "Negative" = "darkblue",
      "Positive" = "darkred",
      "Non-sig" = "#999999"
    ),
    drop = FALSE
  ) +
  labs(
    x = "β",
    y = expression(-log[10](italic(P))),
    title = "(A) Zn"
  ) +
  theme_classic() +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Positive: ", Positive_num),
           color = "darkred",
           hjust = 1.1, vjust = 1.5,
           size = 4, family = "Times New Roman", fontface = "bold"
  ) +
  annotate("text",
           x = -Inf, y = Inf,
           label = paste0("Negative: ", Negative_num),
           color = "darkblue",
           hjust = -0.1, vjust = 1.5,
           size = 4, family = "Times New Roman", fontface = "bold"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(family = "Times New Roman"),
    text = element_text(family = "Times New Roman")
  ) +
  scale_y_continuous(limits = c(NA, 90),labels = scales::number_format(accuracy = 0.1))



dataplot <- MWAS_LSM_V2 %>%
  mutate(
    point_color = factor(
      ifelse(PFDR > 0.2, "Non-sig",
             ifelse(Estimate < 0, "Negative", "Positive")
      ),
      levels = c("Negative", "Positive", "Non-sig") # 确保图例顺序
    )
  )

Positive_num <- dataplot |>
  filter(point_color == "Positive") |>
  nrow()
Negative_num <- dataplot |>
  filter(point_color == "Negative") |>
  nrow()

plot_MWAS2 <- dataplot %>%
  ggplot(aes(x = Estimate, y = -log10(`Pr(>|t|)`), colour = point_color)) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Significance",
    values = c(
      "Negative" = "darkblue",
      "Positive" = "darkred",
      "Non-sig" = "#999999"
    ),
    drop = FALSE
  ) +
  labs(
    x = "β",
    y = expression(-log[10](italic(P))),
    title = "(B) LSM"
  ) +
  theme_classic() +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Positive: ", Positive_num),
           color = "darkred",
           hjust = 1.1, vjust = 1.5,
           size = 4, family = "Times New Roman", fontface = "bold"
  ) +
  annotate("text",
           x = -Inf, y = Inf,
           label = paste0("Negative: ", Negative_num),
           color = "darkblue",
           hjust = -0.1, vjust = 1.5,
           size = 4, family = "Times New Roman", fontface = "bold"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(family = "Times New Roman"),
    text = element_text(family = "Times New Roman")
  ) +
  scale_y_continuous(limits = c(NA, 5),labels = scales::number_format(accuracy = 0.1))

plots <- gridExtra::grid.arrange(plot_MWAS1, plot_MWAS2, nrow = 1)
ggsave(filename = "result/Figure/Figure 2-1.pdf", plot = plots, device = cairo_pdf,width = 8, height = 4, dpi = 300)  

# Class分布 ####
library(dplyr)

df_zn <- MWAS_Zn_V2 %>%
  select(Metabolite_ID, Class, PFDR_zn = PFDR)

df_lsm <- MWAS_LSM_V2 %>%
  select(Metabolite_ID, Class, PFDR_lsm = PFDR)

df_merged <- inner_join(df_zn, df_lsm, by = c("Metabolite_ID","Class"))

# 4. 生成作图需要的数据集
df_plot <- df_merged %>%
  # 处理缺失的分类：你可以选择删掉 filter(!is.na(Class))，或者把 NA 替换为 "Unknown"
  mutate(Class = ifelse(is.na(Class), "Unclassified metabolites", as.character(Class))) %>%
  # 判定显著性
  mutate(
    zn_sig = PFDR_zn < 0.2,
    lsm_sig = PFDR_lsm < 0.2,
    Group = case_when(
      zn_sig & !lsm_sig ~ "Only associated with Zn",
      !zn_sig & lsm_sig ~ "Only associated with LSM",
      zn_sig & lsm_sig ~ "Both associated with Zn and LSM",
      TRUE ~ "Not_Sig"
    )
  ) %>%
  # 过滤掉都不显著的代谢物，只保留有意义的数据用于作图
  filter(Group != "Not_Sig")

library(circlize)
library(dplyr)

chord_data <- table(df_plot$Group, df_plot$Class)

# ========== 用 base R 生成 72 种颜色 ==========
group_names <- rownames(chord_data)
class_names <- colnames(chord_data)
all_sectors <- c(group_names, class_names)
n_all <- length(all_sectors)

# 在色相环上均匀取 72 种颜色（c=80 饱和度适中，l=65 亮度适中）
grid_col <- setNames(
  hcl(h = seq(0, 360, length.out = n_all + 1)[-1], c = 80, l = 65),
  all_sectors
)

# 单独控制三个特殊扇区的颜色
special_sectors <- c("Only associated with Zn",
                     "Only associated with LSM",
                     "Both associated with Zn and LSM")
special_colors <- c("#FF7F00",   # 橙色
                    "#377EB8",   # 蓝色
                    "#E41A1C")   # 红色
names(special_colors) <- special_sectors

for (s in special_sectors) {
  if (s %in% names(grid_col)) {
    grid_col[s] <- special_colors[s]
  }
}

# 图形输出
pdf(file = "result/Figure/Figure 2-2.pdf", 
    width = 7, height = 7)

circos.clear()
par(mar = c(1, 1, 1, 1))   # 请先定义 my_font

chordDiagram(
  chord_data,
  transparency = 0.5,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.25), 
  link.lwd = 1,
  link.lty = 1,
  link.border = NA,
  grid.col = grid_col
)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(
    x = CELL_META$xcenter, 
    y = CELL_META$ylim[1] + 0.1, 
    labels = CELL_META$sector.index, 
    facing = "clockwise", 
    niceFacing = TRUE, 
    adj = c(0, 0.5), 
    cex = 0.35            # 字体已调小
  )
}, bg.border = NA) 

dev.off()




library(dplyr)
library(circlize)

# 1. 提取数据（注意增加 Estimated）
df_zn <- MWAS_Zn_V2 %>%
  select(Metabolite_ID, Class, PFDR_zn = PFDR, Estimated_zn = Estimate)

df_lsm <- MWAS_LSM_V2 %>%
  select(Metabolite_ID, Class, PFDR_lsm = PFDR, Estimated_lsm = Estimate)

# 2. 合并
df_merged <- inner_join(df_zn, df_lsm, by = c("Metabolite_ID", "Class"))

# 3. 生成四个重叠分组并统计与 Class 的交叉表
df_plot <- df_merged %>%
  mutate(Class = ifelse(is.na(Class), "Unclassified metabolites", as.character(Class))) %>%
  mutate(
    Zn_pos  = PFDR_zn < 0.2 & Estimated_zn > 0,
    Zn_neg  = PFDR_zn < 0.2 & Estimated_zn < 0,
    LSM_pos = PFDR_lsm < 0.2 & Estimated_lsm > 0,
    LSM_neg = PFDR_lsm < 0.2 & Estimated_lsm < 0
  )

# 转换为长格式，计数
df_long <- df_plot %>%
  tidyr::pivot_longer(
    cols = c(Zn_pos, Zn_neg, LSM_pos, LSM_neg),
    names_to = "Group",
    values_to = "is_sig"
  ) %>%
  filter(is_sig) %>%
  mutate(Group = recode(Group,
                        Zn_pos  = "Zn positive",
                        Zn_neg  = "Zn negative",
                        LSM_pos = "LSM positive",
                        LSM_neg = "LSM negative"
  )) %>%
  count(Group, Class)

# 生成矩阵
chord_data <- as.data.frame(df_long) %>%
  tidyr::pivot_wider(names_from = Class, values_from = n, values_fill = 0) %>%
  tibble::column_to_rownames("Group") %>%
  as.matrix()

# ---------- 按代谢物类别百分比从大到小排序 ----------
class_order <- names(sort(colSums(chord_data), decreasing = TRUE))
chord_data <- chord_data[, class_order, drop = FALSE]

# 4. 颜色定义
group_names <- rownames(chord_data)
class_names <- colnames(chord_data)
all_sectors <- c(group_names, class_names)
n_all <- length(all_sectors)

# 基础颜色（色相环均匀分布）
grid_col <- setNames(
  hcl(h = seq(0, 360, length.out = n_all + 1)[-1], c = 80, l = 65),
  all_sectors
)

# ---------- 更换为高对比度的特殊分组颜色 ----------
# Zn 正：深红， Zn 负：亮橙
# LSM 正：深蓝， LSM 负：青蓝
special_sectors <- c("Zn positive", "Zn negative", "LSM positive", "LSM negative")
special_colors <- c("#42B540FF",  # 深红，Zn 正
                    "#0099B4FF",  # 亮橙，Zn 负
                    "#925E9FFF",  # 深蓝，LSM 正
                    "#FDAF91FF")  # 青蓝，LSM 负

names(special_colors) <- special_sectors

for (s in special_sectors) {
  if (s %in% names(grid_col)) {
    grid_col[s] <- special_colors[s]
  }
}

# 5. 计算扇区百分比
total_flow <- sum(chord_data)
sector_flow <- c(rowSums(chord_data), colSums(chord_data))
names(sector_flow) <- all_sectors
sector_perc <- sector_flow / total_flow * 100

# 6. 绘图 —— 设置全局字体为 Times New Roman
pdf(file = "result/Figure/Figure 2-2.pdf", 
    width = 7.5, height = 7,
    family = "serif")   # serif 在多数系统中映射为 Times New Roman
par(family = "serif")

circos.clear()
par(mar = c(1, 1, 1, 1))

chordDiagram(
  chord_data,
  transparency = 0.5,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.25),
  link.lwd = 1,
  link.lty = 1,
  link.border = NA,
  grid.col = grid_col
)

# 扇区标签（自动继承 par 中的 serif 字体）
circos.track(track.index = 1, panel.fun = function(x, y) {
  sec <- CELL_META$sector.index
  pct <- round(sector_perc[sec], 1)
  circos.text(
    x = CELL_META$xcenter,
    y = CELL_META$ylim[1] + 0.1,
    labels = paste0(sec, " (", pct, "%)"),
    facing = "clockwise",
    niceFacing = TRUE,
    adj = c(0, 0.5),
    cex = 0.35
  )
}, bg.border = NA)

dev.off()