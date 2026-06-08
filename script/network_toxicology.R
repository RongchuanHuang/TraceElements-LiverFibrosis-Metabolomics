library(tidyverse)
library(bruceR)

# 网络毒理学 ####
## 导入相关基因 ####

Zn_genes_CTD <- import("data/row/CTD_D015032_genes_20260511051111_Zn.csv")
Zn_genes_disgenet <- import("data/row/search_CHEMBL1201279_Zn/search_result_CHEMBL1201279.xlsx")
Zn_genes_genecard <- import("data/row/GeneCards-SearchResults_Zn.csv")

LSM_genes_genecard <- import("data/row/GeneCards-SearchResults-liver fibrosis.csv")

## 质控 ####
Zn_genes_CTD <- Zn_genes_CTD %>% filter(`Interaction Count` > 3)
Zn_genes_genecard <- Zn_genes_genecard %>% filter(`Relevance score` > 5 * median(`Relevance score`, na.rm = TRUE))
Zn_genes_disgenet <- Zn_genes_disgenet %>% filter(score > 0.1)



LSM_genes_genecard <- LSM_genes_genecard %>% filter(`Relevance score` > 5 * median(`Relevance score`, na.rm = TRUE))

## 寻找重叠基因 ####
Zn_genes <- unique(c(Zn_genes_CTD$`Gene Symbol`, Zn_genes_genecard$`Gene Symbol`,Zn_genes_disgenet$Gene))
LSM_genes <- unique(LSM_genes_genecard$`Gene Symbol`)

overlap_genes <- intersect(Zn_genes, LSM_genes)

export(as.data.frame(overlap_genes),file = "result/draft/overlap_gene.csv")

## 重叠基因绘图
library(VennDiagram)

# 2. 打开 PDF 设备（现在路径存在了）
pdf("result/Figure/Figure 5-1-v1.pdf", width = 5, height = 5)

# 3. 绘制韦恩图
venn.diagram(
  x = list("Zn" = Zn_genes, "LSM" = LSM_genes),
  category.names = c("Zn", "LSM"),
  filename = NULL,               # 画到已打开的 PDF 设备上
  lwd = 2,
  col = c("#E41A1C", "#377EB8"),
  fill = c(alpha("#E41A1C", 0.3), alpha("#377EB8", 0.3)),
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = 0.05,
  cat.fontface = "bold",
  cex = 1.2
)

# 4. 关闭设备，保存文件
dev.off()


## GO富集分析 ####
library(clusterProfiler)
library(org.Hs.eg.db) # 人类基因注释数据库
library(enrichplot)
library(DOSE)
library(tidyverse)
library(ggrepel)
library(viridis)
library(ggnewscale)
library(RColorBrewer)

gene_ids <- bitr(
  overlap_genes,
  fromType = "SYMBOL", # 输入类型为基因符号
  toType = c("ENTREZID", "ENSEMBL"), # 转换为Entrez ID和Ensembl ID
  OrgDb = org.Hs.eg.db # 人类基因数据库
)

entrez_ids <- gene_ids$ENTREZID

### 生物过程 (Biological Process) ####
go_bp <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP", # 生物过程
  pvalueCutoff = 0.05, # p值阈值
  qvalueCutoff = 0.05, # q值阈值
  pAdjustMethod = "BH", # 多重检验校正方法
  readable = TRUE # 将Entrez ID转换为基因符号
)

### 细胞组分 (Cellular Component) ####
go_cc <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC", # 细胞组分
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE
)

### 分子功能 (Molecular Function) ####
go_mf <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF", # 分子功能
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE
)

### 可视化 ####
# 设置输出目录

lancet_colors <- c(
  "BP" = "#00468BFF",  # 深蓝色
  "CC" = "#ED0000FF",  # 红色
  "MF" = "#42B540FF"   # 绿色
)

mf_results <- head(go_mf@result[order(go_mf@result$p.adjust), ], 10) %>% mutate(color_c = "MF")
bp_results <- head(go_bp@result[order(go_bp@result$p.adjust), ], 10) %>% mutate(color_c = "BP")
cc_results <- head(go_cc@result[order(go_cc@result$p.adjust), ], 10) %>% mutate(color_c = "CC")

GO_results <- bind_rows(mf_results, bp_results, cc_results)

GO_results <- GO_results %>%
  mutate(
    GeneRatio_num = as.numeric(sub("/\\d+", "", GeneRatio)) /
      as.numeric(sub("\\d+/", "", GeneRatio)) * 100,
    Description = str_to_sentence(Description),
    GeneRatio_label = paste0(round(GeneRatio_num, 1), "%"),
    # 按类别排序：先按类别排序，再按GeneRatio排序
    color_c = factor(color_c, levels = c("BP", "CC", "MF")),
    # 添加类别标签到描述中，便于识别
    Description_with_category = paste0(Description, " (", color_c, ")")
  ) %>% 
  mutate(Description = gsub("\\bDna\\b", "DNA", Description),
         Description = gsub("\\bRna\\b", "RNA", Description))

x_max <- max(GO_results$GeneRatio_num, na.rm = TRUE) * 1.2

# 绘制图形
# 水平分面排列
p_pathname <- ggplot(
  GO_results,
  aes(
    y = reorder(Description, GeneRatio_num),
    x = GeneRatio_num,
    fill = color_c
  )
) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
  geom_text(
    aes(label = GeneRatio_label),
    hjust = -0.2,
    size = 3.2,
    color = "black",
    fontface = "bold",
    family = "Times New Roman"
  ) +
  scale_fill_manual(
    name = "GO Category",
    values = lancet_colors,
    labels = c(
      "BP" = "Biological Process",
      "CC" = "Cellular Component", 
      "MF" = "Molecular Function"
    )
  ) +
  labs(
    x = "Gene Ratio (%)",
    y = ""
  ) +
  scale_x_continuous(
    limits = c(0, x_max),
    expand = expansion(mult = c(0, 0.1))
  ) +
  # 水平分面
  facet_wrap(~ color_c, ncol = 1, scales = "free_y") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "gray"),
    axis.line = element_line(color = "gray", linewidth = 0.5),
    axis.ticks = element_line(color = "gray"),
    text = element_text(family = "Times New Roman"),
    legend.position = "bottom",  # 水平排列时隐藏图例
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3, linetype = "dashed"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_cartesian(clip = "off")


p_GOcode <- ggplot(
  GO_results,
  aes(
    y = reorder(ID, GeneRatio_num),
    x = GeneRatio_num,
    fill = color_c
  )
) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
  geom_text(
    aes(label = GeneRatio_label),
    hjust = -0.2,
    size = 3.2,
    color = "black",
    fontface = "bold",
    family = "Times New Roman"
  ) +
  scale_fill_manual(
    name = "GO Category",
    values = lancet_colors,
    labels = c(
      "BP" = "Biological Process",
      "CC" = "Cellular Component", 
      "MF" = "Molecular Function"
    )
  ) +
  labs(
    x = "Gene Ratio (%)",
    y = ""
  ) +
  scale_x_continuous(
    limits = c(0, x_max),
    expand = expansion(mult = c(0, 0.1))
  ) +
  # 水平分面
  facet_wrap(~ color_c, ncol = 1, scales = "free_y") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "gray"),
    axis.line = element_line(color = "gray", linewidth = 0.5),
    axis.ticks = element_line(color = "gray"),
    text = element_text(family = "Times New Roman"),
    legend.position = "bottom",  # 水平排列时隐藏图例
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3, linetype = "dashed"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_cartesian(clip = "off")

ggsave(filename = "result/Figure/Figure 5-2.pdf", plot = p_pathname, device = cairo_pdf,
       height = 7.5, width = 10, dpi = 300)


## 绘制表格 
table_genes <- GO_results %>% dplyr::select(color_c,Description,geneID)




## KEGG ####
library(clusterProfiler)
library(org.Hs.eg.db)
options(timeout = 600)
kegg_enrich <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa", # 人类
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05, # 可以适当放宽阈值
  minGSSize = 5, # 减小最小基因集大小
  maxGSSize = 500
)

kegg_data <- kegg_enrich@result


# 按照Count降序排列数据
kegg_data_sorted <- kegg_data %>%
  mutate(
    GeneRatio_num = as.numeric(sub("/\\d+", "", GeneRatio)) /
      as.numeric(sub("\\d+/", "", GeneRatio)) * 100,
    log10_padj = -log10(p.adjust),
    Significance = cut(p.adjust,
                       breaks = c(0, 0.001, 0.01, 0.05, 1),
                       labels = c("*** p<0.001", "** p<0.01", "* p<0.05", "NS"),
                       include.lowest = TRUE
    )
  ) %>%
  # 按Count降序排列
  arrange(desc(Count)) %>%
  # 将Description转换为因子，按Count降序排序
  mutate(Description = fct_reorder(Description, Count))

plot_data <- head(kegg_data_sorted,10)

library(scales)
# 创建符合SCI要求的散点图
p_sci <-
  ggplot() +
  # 1. 添加散点
  geom_point(
    data = plot_data,
    aes(
      x = GeneRatio_num,
      y = Description,
      size = Count,
      color = log10_padj
    ),
    shape = 16,
    stroke = 0.5
  ) +
  
  # 2. 添加百分比标签
  geom_text(
    data = plot_data,
    aes(
      x = GeneRatio_num,
      y = Description,
      label = sprintf("%.1f%%", GeneRatio_num)
    ),
    size = 3,
    hjust = -0.5,
    color = "black",
    family = "Times New Roman"
  ) +
  
  # 3. 刻度设置
  scale_x_continuous(
    name = "Gene Ratio (%)",
    limits = c(0, max(plot_data$GeneRatio_num) * 1.25),
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  # 4. 大小范围
  scale_size_continuous(
    name = "Gene Count",
    range = c(3, 8),
    breaks = pretty_breaks(n = 4)(plot_data$Count)
  ) +
  
  # 5. 颜色映射
  scale_color_gradient(
    name = expression(-log[10](bolditalic(P)[FDR])),
    low = "steelblue",
    high = "firebrick",
    guide = guide_colorbar(
      order = 1,
      barwidth = unit(0.5, "cm"),
      barheight = unit(3, "cm")
    )
  ) +
  
  # 6. 主题设置（符合SCI要求）
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 11,face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.spacing.y = unit(0.2, "cm"),
    legend.key = element_blank(),
    text = element_text(family = "Times New Roman", colour = "black")
  ) +
  
  # 8. 标签
  labs(
    title = "",
    y = NULL,
    x = "Gene Ratio (%)"
  )

# 显示图表
print(p_sci)


# 保存为PDF格式
ggsave(
  filename = paste0("result/Figure/Figure 5-3.pdf"),
  plot = p_sci,
  width = 8,
  height = 6,
  device = cairo_pdf
)


## 通路和基因相对应
table_df <- plot_data %>% 
  mutate(color_c = "kegg") %>% 
  dplyr::select(color_c, Description, geneID)

entrez_all <- table_df$geneID %>%
  strsplit("/") %>%
  unlist() %>%
  unique()

id_mapping <- bitr(entrez_all,
                   fromType = "ENTREZID",
                   toType = "SYMBOL",
                   OrgDb = org.Hs.eg.db)

convert_ids <- function(id_string, mapping) {
  ids <- unlist(strsplit(id_string, "/"))
  symbols <- mapping$SYMBOL[match(ids, mapping$ENTREZID)]
  symbols[is.na(symbols)] <- ids[is.na(symbols)]  # 未匹配的保留原ID
  paste(symbols, collapse = "/")
}

# 添加新列
table_df <- table_df %>%
  mutate(geneSymbol = sapply(geneID, convert_ids, mapping = id_mapping)) %>% 
  dplyr::select(color_c, Description, geneSymbol) 

table_df$geneID <- table_df$geneSymbol
table_df <- table_df[,c(1,2,4)]

table_genes <- bind_rows(table_genes, table_df)
export(table_genes, file = "result/Table/Table S2.xlsx")




## 图####
# hub 基因筛选 
library(bruceR)
MCC_df <- read_csv("data/row/string_interactions_short_260526.tsv_MCC_top10 default node.csv")
Degree_df <- read_csv("data/row/string_interactions_short_260526.tsv_Degree_top10 default node.csv")
Betweenness_df <- read_csv("data/row/string_interactions_short_260526.tsv_Betweenness_top10 default node.csv")


# 提取基因名列（如果列名不同，请将下面的 ‘name‘ 替换为实际的列名，例如 ‘node‘ 或 ‘Gene‘）
MCC_genes <- MCC_df$name
Degree_genes <- Degree_df$name
Betweenness_genes <- Betweenness_df$name

# 计算三个列表的交集
common_genes <- Reduce(intersect, list(MCC_genes, Degree_genes, Betweenness_genes))
export(as.data.frame(common_genes), file = "D:\\科研\\项目_ing\\恩施\\肝弹\\网络毒理学\\PPI\\common_genes.csv")

all_unique_genes <- unique(c(MCC_genes, Degree_genes, Betweenness_genes))

# 创建一个二进制矩阵，行是基因，列是算法集合
binary_matrix <- data.frame(
  row.names = all_unique_genes, # 行名为基因
  MCC = as.numeric(all_unique_genes %in% MCC_genes),
  Degree = as.numeric(all_unique_genes %in% Degree_genes),
  Betweenness = as.numeric(all_unique_genes %in% Betweenness_genes)
)
library(rio)
library(UpSetR)

pdf(file = "result/Figure/Figure 5-5.pdf", 
    width = 6, 
    height = 5)
# 绘制图形
upset(binary_matrix,
      sets = c("MCC", "Degree", "Betweenness"),
      sets.bar.color = c("#00468B", "#ED0000", "#42B540"),
      main.bar.color = "darkgrey",
      order.by = "freq",
      text.scale = c(1.2, 1.2, 1, 1, 1.5, 1.2),
      keep.order = TRUE)
dev.off() # 关键！关闭图形设备，完成文件写入


