library(tidyverse)
library(bruceR)

data1 <- import("data/row/鹤峰问卷712人次_20250414_终版0527改.sav")
data2 <- import("data/row/四季体检数据712人次_20250414_终版.sav")
data3 <- import("data/row/四季微量营养素数据712人次LOD替换_250414.sav")
data4 <- import("data/row/四季血浆代谢组数据712人次_250414.sav")
data_list <- list(data1, data2, data3, data4)

# 合并数据集
data5 <- reduce(data_list, function(x, y) full_join(x, y, by = c("ID", "sampleID")))

export(data5,file = "data/processed/data5.rds")
