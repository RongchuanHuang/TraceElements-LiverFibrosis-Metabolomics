# ===========================================================================
# Project: Trace Elements, Metabolomics, and Liver Fibrosis
# Script: Methodological Reference for Statistical Modeling
# Description: This script outlines the core analytical pipeline and model 
#              specifications used in the study. Researchers can adapt this 
#              framework by substituting 'analysis_data' with their own 
#              longitudinal cohort dataset to replicate the methodologies.
# ===========================================================================

# ---------------------------------------------------------------------------
# 0. 加载核心工具包 (Load Core Packages)
# ---------------------------------------------------------------------------
library(tidyverse)
library(lme4)
library(lmerTest)  
library(mgcv)
library(bkmr)
library(qgcomp)
library(mediation)
library(cit)
library(glue)

# ---------------------------------------------------------------------------
# 1. 分析框架与变量集定义 (Analytical Framework & Variable Definitions)
# ---------------------------------------------------------------------------
# 核心协变量集 (Covariates strictly adjusted in all multivariable models)
covariates <- c("Age", "Gender", "BMI", "Smoke", "Drink", "Exercises",
                "Education_new", "diet_score", "Hypertension_new",
                "Hyperlipidemia_new", "Diabetes_new", "eGFR")

# 暴露变量: 对数转换并标准化的微量元素 (Log-transformed & scaled exposures)
exposures_name       <- c("Fe", "Zn", "Cu", "Mn", "Se", "V", "Cr", "Ni", "Mo", "Co")
exposures_log_scaled <- paste0(exposures_name, "_log_scaled")

# 代谢组学特征 (Metabolomic features)
n_metab                <- 20
metabolites_name       <- paste0("Metab", sprintf("%03d", 1:n_metab))
metabolites_log_scaled <- paste0(metabolites_name, "_log_scaled")

# 纵向结局指标 (Longitudinal outcomes: LSM and CAP)
outcomes <- c("liverhard_log", "fat_decay_log")

# ---------------------------------------------------------------------------
# 2. 数据预处理接口 (Data Import & Preprocessing Interface)
#    (需在此处读入用户的真实纵向队列数据)
# ---------------------------------------------------------------------------
# analysis_data <- read.csv("your_longitudinal_cohort_data.csv")

# 确保分类变量已被正确定义为因子 (Ensure categorical variables are factorized)
# factor_vars <- c("Gender", "Smoke", "Drink", "Exercises", "Education_new",
#                  "Hypertension_new", "Hyperlipidemia_new", "Diabetes_new", "ID")
# analysis_data <- analysis_data %>% mutate(across(all_of(factor_vars), as.factor))


# ---------------------------------------------------------------------------
# 3. 线性混合效应模型 (LME: Linear Mixed-Effects Models)
#    评估单一暴露与结局纵向关联 (Exposure-Outcome Associations)
# ---------------------------------------------------------------------------
LME_table1 <- tibble()

for (out in outcomes) {
  for (expo in exposures_log_scaled) {
    # 构建包含随机截距的主效应模型 (Main effect model with random intercept for ID)
    formular <- paste0(out, " ~ ", expo, " + ", paste(covariates, collapse = " + "), " + (1|ID)")
    model <- lmer(formula = as.formula(formular), data = analysis_data)
    
    # 提取系数并整理
    df <- as.data.frame(summary(model)$coefficients)
    df <- df %>% mutate(Exposure = expo, Outcome = out)
    LME_table1 <- bind_rows(LME_table1, df[2, ])   # 提取目标暴露项
  }
}

# 多重比较校正 (FDR adjustment using Benjamini-Hochberg)
LME_results <- LME_table1 %>%
  group_by(Outcome) %>% 
  mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))


# ---------------------------------------------------------------------------
# 4. 广义加性混合模型 (GAMM: Generalized Additive Mixed Models)
#    评估非线性剂量-反应关系 (Non-linear Dose-Response Relationships)
# ---------------------------------------------------------------------------
# 以锌 (Zn) 为例进行平滑样条拟合
expo <- "Zn_log_scaled"

# 剔除极端值以提升样条回归的稳健性 (Exclude 2.5% extreme quantiles)
q_low  <- quantile(analysis_data[[expo]], 0.025, na.rm = TRUE)
q_high <- quantile(analysis_data[[expo]], 0.975, na.rm = TRUE)

data_subset <- analysis_data %>% filter(!!sym(expo) >= q_low & !!sym(expo) <= q_high)

# 构建带有惩罚三次样条和随机效应的公式 (Cubic regression splines with subject RE)
formula_terms <- c(
  paste0("s(", expo, ", bs='cr', k=3)"),   
  covariates,
  "s(ID, bs='re')"                         
)
formula_gam <- reformulate(formula_terms, response = "liverhard_log")

# 采用 REML 方法拟合 (Restricted Maximum Likelihood)
model_gam <- gam(formula_gam, data = data_subset, family = gaussian(), method = "REML")


# ---------------------------------------------------------------------------
# 5. 代谢组学广关联分析 (MWAS: Metabolome-Wide Association Studies)
#    识别暴露驱动的代谢物及其临床结局映射
# ---------------------------------------------------------------------------
MWAS_Zn <- tibble()
MWAS_LSM <- tibble()

for (medi in metabolites_log_scaled) {
  # 路径 A: 暴露至代谢物 (Exposure -> Metabolite)
  formula_A <- reformulate(c("Zn_log_scaled", covariates, "(1|ID)"), response = medi)
  model_A <- lmer(formula = formula_A, data = analysis_data)
  
  df_A <- as.data.frame(summary(model_A)$coefficients) %>% mutate(Exposure = "Zn", Medi = medi)
  MWAS_Zn <- bind_rows(MWAS_Zn, df_A[2, ])
  
  # 路径 B: 代谢物至结局 (Metabolite -> Outcome)
  formula_B <- reformulate(c(medi, covariates, "(1|ID)"), response = "liverhard_log")
  model_B <- lmer(formula = formula_B, data = analysis_data)
  
  df_B <- as.data.frame(summary(model_B)$coefficients) %>% mutate(Medi = medi, Outcome = "LSM")
  MWAS_LSM <- bind_rows(MWAS_LSM, df_B[2, ])
}

# 独立应用 FDR 校正
MWAS_Zn <- MWAS_Zn %>% mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))
MWAS_LSM <- MWAS_LSM %>% mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))


# ---------------------------------------------------------------------------
# 6. 贝叶斯核机器学习 (BKMR: Bayesian Kernel Machine Regression)
#    解析多金属混合暴露模式与变量筛选
# ---------------------------------------------------------------------------
bkmr_expo <- analysis_data[, exposures_log_scaled] %>% as.matrix()
bkmr_cova <- as.matrix(model.matrix(~ . - 1, data = analysis_data[, covariates]))
bkmr_out  <- analysis_data[, "liverhard_log"] %>% as.matrix()
ID_vec    <- analysis_data[, "ID"] %>% as.matrix()

set.seed(2026)
bkmr_model <- kmbayes(
  y = bkmr_out, 
  Z = bkmr_expo, 
  X = bkmr_cova,
  iter = 50000,          # 标准 MCMC 迭代次数 (MCMC iterations)
  id = ID_vec,           # 引入受试者 ID 处理重复测量
  varsel = TRUE          # 启用变量选择 (Enable variable selection)
)

# 提取后验包含概率与整体风险 (Extract PIPs and Overall Risk)
PIP_table <- ExtractPIPs(bkmr_model)
risks.overall <- OverallRiskSummaries(
  fit = bkmr_model, qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5
)


# ---------------------------------------------------------------------------
# 7. 分位数G计算 (Quantile g-computation)
#    评估混合暴露的总体效应幅度与各组分权重
# ---------------------------------------------------------------------------
formular_qg <- formula(paste0(
  "liverhard_log ~ ", paste0(exposures_log_scaled, collapse = " + "), " + ",
  paste0(covariates, collapse = " + ")
))

QG_model <- qgcomp.glm.boot(
  f = formular_qg, 
  data = analysis_data,
  expnms = exposures_log_scaled, 
  q = 4,              # 四分位数转化 (Quartiles)
  B = 1000,           # Bootstrap 重抽样次数 (Bootstrap iterations)
  seed = 2026
)


# ---------------------------------------------------------------------------
# 8. 中介效应评估与因果推断检验 (Mediation & Causal Inference)
#    分析路径: Exposure -> Metabolite -> Outcome
# ---------------------------------------------------------------------------
medi_example <- metabolites_log_scaled[1]

# 8a. 传统中介效应分析 (Counterfactual mediation analysis)
set.seed(2024)
model_med_1 <- lmer(data = analysis_data, formula = as.formula(paste0(
  medi_example, " ~ Zn_log_scaled + ", paste0(covariates, collapse = " + "), " + (1|ID)")))
model_med_2 <- lmer(data = analysis_data, formula = as.formula(paste0(
  "liverhard_log ~ ", medi_example, " + Zn_log_scaled + ", paste0(covariates, collapse = " + "), " + (1|ID)")))

mediation_results <- mediate(
  model_med_1, model_med_2, 
  treat = "Zn_log_scaled", mediator = medi_example, 
  sims = 1000       # Bootstrap 次数
)

# 8b. 因果推断统计检验 (CIT: Causal Inference Test)
L_data <- analysis_data[["Zn_log_scaled"]]
G_data <- analysis_data[[medi_example]]
T_data <- analysis_data[["liverhard_log"]]
C_data <- analysis_data[, covariates]

cit_results <- cit.cp(L_data, G_data, T_data, C_data, rseed = 2026)


# ---------------------------------------------------------------------------
# 9. 敏感性分析与效应修饰 (Sensitivity Analysis: Effect Modification)
#    以年龄作为分层交互因素的检验范式
# ---------------------------------------------------------------------------
analysis_data <- analysis_data %>% mutate(Age_num = as.numeric(as.character(Age)))
LME_interactions <- tibble()

for (expo in exposures_log_scaled) {
  # 构建包含交互项的模型 (Interaction term: Exposure * Age)
  f_inter <- paste0("liverhard_log ~ ", expo, " * Age_num + ",
                    paste0(setdiff(covariates, "Age"), collapse = " + "), " + (1|ID)")
  model_inter <- lmer(as.formula(f_inter), data = analysis_data)
  
  df_inter <- as.data.frame(summary(model_inter)$coefficients) %>% mutate(Exposure = expo)
  LME_interactions <- bind_rows(LME_interactions, tail(df_inter, 1)) # 提取末行的交互项系数
}

LME_interactions <- LME_interactions %>% mutate(PFDR = p.adjust(p = `Pr(>|t|)`, method = "BH"))