

## 饮食转化函数 ####
diet_function <- function(x) {
  case_when(
    x == 0 ~ 0,
    x %in% 1:4 ~ x / 4,
    x %in% 5:10 ~ x - 4,
    x %in% 11:14 ~ (x - 10) * 7,
    x == 15 ~ 15,
    TRUE ~ 0
  )
}

## 自定义众数函数（处理因子/字符）
get_mode <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  ux <- unique(na.omit(x))
  if (length(ux) == 0) {
    return(NA)
  }
  ux[which.max(tabulate(match(x, ux)))]
}



## 百分比转化函数 ####
LME_percent_format <- function(df) {
  df %>%
    mutate(
      LOW_CI = Estimate - 1.96 * `Std. Error`,
      UP_CI = Estimate + 1.96 * `Std. Error`
    ) %>%
    mutate(across(
      c(Estimate, LOW_CI, UP_CI),
      ~ (10^.x - 1) * 100
    )) %>%
    mutate(Exposure = str_split_i(Exposure, "_", 1)) %>%
    mutate(
      across(c(Estimate, LOW_CI, UP_CI), ~ sprintf("%.2f", .x)),
      `Pr(>|t|)` = sprintf("%.3f", `Pr(>|t|)`),
      PFDR = sprintf("%.3f", PFDR)
    ) %>%
    mutate(`Percent change (95%CI)` = glue::glue("{Estimate} ({LOW_CI}, {UP_CI})"))
}



## 中介分析结果处理函数 ####
process_mediation_results <- function(mediation_table, effect_types = c("ACME", "ADE", "Total Effect", "Prop. Mediated"),
                                      id_vars = c("exposure_name", "medi_name", "outcome_name"),
                                      convert_effects = c("ACME", "ADE", "Total Effect")) {
  processed_dfs <- list()
  
  # 处理每个效应类型
  for (effect in effect_types) {
    # 过滤并处理数据（修复1: 处理列名中的空格）
    df <- mediation_table %>%
      filter(name == effect) %>%
      rename_with(~ paste0(., "_", effect), # 用一个点替换所有空格
                  .cols = -all_of(c(id_vars, "name"))
      ) %>%
      dplyr::select(-name)
    
    processed_dfs[[effect]] <- df
  }
  
  final_df <- Reduce(function(x, y) merge(x, y, by = id_vars), processed_dfs)
  
  final_df <- final_df |>
    mutate(
      # 格式化 ACME
      `ACME (β, 95%CI)` = sprintf(
        "%.3f (%.3f, %.3f)",
        Estimate_ACME,
        `95% CI Lower_ACME`,
        `95% CI Upper_ACME`
      ),
      `p-value_ACME` = round(`p-value_ACME`, 3),
      
      # 格式化 ADE
      `ADE (β, 95%CI)` = sprintf(
        "%.3f (%.3f, %.3f)",
        Estimate_ADE,
        `95% CI Lower_ADE`,
        `95% CI Upper_ADE`
      ),
      `p-value_ADE` = round(`p-value_ADE`, 3),
      
      # 格式化 Total Effect
      `Total Effect (β, 95%CI)` = sprintf(
        "%.3f (%.3f, %.3f)",
        `Estimate_Total Effect`,
        `95% CI Lower_Total Effect`,
        `95% CI Upper_Total Effect`
      ),
      `p-value_Total Effect` = round(`p-value_Total Effect`, 3),
      
      # 格式化 Total Effect
      `Prop. Mediated` = sprintf("%.2f%%", `Estimate_Prop. Mediated` * 100)
    ) |>
    # 选择最终需要的列（按逻辑顺序排列）
    dplyr::select(
      exposure_name, medi_name, outcome_name,
      `ACME (β, 95%CI)`, `p-value_ACME`,
      `ADE (β, 95%CI)`, `p-value_ADE`,
      `Total Effect (β, 95%CI)`, `p-value_Total Effect`,
      starts_with("Prop. Mediated"), # 保留 Prop. Mediated 相关列
      starts_with("p-value_Prop. Mediated")
    )
  
  return(final_df)
}
