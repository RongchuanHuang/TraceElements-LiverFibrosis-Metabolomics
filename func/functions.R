
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
    mutate(`Percent change (95%CI)` = glue("{Estimate} ({LOW_CI}, {UP_CI})"))
}