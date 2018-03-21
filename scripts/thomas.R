library(tidyverse)
library(broom)
library(psych)

#IV = 24h infarct volume
#ICH = 24h Haemorrhage
#ICHV = ICH volume
#ICH_exclude = exclusion because big ICH

format_correlation <- function(corr_mat) {
  precision <- 3
  format_str <- paste("%.", precision, "f", sep="")
  var_names <- list(colnames(corr_mat$r))
  r <- sprintf(eval(format_str), corr_mat$r)
  p <- sprintf(eval(format_str), corr_mat$p)
  p <- gsub("0.000", "< 0.001", p)
  tab <- matrix(paste(r, " (", p, ")", sep=""),
                nrow=nrow(corr_mat$r), dimnames=c(var_names, var_names))
  return(tab)
}

bin <- c("Sex", "Side", "ICH", "Randomise", "Recanal", "Recanalandrandom")
ord <- c("mRS_3m", "occl4", "pig3")
outVars <- c("dead", "poor", "moderate", "good", "excellent")
cont <- c("Age",
          "BslNIH",
          "Core_no_zero",
          "Penumbra",
          "Mismatch",
          "IV")

vars <- c(bin, ord, cont, outVars)

read_csv("../data/data_final.csv", na = c("", "?")) %>%
  mutate(
    full = mRS_3m,
    dead = mRS_3m == 6,
    poor = mRS_3m >= 5,
    moderate = mRS_3m <= 3,
    good = mRS_3m <= 2,
    excellent = mRS_3m <= 1,
    Mismatch = Total/Core_no_zero,
    Recanal = IV < Total / 2,
    Randomise = Randomise == 1,
    ICH = recode_factor(ICH,
                        "Petechial haemorrhage" = 0,
                        "Small SAH" = 0,
                        "Yes" = 1,
                        "No" = 0),
    ICH = ICH == 1,
    Recanalandrandom = (Recanal == T) & (Randomise == 1),
    Sex = recode_factor(Sex, "M" = 0, "F" = 1),
    Side = recode_factor(Side, "R" = 0, "L" = 1)) %>%
  mutate_at(ord, funs(ordered(.))) %>%
  select(vars) %>%
  drop_na ->
  df

ord = ord[!ord %in% "mRS_3m"] # gracelessly pop off mRS
catVars <- c(bin, ord)

print(format_correlation(corr.test(df[cont])))

catResults <- matrix(, nrow = length(catVars), ncol = 0)
contResults <- matrix(, nrow = length(cont), ncol = 0)
rownames(catResults) <- catVars
rownames(contResults) <- cont

for (outcome in outVars) {
  map(df[catVars], function(x)
    table(x, df[[outcome]])) %>%
    map(., function(x)
      chisq.test(x)) %>%
    map_df(., magrittr::extract, c("statistic", "p.value")) -> x2
  catResults <- cbind(catResults, x2)

  map(df[cont], function(x)
    t.test(x ~ df[[outcome]])) %>%
    map_df(., magrittr::extract, c("statistic", "p.value")) -> f
  contResults <- cbind(contResults, f)
}

write.table(round(catResults, 3), file = "cat_univariate.txt")
write.table(round(contResults, 3), file = "cont_univariate.txt")

for (outcome in outVars) {
  paste(outcome,
        "Age + BslNIH + Core_no_zero + Mismatch + occl4 + pig3 + Recanal",
        sep = " ~ ") %>%
    glm(., data = df, family = binomial) ->
    logit

  logit %>%
    confint %>%
    exp ->
    ci

  logit %>%
    tidy() %>%
    mutate(or = exp(estimate)) %>%
    select(or, p.value) %>%
    cbind(., ci) ->
    results

  print(outcome)
  print(round(results, 3))
}

fit <- lm(IV ~ Age + BslNIH + Core_no_zero + Mismatch + occl4 + pig3, data=df)
print(summary(fit))
