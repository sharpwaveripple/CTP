library(tidyverse)
library(broom)
library(rpart)
library(pROC)

correct_p <- function(pMat, method) {
  pMat %>%
    unlist %>%
    p.adjust(method = method) -> pVec
  dim(pVec) <- dim(pMat)
  dimnames(pVec) <- dimnames(pMat)
  return(pVec)
}

bin <- c("Sex", "Diabetes", "HF", "MI/Angina", "Side")
ord <- c("mRS_3m", "occl4", "pig3")
cont <- c("Age",
          "BslNIH",
          "Core",
          "Penumbra",
          "DT8",
          "24h Infarct Vol (cm3)")
vars <- c(bin, ord, cont)
sig <- 3

read_csv('../data/data_final.csv', na = "?") %>%
  select(vars) %>%
  drop_na(., BslNIH) %>%
  mutate(
    full = mRS_3m,
    dead = mRS_3m == 6,
    poor = mRS_3m >= 5,
    moderate = mRS_3m <= 3,
    good = mRS_3m <= 2,
    excellent = mRS_3m <= 1
  ) %>%
  mutate_at(bin, funs(factor(.))) %>%
  mutate_at(ord, funs(ordered(.))) ->
  df

df %>%
  summarise_at(cont, c("median", "IQR"), na.rm = T)

df %>%
  count(mRS_3m) %>%
  mutate(per = round(n / sum(n), sig) * 100) -> clinical


ord = ord[!ord %in% 'mRS_3m'] # gracelessly pop off mRS

catVars <- c(bin, ord)
outVars <- c('dead', 'poor', 'moderate', 'good', 'excellent')


catResults <- matrix(, nrow = length(catVars), ncol = 0)
contResults <- matrix(, nrow = length(cont), ncol = 0)
rownames(catResults) <- catVars
rownames(contResults) <- cont

for (outcome in outVars) {
  map(df[catVars], function(x)
    table(x, df[[outcome]])) %>%
    map(., function(x)
      chisq.test(x)) %>%
    map_df(., magrittr::extract, c('statistic', 'p.value')) -> x2
  catResults <- cbind(catResults, x2)

  map(df[cont], function(x)
    t.test(x ~ df[[outcome]])) %>%
    map_df(., magrittr::extract, c('statistic', 'p.value')) -> f
  contResults <- cbind(contResults, f)
  print(outcome)
  print(contResults)
}

map(df[catVars], function(x)
  table(x, df$mRS_3m)) %>%
  map(., function(x)
    chisq.test(x)) %>%
  map_df(., magrittr::extract, c('statistic', 'p.value')) -> x2
catResults <- cbind(catResults, x2)

map(df[cont], function(x)
  aov(x ~ df$mRS_3m)) %>%
  map(., function(x)
    tidy(x)) %>%
  map_df(., magrittr::extract, c('statistic', 'p.value')) %>%
  drop_na -> f
contResults <- cbind(contResults, f)

uni <- rbind(catResults, contResults)
uni_stats <- uni[seq(1, length(uni), 2)]
uni_ps <- uni[seq(2, length(uni), 2)]
colnames(uni_ps) <- c(outVars, 'full')
x = correct_p(uni_ps, method = 'holm')

for (outcome in outVars) {
  paste(outcome, "Age + BslNIH + 24h Infarct Vol (cm3)",
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
  print(results)

  #dplyr rearrange columns
}

cartVars <- c("CoreCBF30", "CoreCBV60", "Age", "BslNIH")
outVars <- c(outVars, 'full')

for (outcome in outVars) {
  paste(outcome, cartVars, sep="~") %>%
    map(as.formula) %>%
    map(rpart, data=df) %>%
    map(magrittr::extract('splits')) %>%
    map(as.tibble) %>%
    map(., function(x) select(x, 'index')) %>%
    map(., function(x) slice(x, 1)) %>%
    map(as.numeric) %>%
    unlist ->
    splits

  newnames <- paste(cartVars, "bin", sep="_")

  df %>%
    select(cartVars) < splits ->
    splitVars

  colnames(splitVars) <- newnames

  temp <- cbind(df, splitVars)

  paste(outcome, "CoreCBV60_bin + Age_bin + BslNIH_bin",
        sep="~") %>%
    as.formula %>%
    glm(., data=temp, family=binomial) ->
    logit

    predict(., type="response") ->
    predVals

  temp$prob <- predVals
  roc(as.factor(temp$poor), temp$prob)

}

x <- rpart(poor ~ CoreCBF30, data = df)
x$splits %>% as.tibble() %>% select(index) -> y
slice(y, 1) %>% as.numeric

df %>% mutate(binCore = CoreCBF30 < 66.85) -> test
x2 <- table(test$binCore, test$poor)
# plot his of df$occl4
