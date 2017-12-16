library(tidyverse)
library(broom)

extractP <- function(results) {
  p <- results[seq(2, length(results), 2)]
  p <- unlist(p)
  corr <- p.adjust(p)
  return(corr)
}

calc.odds.ratios <- function(logistic.formula, df) {
  sig.digits <- 3
  format.str <- paste("%.", sig.digits, "f", sep = "")
  logit <- glm(formula(logistic.formula), data = df, family = binomial)
  logit.df <- tidy(logit)
  logit.df %>%
    mutate(or = exp(estimate),
           var.diag = diag(vcov(logit)),
           or.se = sqrt(or^2 * var.diag)) -> logit.df
  results <- round(cbind(logit.df$or, 
                         logit.df$or.se, 
                         exp(confint(logit)), 
                         logit.df$p.value), sig.digits)
  results.format <- matrix(data = sprintf(format.str, results),
                           nrow = dim(results)[1],
                           ncol = dim(results)[2],
                           dimnames = list(rownames(results),
                                           c("Odds ratio",
                                             "Standard error",
                                             "2.5%", "97.5%", "p")))
  results.table <- as.table(results.format)
  return(results.table)
}

bin <- c("Sex", "Diabetes", "HF", "MI/Angina", "Side")
ord <- c("mRS_3m", "occl4", "pig3")
cont <- c("Age", "BslNIH", "tPAdelay",
          "CoreCBF30", "MismCBF30", "CoreCBV60", "MismCBV60")
vars <- c(bin, ord, cont)

read_csv('../data/data_final.csv', na = "?") %>%
  select(vars) %>%
  drop_na(., BslNIH) %>%
  mutate(poor = mRS_3m >= 5,
         moderate = mRS_3m <= 3,
         good = mRS_3m <= 2,
         excellent = mRS_3m <= 1) %>%
  mutate_at(bin, funs(factor(.))) %>%
  mutate_at(ord, funs(ordered(.))) ->
  df
  
  
summarise_at(cont, c("median", "IQR"), na.rm=T)

df %>%
  count(mRS_3m) %>%
  mutate(per=round(n/sum(n), 3)*100) -> clinical

catVars <- c(bin, ord)
outVars <- c('poor', 'moderate', 'good', 'excellent')

catResults <- matrix(, nrow=length(catVars), ncol = 0)
contResults <- matrix(, nrow=length(cont), ncol = 0)
rownames(catResults) <- catVars
rownames(contResults) <- cont

for (outcome in outVars) {
  map(df[catVars], function(x) table(x, df[[outcome]])) %>%
  map(., function(x) chisq.test(x)) %>%
  map_df(., magrittr::extract, c('statistic', 'p.value')) %>%
  round(., 3) -> x2
  catResults <- cbind(catResults, x2)

  map(df[cont], function(x) t.test(x ~ df[[outcome]])) %>%
  map_df(., magrittr::extract, c('statistic', 'p.value')) %>%
  round(., 3) -> t
  contResults <- cbind(contResults, t)
  
  df %>%
    when(
      outcome=
    )
}


map(df[catVars], function(x) table(x, df$mRS_3m)) %>%
map(., function(x) chisq.test(x)) %>%
map_df(., magrittr::extract, c('statistic', 'p.value')) %>%
round(., 3) -> x2

map(df[contVars], function(x) aov(x ~ df$mRS_3m)) %>%
map_df(., magrittr::extract, c('statistic', 'p.value')) %>%
round(., 3) -> t

x2p <- extractP(catResults)
tp <- extractP(contResults)


for (outcome in outVars) {
  logistic.formula <- paste(outcome, "Age + BslNIH + occl4 + pig3 + CoreCBF30 + MismCBF30",
                            sep = " ~ ")
  results <- calc.odds.ratios(logistic.formula, df)
}
