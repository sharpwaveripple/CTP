library(tidyverse)
library(ordinal)

datafile <- "../data/data_final.csv"
dv <- "mRS_3m"
static.cov <- c("Age", "BslNIH",
                "occl4", "pig3")
iv <- c("CoreCBF30", "CoreCBV60")
mism.cov <- c("MismCBF30", "MismCBV60")
all.vars <- c(iv, static.cov, dv, mism.cov)

ord <- c("mRS_3m", "occl4", "pig3")
read_csv(datafile, na="?") %>%
  select(all.vars) %>%
  drop_na -> df

df[ord] <- lapply(df[ord], ordered)

for (i in iv) {
 for (j in mism.cov) {
    model.vars <- c(i, j, static.cov)
    ordinal.formula <- paste("mRS_3m ~", paste(model.vars, collapse = " + "))
    fit <- clm(formula(ordinal.formula), data=df, Hess=T)
    print(exp(coef(fit)[-(1:6)]))
    print(exp(confint(fit)))
    print("=========")
 }
}

