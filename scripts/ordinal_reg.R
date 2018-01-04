library(ordinal)

datafile <- "../data/data_final.csv"
dv <- "mRS_3m"
static.cov <- c("Age", "BslNIH",
                "occl4", "pig3")
iv <- c("CoreCBF30", "CoreCBV60")
mism.cov <- c("MismCBF30", "MismCBV60")
all.vars <- c(iv, static.cov, dv, mism.cov)



data.raw <- read.csv(datafile, header = T, na.strings = "?")
df <- data.raw[all.vars]
df.final <- df[complete.cases(df), ]
df.final$mRS_3m <- ordered(df.final$mRS_3m)
df.final$occl4 <- ordered(df.final$occl4)
df.final$pig3 <- ordered(df.final$pig3)

for (i in iv) {
 for (j in mism.cov) {
    model.vars <- c(i, j, static.cov)
    ordinal.formula <- paste("mRS_3m ~", paste(model.vars, collapse = " + "))
    fit <- clm(formula(ordinal.formula), data=df.final, Hess=T)
    print(summary(fit))
 }
}
