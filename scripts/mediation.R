library(psych)
library(lavaan)

datafile <- "../data/data_final.csv"
dv <- "mRS_3m"
all.vars <- c(dv, "CoreCBF30", "pig3", "CoreCBV60")

data.raw <- read.csv(datafile, header = T, na.strings = "?")
df <- data.raw[all.vars]
df.final <- df[complete.cases(df), ]
df.final$mRS_3m <- ordered(df.final$mRS_3m)
df.final$pig3 <- ordered(df.final$pig3)

model <- '# direct effect
          mRS_3m ~ c*pig3
          # mediator
          CoreCBF30 ~ a*pig3
          mRS_3m ~ b*CoreCBF30
          # indirect effect
          ab := a*b
          # total effect
          total := c + (a*b)
'

model2 <- '# direct effect
           mRS_3m ~ c*pig3
           # mediator
           CoreCBV60 ~ a*pig3
           mRS_3m ~ b*CoreCBV60
           # indirect effect
           ab := a*b
           # total effect
           total := c + (a*b)
'

model3 <- 'mRS_3m ~ pig3
'

fit <- sem(model, data=df.final)
fit2 <- sem(model2, data=df.final)
fit3 <- sem(model3, data=df.final)

summary(fit, standardized = T)
#summary(fit2, standardized = T)
#summary(fit3, standardized=T)
