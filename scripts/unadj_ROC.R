library(pROC)
library(tidyverse)
library(broom)
library(MASS)


datafile <- "../data/ctp_data.csv"
dv <- "mRS_3m"
cov.list <- c("CoreCBV60","MismCBV60",
              "Age", "BslNIH", "Side",
              "Occl", "corticalscore")


data.raw <- read.csv(datafile, header = T, na.strings = "?")
df <- data.raw[cov.list]
df.drop <- df[complete.cases(df), ]
varNames <- as.vector(outer(iv.prefix, iv.thresh, paste, sep=""))
sig.digits = 2
n.bootstraps = 10000 
format.str = paste("%.", sig.digits, "f", sep = "")

poorCut <- 5
goodCut <- 1

data.raw["poor"] <- as.factor(as.numeric(data.raw[dv] >= poorCut))
data.raw["good"] <- as.factor(as.numeric(data.raw[dv] <= goodCut))
occl.cov <- c("occl2", "occl3", "occl4")
pig.cov <- c("pig3", "pig2p", "pig2g")

for (thresh in iv.thresh) {
  vars <- as.vector(outer(iv.prefix, thresh, paste, sep=""))
  rocForm1 <- paste("poor ~", vars[1])
  rocForm2 <- paste("poor ~", vars[2])
  rocForm3 <- paste("good ~", vars[1])
  rocForm4 <- paste("good ~", vars[2])
  roc1 <- roc(formula(rocForm1), data=data.raw)
  roc2 <- roc(formula(rocForm2), data=data.raw)
  roc3 <- roc(formula(rocForm3), data=data.raw)
  roc4 <- roc(formula(rocForm4), data=data.raw)
}

