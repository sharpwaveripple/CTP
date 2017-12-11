library(pROC)
library(tidyverse)
library(broom)


calc.odds.ratios <- function(logistic.formula, df) {
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
goodCut <- 2

data.raw["poor"] <- as.factor(as.numeric(data.raw[dv] >= poorCut))
data.raw["good"] <- as.factor(as.numeric(data.raw[dv] <= goodCut))
data.raw$Occl <- ordered(data.raw$Occl)
data.raw$pig_all <- as.factor(data.raw$pig_all)

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



logistic.formula <- paste("poor ~", paste(cov.list, collapse = " + "))
results <- calc.odds.ratios(logistic.formula, df.drop)
results.table <- vanilla.table(results, add.rownames = T)
doc <- addParagraph(doc, c(sprintf("Poor")))
doc <- addFlexTable(doc, results.table)
doc <- addParagraph(doc, c("", ""))

logistic.formula <- paste("good ~", paste(cov.list, collapse = " + "))
results <- calc.odds.ratios(logistic.formula, df.drop)
results.table <- vanilla.table(results, add.rownames = T)
doc <- addParagraph(doc, c(sprintf("Good")))
doc <- addFlexTable(doc, results.table)
doc <- addParagraph(doc, c("", ""))

writeDoc(doc, file = output)



logistic.formula <- paste("poor ~", paste(cov.list2, collapse = " + "))
results <- calc.odds.ratios(logistic.formula, data.raw)
results.table <- vanilla.table(results, add.rownames = T)
doc <- addParagraph(doc, c(sprintf("Poor")))
doc <- addFlexTable(doc, results.table)
doc <- addParagraph(doc, c("", ""))

logistic.formula <- paste("good ~", paste(cov.list2, collapse = " + "))
results <- calc.odds.ratios(logistic.formula, data.raw)
results.table <- vanilla.table(results, add.rownames = T)
doc <- addParagraph(doc, c(sprintf("Good")))
doc <- addFlexTable(doc, results.table)
doc <- addParagraph(doc, c("", ""))

writeDoc(doc, file = output)

