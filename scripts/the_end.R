library(pROC)
library(tidyverse)
library(broom)
library(ReporteRs)


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


datafile <- "../data/data_final.csv"
dv <- "mRS_3m"

core.iv <- c("CoreCBF30", "CoreCBV60")
mism.cov <- c("MismCBF30", "MismCBV60")
static.cov <- c("Age", "BslNIH", "Side")
occl.cov <- c("occl2", "occl3", "occl4")
pig.cov <- c("pig3", "pig2p", "pig2g")
all.vars <- c(core.iv, mism.cov, static.cov, occl.cov, pig.cov, dv)
poorCut <- 5
goodCut <- 2
sig.digits = 3
format.str = paste("%.", sig.digits, "f", sep = "")


data.raw <- read.csv(datafile, header = T, na.strings = "?")
df <- data.raw[all.vars]
df.final <- df[complete.cases(df), ]
df.final$Side <- as.factor(df.final$Side)
df.final$occl2 <- as.factor(df.final$occl2)
df.final$pig2p <- as.factor(df.final$pig2p)
df.final$pig2g <- as.factor(df.final$pig2g)
df.final$occl3 <- ordered(df.final$occl3)
df.final$occl4 <- ordered(df.final$occl4)
df.final$pig3 <- ordered(df.final$pig3)

df.final["poor"] <- as.factor(as.numeric(df.final[dv] >= poorCut))
df.final["good"] <- as.factor(as.numeric(df.final[dv] <= goodCut))

doc <- docx()

for (i in core.iv) {
  for (j in mism.cov) {
    for (k in occl.cov) {
      for (l in pig.cov) {
        model.vars <- c(i, j, k, l, static.cov)

        logistic.formula <- paste("poor ~", paste(model.vars, collapse = " + "))
        results <- calc.odds.ratios(logistic.formula, df.final)
        results.table <- vanilla.table(results, add.rownames = T)
        doc <- addParagraph(doc, c(logistic.formula, ""))
        doc <- addFlexTable(doc, results.table)
        doc <- addParagraph(doc, c("", ""))

        logistic.formula <- paste("good ~", paste(model.vars, collapse = " + "))
        results <- calc.odds.ratios(logistic.formula, df.final)
        results.table <- vanilla.table(results, add.rownames = T)
        doc <- addParagraph(doc, c(logistic.formula, ""))
        doc <- addFlexTable(doc, results.table)
        doc <- addParagraph(doc, c("", ""))
      }
    }
  }
}
writeDoc(doc, file="final_results.docx")

doc <- docx()
for (i in core.iv) {
    for (k in occl.cov) {
      for (l in pig.cov) {
        model.vars <- c(i, k, l, static.cov)
        
        logistic.formula <- paste("poor ~", paste(model.vars, collapse = " + "))
        results <- calc.odds.ratios(logistic.formula, df.final)
        results.table <- vanilla.table(results, add.rownames = T)
        doc <- addParagraph(doc, c(logistic.formula, ""))
        doc <- addFlexTable(doc, results.table)
        doc <- addParagraph(doc, c("", ""))
        
        logistic.formula <- paste("good ~", paste(model.vars, collapse = " + "))
        results <- calc.odds.ratios(logistic.formula, df.final)
        results.table <- vanilla.table(results, add.rownames = T)
        doc <- addParagraph(doc, c(logistic.formula, ""))
        doc <- addFlexTable(doc, results.table)
        doc <- addParagraph(doc, c("", ""))
      }
    }
  }

writeDoc(doc, file="final_results_nomism.docx")


doc <- docx()
for (i in core.iv) {
      model.vars <- c(i, static.cov)
      logistic.formula <- paste("poor ~", paste(model.vars, collapse = " + "))
      results <- calc.odds.ratios(logistic.formula, df.final)
      results.table <- vanilla.table(results, add.rownames = T)
      doc <- addParagraph(doc, c(logistic.formula, ""))
      doc <- addFlexTable(doc, results.table)
      doc <- addParagraph(doc, c("", ""))
      
      logistic.formula <- paste("good ~", paste(model.vars, collapse = " + "))
      results <- calc.odds.ratios(logistic.formula, df.final)
      results.table <- vanilla.table(results, add.rownames = T)
      doc <- addParagraph(doc, c(logistic.formula, ""))
      doc <- addFlexTable(doc, results.table)
      doc <- addParagraph(doc, c("", ""))
}
writeDoc(doc, file="final_results_only_sig.docx")
