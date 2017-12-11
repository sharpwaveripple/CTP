library(psych)

root <- "C:/Users/jonat/Desktop/ctp/"
datafile <- "../data/data_final.csv"
dv <- "mRS_3m"

core.iv <- c("CoreCBF30", "CoreCBV60")
mism.cov <- c("MismCBF30", "MismCBV60")
static.cov <- c("Age", "BslNIH")
occl.cov <- c("occl2", "occl3", "occl4")
pig.cov <- c("pig3", "pig2p", "pig2g")
all.vars <- c(core.iv, mism.cov, static.cov, occl.cov, pig.cov, dv)

data.raw <- read.csv(datafile, header = T, na.strings = "?")
df <- data.raw[all.vars]

df.final <- df[complete.cases(df), ]
#df.final$Side <- as.factor(df.final$Side)
df.final$occl2 <- as.factor(df.final$occl2)
df.final$pig2p <- as.factor(df.final$pig2p)
df.final$pig2g <- as.factor(df.final$pig2g)
df.final$occl3 <- ordered(df.final$occl3)
df.final$occl4 <- ordered(df.final$occl4)
df.final$pig3 <- ordered(df.final$pig3)

poorCut <- 5
goodCut <- 1

df.final["poor"] <- as.factor(as.numeric(df.final[dv] >= poorCut))
df.final["good"] <- as.factor(as.numeric(df.final[dv] <= goodCut))
write.csv(df.final, row.names=FALSE, file="rank.csv")

#df.final <- lapply(df.final, as.numeric)
#df.final <- as.matrix(df.final)
#is.factor(df.final$occl3)
#test <- cor.test(x = df.final['occl3'],y= df.final['pig3'], method="spearman")

#r <- corr.test(df.final, method="kendall")
#print(r, short=F, digits=3)

