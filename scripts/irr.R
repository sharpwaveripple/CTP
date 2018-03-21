library(tidyverse)

df <- read.csv('../data/ctp_irr.csv')
psych::cohen.kappa(df[c('tom_vol', 'dan_vol')])
