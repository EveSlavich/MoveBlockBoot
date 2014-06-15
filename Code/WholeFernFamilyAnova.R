setwd("~/MovingBlockBootstrap")
load("Data/MovingBlockBootstrap.RData")
#dat$SiteNumber_new[order(dat$SiteNumber_new)]
dat = dat[order(dat$SiteNumber_new),]
Sites$names =Sites$SiteNumber_new
Sites=Sites[order(Sites$names),]
dat = cbind(dat,Sites[, c("x","y","names")])
#head(dat)
which (dat[,"SiteNumber_new"] !=dat[,"names"])
source("~/MovingBlockBootstrap/Functions.R")
names(dat)
ferns=dat[,24:28]
library(mvabund)
presence.ferns=mvabund(ferns)
lookup_table_20k = create_moving_block_lookup_RLAB(dat$x, dat$y,block_L = 20000, scale1 = 5000)
lookup_table_10k = create_moving_block_lookup_RLAB(dat$x, dat$y,block_L = 10000, scale1 = 5000)
lookup_table_5k = create_moving_block_lookup_RLAB(dat$x, dat$y,block_L = 5000, scale1 = 2500)

save(lookup_table_20k, file="DataMod/lookup_table_20k.RData")
save(lookup_table_10k, file="DataMod/lookup_table_10k.RData")
save(lookup_table_5k, file="DataMod/lookup_table_5k.RData")



load ("DataMod/lookup_table_20k.RData")
load ("DataMod/lookup_table_10k.RData")
load ("DataMod/lookup_table_5k.RData")
set.seed(1901)
sample1 = resample_blocks( dat$x, dat$y, NBoot = 1000, lookup_table = lookup_table_20k) 
sample2 = resample_blocks( dat$x, dat$y, NBoot = 1000, lookup_table = lookup_table_10k) 
sample3 = resample_blocks( dat$x, dat$y, NBoot = 1000, lookup_table = lookup_table_5k) 

save(sample1,sample2, sample3, file="Datamod/newsample.RData")
mod.1 = manyglm(presence.ferns~ poly(AnnualMeanTemp,2) + poly(AnnualPrec,2)+ poly(mint5,2)+poly(maxt95,2)+poly(maxt5,2)+poly(mint95,2)+poly(PrecipSeas,2)+poly(PrecipWettMon,2)+poly(PrecipDieMon,2)+poly(MaxTempWarMon,2)+poly(MinTempColMon,2), data=dat, family="binomial")

an1.naive.allferns = anova.manyglm(mod.1, nBoot = 1000)
save(an1.naive.allferns, file="Datamod/an1.naive.allferns.RData")
an1.20k.allferns = anova.manyglm(mod.1, bootID = t(sample1[,1:1000]))
save(an1.20k.allferns, file="Datamod/an1.20k.allferns.RData")
an1.10k.allferns = anova.manyglm(mod.1, bootID = t(sample2[,1:1000]))
save(an1.10k.allferns, file="Datamod/an1.10k.allferns.RData")
an1.5k.allferns=anova.manyglm(mod.1, bootID = t(sample3[,1:1000]))
save(an1.5k.allferns, file="Datamod/an1.5k.allferns.RData")
load("Datamod/an1.5k.allferns.RData")
load("Datamod/an1.10k.allferns.RData")
load("Datamod/an1.20k.allferns.RData")
load("Datamod/an1.naive.allferns.RData")
Results_table = data.frame(an1.naive.allferns$table[,3:4],an1.5k.allferns$table[,4], an1.10k.allferns$table[,4], an1.20k.allferns$table[,4])
names(Results_table)= c("Deviance","Naive", "L = 5km", "L = 10km", "L=20km")
Results_table=Results_table[-1,]
round(Results_table,3)
library(xtable)
table_anova_all_ferns = xtable(Results_table)
caption(table_anova_all_ferns) = ""
print.xtable(table_anova_all_ferns, file = "Tex/table_anova_all_ferns.tex")


#these didn't really show interesting changes
mod.2 = manyglm(presence.ferns~poly(AnnualMeanTemp,2) + poly(AnnualPrec,2)+ poly(mint5,2)+poly(maxt95,2)+poly(maxt5,2)+poly(mint95,2)+poly(PrecipSeas,2)+poly(PrecipWettMon,2)+poly(PrecipDieMon,2)+poly(MaxTempWarMon,2)+poly(MinTempColMon,2)+poly(minh5,2)+poly(maxh95,2) , data=dat, family="binomial")
an1.naive.allferns2 = anova.manyglm(mod.2, nBoot = 100)
#save(an1.naive.allferns, file="Datamod/an1.naive.allferns.RData")
an1.20k.allferns2 = anova.manyglm(mod.2, bootID = t(sample1[,1:100]))

mod.3 = manyglm(presence.ferns~poly(mint5,2)+poly(maxt95,2)+poly(PrecipWettMon,2)+poly(PrecipDieMon,2)+poly(MaxTempWarMon,2)+poly(MinTempColMon,2)+poly(minh5,2)+poly(maxh95,2) , data=dat, family="binomial")
an1.naive.allferns3 = anova.manyglm(mod.3, nBoot = 100)
#save(an1.naive.allferns, file="Datamod/an1.naive.allferns.RData")
an1.20k.allferns3 = anova.manyglm(mod.3, bootID = t(sample1[,1:100]))


