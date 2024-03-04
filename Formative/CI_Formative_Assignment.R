## Clear the R environment
rm(list = ls())

## Load packages
library(foreign)
library(ggplot2)
library(scales)
library(dplyr)
library(plyr)
library(car)
library(stargazer)
library(lmtest)
library(sandwich)
library(plyr)
library(reshape)
library(Amelia)
library(AER)

## Cluster SE function
cl = function(dat, fm, cluster) {
  
  require(sandwich, quietly = TRUE)
  require(lmtest, quietly = TRUE)
  M = length(unique(cluster))
  N = length(cluster)
  K = fm$rank
  dfc = (M / (M - 1)) * ((N - 1) / (N - K))
  uj  = apply(estfun(fm), 2, function(x) tapply(x, cluster, sum));
  vcovCL = dfc * sandwich(fm, meat = crossprod(uj) / N)
  coeftest(fm, vcovCL)
  
}

## Set working directory
setwd("/Users/alasd/Documents/Git_Projects/ci-assignment/dataverse_files")

## Read csv data
police.data.save = read.csv("HuCon_ISQ_Data.csv")

police.data = read.csv("HuCon_ISQ_Data.csv")

## Delete DAMAN & DIU 2001
police.data = police.data[-which(is.na(police.data$death_not_remanded)), ]

##############
###Table A1###
##############
stargazer(police.data, median = T)

###############
###Figure A1###
###############

## min max for death
min(police.data$death_remanded, na.rm = T)
max(police.data$death_remanded, na.rm = T)

min(police.data$death_not_remanded, na.rm = T)
max(police.data$death_not_remanded, na.rm = T)

sum(police.data$death_remanded == 0, na.rm = T) 
sum(police.data$death_not_remanded == 0, na.rm = T) 

death.state = ddply(police.data, .(state_ut), summarise, sum = sum(death_remanded))
death.state$not_remanded = ddply(police.data, .(state_ut), summarise, sum = sum(death_not_remanded))

death.year = ddply(police.data, .(year), summarise, sum = sum(death_remanded, na.rm = T))
death.year.not = ddply(police.data, .(year), summarise, sum.not = sum(death_not_remanded, na.rm = T))

merge = merge(death.year, death.year.not, by = "year")

merge.long = melt(merge, id = "year")
names(merge.long)[2] = "Variable"
f.a1 = ggplot(merge.long, aes(year, value, colour = Variable)) + geom_line() + scale_x_continuous(breaks = c(2002, 2004, 2006, 2008, 2010, 2012, 2014, 2016)) + scale_color_manual(labels = c("Death remanded", "Death not remanded"), values = c("#F8766D", "#00BFC4")) + ylab("Count") + xlab("Year")

f.a1

ggsave("death_time.pdf", f.a1, width = 6, height = 4)

##############
###Table A3###
##############
police.data.2006.all = subset(police.data, year <= 2006)

death.state.2006.all = ddply(police.data.2006.all, .(state_ut), summarise, remanded.2006.all = sum(death_remanded, na.rm = T))
death.state.not.2006.all = ddply(police.data.2006.all, .(state_ut), summarise, notremanded.2006.all = sum(death_not_remanded, na.rm = T))

death.state.2006.all$notremanded.2006.all = death.state.not.2006.all$notremanded.2006.all

death.state.2006.all

##############
###Table A4###
##############
police.data.2006 = subset(police.data, year == 2006)

death.state.2006 = ddply(police.data.2006, .(state_ut), summarise, remanded.2006 = sum(death_remanded, na.rm = T))
death.state.not.2006 = ddply(police.data.2006, .(state_ut), summarise, notremanded.2006 = sum(death_not_remanded, na.rm = T))

death.state.2006$notremanded.2006 = death.state.not.2006$notremanded.2006

death.state.2006

################
## Imputation ##
################
## Because Multiple Imputation is a random process, results are slightly different every time
## Load data
police.imp = police.data.save[,  c("state_ut", "year", "death_remanded", "death_not_remanded", "state_pca", "district_pca", "type", "sc_order1", "committee1",
                                 "gdp", "religion2", "head_trans")]

## AmeliaView()

## Multiple imputation with settings below
bds.3 = c(3, 0, 100)
bds.4 = c(4, 0, 100)
bds.12 = c(12, 0, 50)
bds = rbind(bds.3, bds.4, bds.12)

a.out = amelia(police.imp, m = 5, idvars = "type", 
                ts = "year", cs = "state_ut", priors = NULL, lags = "gdp", 
                empri = 0, intercs = TRUE, leads = "gdp", splinetime = 0, 
                logs = c("gdp", "head_trans"), sqrts = NULL, 
                lgstc = NULL, ords = NULL, noms = c("state_pca", "district_pca", 
                                                    "sc_order1", "committee1", "religion2"), bounds = bds, max.resample = 1000, 
                tolerance = 1e-04)

#write.amelia(obj = a.out, file.stem = "outdata")

#############
###Table 1###
#############

police.data.t1 = police.data[, c("death_not_remanded", "death_remanded", "state_ut", "year", "state_pca", "t")]

police.data.t1 = na.omit(police.data.t1)

## Lagged state_pca
police.data.t1 = ddply(police.data.t1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))

## fill NA with 0
police.data.t1$l.state_pca = ifelse(is.na(police.data.t1$l.state_pca), 0, police.data.t1$l.state_pca)

## Table 1 model
model.poisson.t1 = glm(death_not_remanded ~ 1 + l.state_pca + state_ut + as.factor(year), data = police.data.t1, family = "poisson")

model.poisson.t1.cl = cl(police.data.t1, model.poisson.t1, police.data.t1$state_ut)

stargazer(model.poisson.t1.cl)

## predict death count if all PCAs are inplemented on time
police.imp.p = police.data.t1
police.imp.p$l.state_pca = ifelse(police.imp.p$year >= 2008, 1, 0)
Y = predict(model.poisson.t1, police.imp.p, type = "response")

sum(Y)
sum(police.data.t1$death_not_remanded) - sum(Y)

## predict death count if no PCA is inplemented.
police.imp.p = police.data.t1
police.imp.p$l.state_pca = 0
Y.2 = predict(model.poisson.t1, police.imp.p, type = "response")

sum(Y.2)
sum(Y.2) - sum(police.data.t1$death_not_remanded)

##############
###Figure 1###
##############

## Lagged state_pca
police.data.f1 = ddply(police.data.t1, .(state_ut), transform, tm1 = lead(t))
police.data.f1 = ddply(police.data.f1, .(state_ut), transform, tm2 = lead(tm1))
police.data.f1 = ddply(police.data.f1, .(state_ut), transform, tm3 = lead(tm2))
police.data.f1 = ddply(police.data.f1, .(state_ut), transform, tm4 = lead(tm3))

police.data.f1 = ddply(police.data.f1, .(state_ut), transform, tp1 = lag(t))
police.data.f1 = ddply(police.data.f1, .(state_ut), transform, tp2 = lag(tp1))
police.data.f1 = ddply(police.data.f1, .(state_ut), transform, tp3 = lag(tp2))
police.data.f1 = ddply(police.data.f1, .(state_ut), transform, tp4 = lag(tp3))

police.data.f1[is.na(police.data.f1)] = 0

## Poisson Placebo Test
model.poisson.plb = glm(death_not_remanded ~ 1 + tm3 + tm2 + tm1 + t + tp1 + tp2 + tp3 + state_ut + as.factor(year), data = police.data.f1, family = "poisson")

model.poisson.plb.cl = cl(police.data.f1, model.poisson.plb, police.data.f1$state_ut)

stargazer(model.poisson.plb.cl)

## Overdispersion test
dispersiontest(model.poisson.plb, trafo = 1)

## Graph Placebo Test Figure 3
## Save Ts Poisson result
graph.f1 = as.data.frame(model.poisson.plb.cl[2:8, ])
graph.f1$time = c(-3, -2, -1, 0, 1, 2, 3)

## Calculate CIs
graph.f1$ci.l = graph.f1[, 1] - qnorm(0.975) * graph.f1[, 2]
graph.f1$ci.u = graph.f1[, 1] + qnorm(0.975) * graph.f1[, 2]

graph.f1$ci.l.90 = graph.f1[, 1] - qnorm(0.95) * graph.f1[, 2]
graph.f1$ci.u.90 = graph.f1[, 1] + qnorm(0.95) * graph.f1[, 2]

## Plot
p.placebo = ggplot(graph.f1, aes(time, Estimate)) +
  #geom_ribbon(aes(ymin=ci.l,ymax=ci.u),alpha=0.3) +
  geom_errorbar(aes(ymin = ci.l, ymax = ci.u), width = 0.3, color = "#999999") +
  #geom_errorbar(aes(ymin=ci.l.90,ymax=ci.u.90),width=0.1, color = "#999999") +
  geom_pointrange(aes(ymin = ci.l.90, ymax = ci.u.90), size = 1.5, shape = 46, color = "#999999") +
  geom_point(size = 2) +
  geom_line() +
  ylim(-1.1, 1.1) +
  xlab("Years from PCA Creation") +
  ylab("Coefficient of PCA Creation") +
  #geom_line(aes(y=ci.l)) +
  #geom_line(aes(y=ci.u)) +
  #geom_line(aes(y=ci.l.90), linetype = "dashed") +
  #geom_line(aes(y=ci.u.90), linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3))

p.placebo

ggsave("p_placebo_good_2016.pdf", plot = p.placebo, height = 4.5, width = 4.5)

#############
###Table 2###
#############

## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + l.state_pca + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:4, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:4, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t2 = matrix(NA, nrow = 3, ncol = 3)
for (i in 1:3) {
  
  result.t2[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

## T2 results 
## Row: State PCA, State Capacity, and State Desire
## Column: Effect, SE, P value
result.t2

#############
###Table 3###
#############
## Add SHRC to police data
police.data.t1$SHRC = police.data$SHRC
  
## Lagged SHRC
police.data.t1 = ddply(police.data.t1, .(state_ut), transform, l.SHRC = c(NA, SHRC[-length(SHRC)]))

## Fill NA with 0
police.data.t1$l.SHRC = ifelse(is.na(police.data.t1$l.SHRC), 0, police.data.t1$l.SHRC)

## Correlation check
cor.test(police.data.t1$state_pca, police.data.t1$SHRC)

## Model with SHRC
model.poisson.SHRC = glm(death_not_remanded ~ 1 + l.SHRC + state_ut + as.factor(year), data = police.data.t1, family = "poisson")

model.poisson.SHRC.cl = cl(police.data.t1, model.poisson.SHRC, police.data.t1$state_ut)

stargazer(model.poisson.SHRC.cl)

## Model SHRC with controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Add SHRC
  police.imp.1.l$l.SHRC = police.data.t1$l.SHRC
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + l.SHRC + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:4, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:4, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 3, ncol = 3)

for (i in 1:3) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

## T3 (2) results 
## Row: State PCA, State Capacity, and State Desire
## Column: Effect, SE, P value
result.t3

##############
###Table A5###
##############
police.imp.d =  police.data.save[,  c("state_ut", "year", "death_remanded", "death_not_remanded", "state_pca", "district_pca", "type", "sc_order1", "committee1",
                                       "gdp", "religion2", "head_trans")]

stargazer(police.imp.d, median = T)

##############
###Table A6###
##############
##  OLS Placebo Test
police.data.f3 = police.data.f1
model.ols.plb = lm(death_not_remanded ~ 1 + tm3 + tm2 + tm1 + t + tp1 + tp2 + tp3 + state_ut + as.factor(year), data = police.data.f3)
model.ols.plb.cl = cl(police.data.f3, model.ols.plb, police.data.f3$state_ut)

## OLS Placebo Test with controls
## Loop models for 5 imputation datasets
i = 1
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add t to outdata
  police.imp.1.l$t = police.data.t1$t
  
  ## lags and leads
  police.data.f3 = ddply(police.imp.1.l, .(state_ut), transform, tm1 = lead(t))
  police.data.f3 = ddply(police.data.f3, .(state_ut), transform, tm2 = lead(tm1))
  police.data.f3 = ddply(police.data.f3, .(state_ut), transform, tm3 = lead(tm2))
  police.data.f3 = ddply(police.data.f3, .(state_ut), transform, tm4 = lead(tm3))
  
  police.data.f3 = ddply(police.data.f3, .(state_ut), transform, tp1 = lag(t))
  police.data.f3 = ddply(police.data.f3, .(state_ut), transform, tp2 = lag(tp1))
  police.data.f3 = ddply(police.data.f3, .(state_ut), transform, tp3 = lag(tp2))
  police.data.f3 = ddply(police.data.f3, .(state_ut), transform, tp4 = lag(tp3))
  
  police.data.f3[is.na(police.data.f3)] = 0
  
  ## Poisson Placebo Test
  imp.1.p = lm(death_not_remanded ~ 1 + tm3 + tm2 + tm1 + t + tp1 + tp2 + tp3 + gdp + 
                   head_trans+ state_ut + as.factor(year), data = police.data.f3)
  
  result.p.1 = cl(police.data.f3, imp.1.p, police.data.f3$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:10, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:10, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t2 = matrix(NA, nrow = 9, ncol = 3)
for (i in 1:9) {
  
  result.t2[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

## Replace results to model result
model.ols.plb.cl.c = result.p.1

model.ols.plb.cl.c[2:10, 1] = result.t2[, 1]
model.ols.plb.cl.c[2:10, 2] = result.t2[, 2]
model.ols.plb.cl.c[2:10, 4] = result.t2[, 3]

model.ols.plb.cl.c

## Poisson Placebo Test with controls
## Loop models for 5 imputation datasets
i = 1
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add t to outdata
  police.imp.1.l$t = police.data.t1$t
  
  ## lags and leads
  police.data.f3 = ddply(police.imp.1.l, .(state_ut), transform, tm1 = lead(t))
  police.data.f3 = ddply(police.data.f3, .(state_ut), transform, tm2 = lead(tm1))
  police.data.f3 = ddply(police.data.f3, .(state_ut), transform, tm3 = lead(tm2))
  police.data.f3 = ddply(police.data.f3, .(state_ut), transform, tm4 = lead(tm3))
  
  police.data.f3 = ddply(police.data.f3, .(state_ut), transform, tp1 = lag(t))
  police.data.f3 = ddply(police.data.f3, .(state_ut), transform, tp2 = lag(tp1))
  police.data.f3 = ddply(police.data.f3, .(state_ut), transform, tp3 = lag(tp2))
  police.data.f3 = ddply(police.data.f3, .(state_ut), transform, tp4 = lag(tp3))
  
  police.data.f3[is.na(police.data.f3)] = 0
  
  ## Poisson Placebo Test
  imp.1.p = glm(death_not_remanded ~ 1 + tm3 + tm2 + tm1 + t + tp1 + tp2 + tp3 + gdp + 
                   head_trans+ state_ut + as.factor(year), data = police.data.f3, family = "poisson")
  
  result.p.1 = cl(police.data.f3, imp.1.p, police.data.f3$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:10, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:10, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t2 = matrix(NA, nrow = 9, ncol = 3)
for (i in 1:9) {
  
  result.t2[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

## Replace results to model result
model.poisson.plb.cl.c = result.p.1

model.poisson.plb.cl.c[2:10, 1] = result.t2[, 1]
model.poisson.plb.cl.c[2:10, 2] = result.t2[, 2]
model.poisson.plb.cl.c[2:10, 4] = result.t2[, 3]

model.poisson.plb.cl.c

## Make table
# ERROR 
#stargazer(model.ols.plb.cl, model.ols.plb.cl.c, model.poisson.plb.cl, model.poisson.plb.cl.c)

###############
###Figure A3###
###############
## Save Ts Poisson result
graph.a1 = as.data.frame(model.poisson.plb.cl.c[2:8, ])
graph.a1$time = c(-3, -2, -1, 0, 1, 2, 3)

## Calculate CIs
graph.a1$ci.l = graph.a1[, 1] - qnorm(0.975) * graph.a1[, 2]
graph.a1$ci.u = graph.a1[, 1] + qnorm(0.975) * graph.a1[, 2]

graph.a1$ci.l.90 = graph.a1[, 1] - qnorm(0.95) * graph.a1[, 2]
graph.a1$ci.u.90 = graph.a1[, 1] + qnorm(0.95) * graph.a1[, 2]

## Plot
p.placebo.a3 = ggplot(graph.a1, aes(time, Estimate)) +
  #geom_ribbon(aes(ymin=ci.l,ymax=ci.u),alpha=0.3) +
  geom_errorbar(aes(ymin = ci.l, ymax = ci.u), width = 0.3, color = "#999999") +
  #geom_errorbar(aes(ymin=ci.l.90,ymax=ci.u.90),width=0.1, color = "#999999") +
  geom_pointrange(aes(ymin = ci.l.90, ymax = ci.u.90), size = 1.5, shape = 46, color = "#999999") +
  geom_point(size = 2) +
  geom_line() +
  ylim(-1.2, 1.2) +
  xlab("Years from PCA Creation") +
  ylab("Coefficient of PCA Creation") +
  #geom_line(aes(y=ci.l)) +
  #geom_line(aes(y=ci.u)) +
  #geom_line(aes(y=ci.l.90), linetype = "dashed") +
  #geom_line(aes(y=ci.u.90), linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3))

p.placebo.a3

ggsave("p_placebo_controls_2016.pdf", plot = p.placebo.a3, height = 4.8, width = 4.5)

##############
###Table A7###
##############
## Load GTD data
police.data.ta5 = police.data.t1

police.data.ta5$l.event = police.data$l.event

#police.data.save = merge(police.data, gtd.sum.l, by = c("state_ut", "year"), all.x = T)

#police.data.save  = subset(police.data.save, select=-c(iyear, provstate))

#write.csv(police.data.save, "final1.csv")

## fill NA with 0
#police.data.save$l.event = ifelse(is.na(police.data.save$l.event), 0, police.data.save$l.event)

##Correlation check
cor.test(police.data.ta5$l.event, police.data.ta5$l.state_pca)

## OLS with GTD
model.ols.GTD = lm(death_not_remanded ~ 1 + l.state_pca + l.event + state_ut + as.factor(year), data = police.data.ta5)

model.ols.GTD.cl = cl(police.data.ta5, model.ols.GTD, police.data.ta5$state_ut)

## Poisson with GTD
model.poisson.GTD = glm(death_not_remanded ~ 1 + l.state_pca + l.event + state_ut + as.factor(year), data = police.data.ta5, family = "poisson")

model.p.GTD.cl = cl(police.data.ta5, model.poisson.GTD, police.data.ta5$state_ut)

## Poisson with GTD and Controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add GTD l.event to outdata
  police.imp.1.l$l.event = police.data.ta5$l.event
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded ~ 1 + l.state_pca + l.event + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t2 = matrix(NA, nrow = 4, ncol = 3)
for (i in 1:4) {
  
  result.t2[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t2

## Replace results to model result
model.ols.GTD.cl.c = result.p.1

model.ols.GTD.cl.c[2:5, 1] = result.t2[, 1]
model.ols.GTD.cl.c[2:5, 2] = result.t2[, 2]
model.ols.GTD.cl.c[2:5, 4] = result.t2[, 3]

model.ols.GTD.cl.c

## Poisson with GTD and Controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add GTD l.event to outdata
  police.imp.1.l$l.event = police.data.ta5$l.event
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + l.state_pca + l.event + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t2 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t2[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t2

## Replace results to model result
model.p.GTD.cl.c = result.p.1

model.p.GTD.cl.c[2:5, 1] = result.t2[, 1]
model.p.GTD.cl.c[2:5, 2] = result.t2[, 2]
model.p.GTD.cl.c[2:5, 4] = result.t2[, 3]

stargazer(model.ols.GTD.cl, model.ols.GTD.cl.c, model.p.GTD.cl, model.p.GTD.cl.c)

##############
###Table A8###
##############
## Add religion to police data
police.data.t1$religion2 = police.data$religion2

## OLS Model with religion
model.ols.religion = lm(death_not_remanded ~ 1 + l.state_pca + religion2 + state_ut + as.factor(year), data = police.data.t1)

model.ols.religion.cl = cl(police.data.t1, model.ols.religion, police.data.t1$state_ut)

## Poisson Model with religion
model.poisson.religion = glm(death_not_remanded ~ 1 + l.state_pca + religion2 + state_ut + as.factor(year), data = police.data.t1, family = "poisson")

model.poisson.religion.cl = cl(police.data.t1, model.poisson.religion, police.data.t1$state_ut)

## OLS Model with religion
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add religion
  police.imp.1.l$religion2 = police.data.t1$religion2
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded ~ 1 + l.state_pca + religion2 + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.ols.religion.cl.c = result.p.1

model.ols.religion.cl.c[2:5, 1] = result.t3[, 1]
model.ols.religion.cl.c[2:5, 2] = result.t3[, 2]
model.ols.religion.cl.c[2:5, 4] = result.t3[, 3]

## Poisson Model with religion and controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add religion
  police.imp.1.l$religion2 = police.data.t1$religion2
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + l.state_pca + religion2 + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.poisson.religion.cl.c = result.p.1

model.poisson.religion.cl.c[2:5, 1] = result.t3[, 1]
model.poisson.religion.cl.c[2:5, 2] = result.t3[, 2]
model.poisson.religion.cl.c[2:5, 4] = result.t3[, 3]

# ERROR 
#stargazer(model.ols.religion.cl, model.ols.religion.cl.c, model.poisson.religion.cl, model.poisson.religion.cl.c)


##############
###Table A9###
##############
## OLS
model.ols = lm(death_not_remanded ~ 1 + l.state_pca + state_ut + 
                      as.factor(year), data = police.data.t1)

model.ols.cl = cl(police.data.t1, model.ols, police.data.t1$state_ut)


## OLS with logged DV
police.data.t1$death_not_remanded_ln = log(police.data.t1$death_not_remanded + 1)

model.ols.log = lm(death_not_remanded_ln ~ 1 + l.state_pca + state_ut + 
                as.factor(year), data = police.data.t1)

model.ols.log.cl = cl(police.data.t1, model.ols.log, police.data.t1$state_ut)

## OLS with Controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded ~ 1 + l.state_pca + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:4, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:4, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t2 = matrix(NA, nrow = 3, ncol = 3)

for (i in 1:3) {
  
  result.t2[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t2

## Replace results to model result
model.ols.cl.c = result.p.1

model.ols.cl.c[2:4, 1] = result.t2[, 1]
model.ols.cl.c[2:4, 2] = result.t2[, 2]
model.ols.cl.c[2:4, 4] = result.t2[, 3]

## OLS with logged DV and Controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Log DV
  police.imp.1$death_not_remanded_ln = log(police.imp.1$death_not_remanded + 1)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded_ln ~ 1 + l.state_pca + gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:4, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:4, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t2 = matrix(NA, nrow = 3, ncol = 3)

for (i in 1:3) {
  
  result.t2[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t2

## Replace results to model result
model.ols.log.cl.c = result.p.1

model.ols.log.cl.c[2:4, 1] = result.t2[, 1]
model.ols.log.cl.c[2:4, 2] = result.t2[, 2]
model.ols.log.cl.c[2:4, 4] = result.t2[, 3]

stargazer(model.ols.cl, model.ols.cl.c, model.ols.log.cl, model.ols.log.cl.c)


###############
###Table A10###
###############
## OLS no lag
model.ols.nl = lm(death_not_remanded ~ 1 + state_pca + state_ut + 
                  as.factor(year), data = police.data.t1)

model.ols.nl.cl = cl(police.data.t1, model.ols.nl, police.data.t1$state_ut)

## Poisson no lag
model.p.nl = glm(death_not_remanded ~ 1 + state_pca + state_ut + 
                     as.factor(year), data = police.data.t1, family = "poisson")

model.p.nl.cl = cl(police.data.t1, model.p.nl, police.data.t1$state_ut)

## OLS no lag with Controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded ~ 1 + state_pca + gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:4, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:4, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t2 = matrix(NA, nrow = 3, ncol = 3)

for (i in 1:3) {
  
  result.t2[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t2

## Replace results to model result
model.ols.nl.cl.c = result.p.1

model.ols.cl.c[2:4, 1] = result.t2[, 1]
model.ols.cl.c[2:4, 2] = result.t2[, 2]
model.ols.cl.c[2:4, 4] = result.t2[, 3]

## Poisson no lag with Controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + state_pca + gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:4, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:4, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t2 = matrix(NA, nrow = 3, ncol = 3)

for (i in 1:3) {
  
  result.t2[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t2

## Replace results to model result
model.p.nl.cl.c = result.p.1

model.ols.log.cl.c[2:4, 1] = result.t2[, 1]
model.ols.log.cl.c[2:4, 2] = result.t2[, 2]
model.ols.log.cl.c[2:4, 4] = result.t2[, 3]

stargazer(model.ols.nl.cl, model.ols.nl.cl.c, model.p.nl.cl, model.p.nl.cl.c)


###############
###Table A11###
###############
## Balanced Pabel
police.data.b = subset(police.data.t1, police.data.t1$state_ut !=  "TELANGANA")
police.data.b = subset(police.data.b, police.data.b$state_ut !=  "Z DAMAN & DIU")

police.data.b$state_ut = as.factor(as.character(police.data.b$state_ut))
levels(police.data.b$state_ut)
length(police.data.b$death_not_remanded)

model.poisson.b = glm(death_not_remanded ~ 1 + l.state_pca + state_ut + as.factor(year), data = police.data.b, family = "poisson")
model.poisson.b.cl = cl(police.data.b, model.poisson.b, police.data.b$state_ut)

##Quasi-poisson
model.qp = glm(death_not_remanded ~ 1 + l.state_pca + state_ut + as.factor(year), data = police.data.t1, family = "quasipoisson")
model.qp.cl = cl(police.data.t1, model.qp, police.data.t1$state_ut)

##Negative binominal
library(MASS)
model.nb = glm.nb(death_not_remanded ~ 1 + l.state_pca + state_ut + as.factor(year), data = police.data.t1)
model.nb.cl = cl(police.data.t1, model.nb, police.data.t1$state_ut)

## Delete new states
police.data.nn = subset(police.data.t1, police.data.t1$state_ut !=  "TELANGANA")
police.data.nn = police.data.nn[!police.data.nn$state_ut == "ANDHRA PRADESH", ]

police.data.nn$state_ut = as.factor(as.character(police.data.nn$state_ut))
levels(police.data.nn$state_ut)
length(police.data.nn$death_not_remanded)

model.poisson.nn = glm(death_not_remanded ~ 1 + l.state_pca + state_ut + as.factor(year), data = police.data.nn, family = "poisson")
model.poisson.nn.cl = cl(police.data.nn, model.poisson.nn, police.data.nn$state_ut)


## Delete MAHARASHTRA
police.data.nm = subset(police.data.t1, police.data.t1$state_ut !=  "MAHARASHTRA")

police.data.nm$state_ut = as.factor(as.character(police.data.nm$state_ut))
levels(police.data.nm$state_ut)
length(police.data.nm$death_not_remanded)

model.poisson.nm = glm(death_not_remanded ~ 1 + l.state_pca + state_ut + as.factor(year), data = police.data.nm, family = "poisson")
model.poisson.nm.cl = cl(police.data.nm, model.poisson.nm, police.data.nm$state_ut)

## Delete MAHARASHTRA and ANDHRA PRADESH
police.data.nma = subset(police.data.t1, police.data.t1$state_ut !=  "MAHARASHTRA")
police.data.nma = police.data.nma[!police.data.nma$state_ut == "ANDHRA PRADESH", ]

police.data.nma$state_ut = as.factor(as.character(police.data.nma$state_ut))
levels(police.data.nma$state_ut)
length(police.data.nma$death_not_remanded)

model.poisson.nma = glm(death_not_remanded ~ 1 + l.state_pca + state_ut + as.factor(year), data = police.data.nma, family = "poisson")
model.poisson.nma.cl = cl(police.data.nma, model.poisson.nma, police.data.nma$state_ut)

##Print
# ERROR 
#stargazer(model.poisson.b, model.qp.cl, model.nb.cl, model.poisson.nn.cl, model.poisson.nm.cl, model.poisson.nma.cl)

###############
###Table A12###
###############

## OLS Model with SHRC
model.ols.SHRCc = lm(death_not_remanded ~ 1 + l.state_pca + l.SHRC + state_ut + as.factor(year), data = police.data.t1)

model.ols.SHRCc.cl = cl(police.data.t1, model.ols.SHRCc, police.data.t1$state_ut)

## Poisson Model with SHRC
model.poisson.SHRCc = glm(death_not_remanded ~ 1 + l.state_pca + l.SHRC + state_ut + as.factor(year), data = police.data.t1, family = "poisson")

model.poisson.SHRCc.cl = cl(police.data.t1, model.poisson.SHRCc, police.data.t1$state_ut)

## OLS Model with SHRC
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add SHRC
  police.imp.1.l$l.SHRC = police.data.t1$l.SHRC
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded ~ 1 + l.state_pca + l.SHRC + gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.ols.SHRCc.cl.c = result.p.1

model.ols.SHRCc.cl.c[2:5, 1] = result.t3[, 1]
model.ols.SHRCc.cl.c[2:5, 2] = result.t3[, 2]
model.ols.SHRCc.cl.c[2:5, 4] = result.t3[, 3]

## Poisson Model with SHRC and controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add SHRC
  police.imp.1.l$l.SHRC = police.data.t1$l.SHRC
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + l.state_pca + l.SHRC + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.poisson.SHRCc.cl.c = result.p.1

model.poisson.SHRCc.cl.c[2:5, 1] = result.t3[, 1]
model.poisson.SHRCc.cl.c[2:5, 2] = result.t3[, 2]
model.poisson.SHRCc.cl.c[2:5, 4] = result.t3[, 3]

# ERROR 
#stargazer(model.ols.SHRCc.cl, model.ols.SHRCc.cl.c, model.poisson.SHRCc.cl, model.poisson.SHRCc.cl.c)


###############
###Table A13###
###############
## Add party_match to data
police.data.t1$party_match = police.data$party_match

## OLS Model with party
model.ols.party = lm(death_not_remanded ~ 1 + l.state_pca + party_match + state_ut + as.factor(year), data = police.data.t1)

model.ols.party.cl = cl(police.data.t1, model.ols.party, police.data.t1$state_ut)

## Poisson Model with party
model.poisson.party = glm(death_not_remanded ~ 1 + l.state_pca + party_match + state_ut + as.factor(year), data = police.data.t1, family = "poisson")

model.poisson.party.cl = cl(police.data.t1, model.poisson.party, police.data.t1$state_ut)

## OLS Model with party
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add party
  police.imp.1.l$party_match = police.data.t1$party_match
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded ~ 1 + l.state_pca + party_match + gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.ols.party.cl.c = result.p.1

model.ols.party.cl.c[2:5, 1] = result.t3[, 1]
model.ols.party.cl.c[2:5, 2] = result.t3[, 2]
model.ols.party.cl.c[2:5, 4] = result.t3[, 3]

## Poisson Model with party and controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add party
  police.imp.1.l$party_match = police.data.t1$party_match
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + l.state_pca + party_match + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.poisson.party.cl.c = result.p.1

model.poisson.party.cl.c[2:5, 1] = result.t3[, 1]
model.poisson.party.cl.c[2:5, 2] = result.t3[, 2]
model.poisson.party.cl.c[2:5, 4] = result.t3[, 3]

# ERROR 
#stargazer(model.ols.party.cl, model.ols.party.cl.c, model.poisson.party.cl, model.poisson.party.cl.c)


###############
###Table A14###
###############
## Add party_match_2006 to data
police.data.t1$party_match_2006 = police.data$party_match_2006

## OLS Model with party 2006
model.ols.party06 = lm(death_not_remanded ~ 1 + l.state_pca + party_match_2006 + state_ut + as.factor(year), data = police.data.t1)

model.ols.party06.cl = cl(police.data.t1, model.ols.party06, police.data.t1$state_ut)

## Poisson Model with party 2006
model.poisson.party06 = glm(death_not_remanded ~ 1 + l.state_pca + party_match_2006 + state_ut + as.factor(year), data = police.data.t1, family = "poisson")

model.poisson.party06.cl = cl(police.data.t1, model.poisson.party06, police.data.t1$state_ut)

## OLS Model with party 2006
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add party 2006
  police.imp.1.l$party_match_2006 = police.data.t1$party_match_2006
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded ~ 1 + l.state_pca + party_match_2006 + gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.ols.party06.cl.c = result.p.1

model.ols.party06.cl.c[2:5, 1] = result.t3[, 1]
model.ols.party06.cl.c[2:5, 2] = result.t3[, 2]
model.ols.party06.cl.c[2:5, 4] = result.t3[, 3]

## Poisson Model with party 2006 and controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add party 2006
  police.imp.1.l$party_match_2006 = police.data.t1$party_match_2006
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + l.state_pca + party_match_2006 + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.poisson.party06.cl.c = result.p.1

model.poisson.party06.cl.c[2:5, 1] = result.t3[, 1]
model.poisson.party06.cl.c[2:5, 2] = result.t3[, 2]
model.poisson.party06.cl.c[2:5, 4] = result.t3[, 3]

# ERROR 
#stargazer(model.ols.party06.cl, model.ols.party06.cl.c, model.poisson.party06.cl, model.poisson.party06.cl.c)

###############
###Table A15###
###############
## Add directives to police data
police.data.t1$ssc = police.data$ssc
police.data.t1$dgp_tenure = police.data$dgp_tenure
police.data.t1$o_tenure = police.data$o_tenure
police.data.t1$invest_law = police.data$invest_law
police.data.t1$peb = police.data$peb
police.data.t1$district_pca = police.data$district_pca

## Lagged ssc
police.data.t1 = ddply(police.data.t1, .(state_ut), transform, l.ssc = c(NA, ssc[-length(ssc)]))

## fill NA with 0
police.data.t1$l.ssc = ifelse(is.na(police.data.t1$l.ssc), 0, police.data.t1$l.ssc)

## Lagged dgp_tenure
police.data.t1 = ddply(police.data.t1, .(state_ut), transform, l.dgp_tenure = c(NA, dgp_tenure[-length(dgp_tenure)]))

## fill NA with 0
police.data.t1$l.dgp_tenure = ifelse(is.na(police.data.t1$l.dgp_tenure), 0, police.data.t1$l.dgp_tenure)

## Lagged o_tenure
police.data.t1 = ddply(police.data.t1, .(state_ut), transform, l.o_tenure = c(NA, o_tenure[-length(o_tenure)]))

## fill NA with 0
police.data.t1$l.o_tenure = ifelse(is.na(police.data.t1$l.o_tenure), 0, police.data.t1$l.o_tenure)

## Lagged invest_law
police.data.t1 = ddply(police.data.t1, .(state_ut), transform, l.invest_law = c(NA, invest_law[-length(invest_law)]))

## fill NA with 0
police.data.t1$l.invest_law = ifelse(is.na(police.data.t1$l.invest_law), 0, police.data.t1$l.invest_law)

## Lagged peb
police.data.t1 = ddply(police.data.t1, .(state_ut), transform, l.peb = c(NA, peb[-length(peb)]))

## fill NA with 0
police.data.t1$l.peb = ifelse(is.na(police.data.t1$l.peb), 0, police.data.t1$l.peb)

## Lagged district_pca
police.data.t1 = ddply(police.data.t1, .(state_ut), transform, l.district_pca = c(NA, district_pca[-length(district_pca)]))

## fill NA with 0
police.data.t1$l.district_pca = ifelse(is.na(police.data.t1$l.district_pca), 0, police.data.t1$l.district_pca)

## Directives correlations
## directives data
directives = police.data.t1[, c("l.ssc", "l.dgp_tenure", "l.o_tenure", "l.invest_law", "l.peb", "l.state_pca", "l.district_pca")]

stargazer(cor(directives))

###############
###Table A16###
###############

## OLS Models
model.ols.dis = lm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + state_ut + as.factor(year), data = police.data.t1)

model.ols.dis.cl = cl(police.data.t1, model.ols.dis, police.data.t1$state_ut)

model.ols.dir = lm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + l.ssc + l.dgp_tenure + l.o_tenure + l.invest_law + l.peb + state_ut + as.factor(year), data = police.data.t1)

model.ols.dir.cl = cl(police.data.t1, model.ols.dir, police.data.t1$state_ut)

model.ols.dir = lm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + l.ssc + l.dgp_tenure + l.o_tenure + l.invest_law + l.peb + l.state_pca * l.invest_law + state_ut + as.factor(year), data = police.data.t1)

model.ols.dir.i1.cl = cl(police.data.t1, model.ols.dir, police.data.t1$state_ut)

model.ols.dir = lm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + l.ssc + l.dgp_tenure + l.o_tenure + l.invest_law + l.peb + l.state_pca * l.o_tenure + state_ut + as.factor(year), data = police.data.t1)

model.ols.dir.i2.cl = cl(police.data.t1, model.ols.dir, police.data.t1$state_ut)

## Poisson Models
model.p.dis = glm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + state_ut + as.factor(year), data = police.data.t1, family = "poisson")

model.p.dis.cl = cl(police.data.t1, model.p.dis, police.data.t1$state_ut)

model.p.dir = glm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + l.ssc + l.dgp_tenure + l.o_tenure + l.invest_law + l.peb + state_ut + as.factor(year), data = police.data.t1, family = "poisson")

model.p.dir.cl = cl(police.data.t1, model.p.dir, police.data.t1$state_ut)

model.p.dir = glm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + l.ssc + l.dgp_tenure + l.o_tenure + l.invest_law + l.peb + l.state_pca * l.invest_law + state_ut + as.factor(year), data = police.data.t1, family = "poisson")

model.p.dir.i1.cl = cl(police.data.t1, model.p.dir, police.data.t1$state_ut)

model.p.dir = glm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + l.ssc + l.dgp_tenure + l.o_tenure + l.invest_law + l.peb + l.state_pca * l.o_tenure + state_ut + as.factor(year), data = police.data.t1, family = "poisson")

model.p.dir.i2.cl = cl(police.data.t1, model.p.dir, police.data.t1$state_ut)

# ERROR 
# stargazer(model.ols.dis.cl, model.ols.dir.cl, model.ols.dir.i1.cl,  model.ols.dir.i2.cl,
#           model.p.dis.cl,  model.p.dir.cl,  model.p.dir.i1.cl,  model.p.dir.i2.cl)


###############
###Table A17###
###############

## OLS Model with controls
## Dis
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add directives
  police.imp.1.l$l.ssc = police.data.t1$l.ssc
  police.imp.1.l$l.dgp_tenure = police.data.t1$l.dgp_tenure
  police.imp.1.l$l.o_tenure = police.data.t1$l.o_tenure
  police.imp.1.l$l.invest_law = police.data.t1$l.invest_law
  police.imp.1.l$l.peb = police.data.t1$l.peb
  police.imp.1.l$l.district_pca = police.data.t1$l.district_pca
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.ols.dis.cl.c = result.p.1

model.ols.dis.cl.c[2:5, 1] = result.t3[, 1]
model.ols.dis.cl.c[2:5, 2] = result.t3[, 2]
model.ols.dis.cl.c[2:5, 4] = result.t3[, 3]


## Dir
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add directives
  police.imp.1.l$l.ssc = police.data.t1$l.ssc
  police.imp.1.l$l.dgp_tenure = police.data.t1$l.dgp_tenure
  police.imp.1.l$l.o_tenure = police.data.t1$l.o_tenure
  police.imp.1.l$l.invest_law = police.data.t1$l.invest_law
  police.imp.1.l$l.peb = police.data.t1$l.peb
  police.imp.1.l$l.district_pca = police.data.t1$l.district_pca
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + l.ssc + l.dgp_tenure + l.o_tenure + l.invest_law + l.peb + gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:10, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:10, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 9, ncol = 3)

for (i in 1:9) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.ols.dir.cl.c = result.p.1

model.ols.dir.cl.c[2:10, 1] = result.t3[, 1]
model.ols.dir.cl.c[2:10, 2] = result.t3[, 2]
model.ols.dir.cl.c[2:10, 4] = result.t3[, 3]


## Dir.i1
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add directives
  police.imp.1.l$l.ssc = police.data.t1$l.ssc
  police.imp.1.l$l.dgp_tenure = police.data.t1$l.dgp_tenure
  police.imp.1.l$l.o_tenure = police.data.t1$l.o_tenure
  police.imp.1.l$l.invest_law = police.data.t1$l.invest_law
  police.imp.1.l$l.peb = police.data.t1$l.peb
  police.imp.1.l$l.district_pca = police.data.t1$l.district_pca
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + l.ssc + l.dgp_tenure + l.o_tenure + l.invest_law + l.peb + 
                  l.state_pca * l.invest_law +gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[c(2:10, 61), 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[c(2:10, 61), 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 10, ncol = 3)

for (i in 1:10) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.ols.dir.i1.cl.c = result.p.1

model.ols.dir.i1.cl.c[c(2:10, 61), 1] = result.t3[, 1]
model.ols.dir.i1.cl.c[c(2:10, 61), 2] = result.t3[, 2]
model.ols.dir.i1.cl.c[c(2:10, 61), 4] = result.t3[, 3]

## Dir.i2
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add directives
  police.imp.1.l$l.ssc = police.data.t1$l.ssc
  police.imp.1.l$l.dgp_tenure = police.data.t1$l.dgp_tenure
  police.imp.1.l$l.o_tenure = police.data.t1$l.o_tenure
  police.imp.1.l$l.invest_law = police.data.t1$l.invest_law
  police.imp.1.l$l.peb = police.data.t1$l.peb
  police.imp.1.l$l.district_pca = police.data.t1$l.district_pca
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + l.ssc + l.dgp_tenure + l.o_tenure + l.invest_law + l.peb + 
                  l.state_pca * l.o_tenure +gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[c(2:10, 61), 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[c(2:10, 61), 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 10, ncol = 3)

for (i in 1:10) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.ols.dir.i2.cl.c = result.p.1

model.ols.dir.i2.cl.c[c(2:10, 61), 1] = result.t3[, 1]
model.ols.dir.i2.cl.c[c(2:10, 61), 2] = result.t3[, 2]
model.ols.dir.i2.cl.c[c(2:10, 61), 4] = result.t3[, 3]


## Poisson Model with controls
## Dis
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add directives
  police.imp.1.l$l.ssc = police.data.t1$l.ssc
  police.imp.1.l$l.dgp_tenure = police.data.t1$l.dgp_tenure
  police.imp.1.l$l.o_tenure = police.data.t1$l.o_tenure
  police.imp.1.l$l.invest_law = police.data.t1$l.invest_law
  police.imp.1.l$l.peb = police.data.t1$l.peb
  police.imp.1.l$l.district_pca = police.data.t1$l.district_pca
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.p.dis.cl.c = result.p.1

model.p.dis.cl.c[2:5, 1] = result.t3[, 1]
model.p.dis.cl.c[2:5, 2] = result.t3[, 2]
model.p.dis.cl.c[2:5, 4] = result.t3[, 3]


## Dir
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add directives
  police.imp.1.l$l.ssc = police.data.t1$l.ssc
  police.imp.1.l$l.dgp_tenure = police.data.t1$l.dgp_tenure
  police.imp.1.l$l.o_tenure = police.data.t1$l.o_tenure
  police.imp.1.l$l.invest_law = police.data.t1$l.invest_law
  police.imp.1.l$l.peb = police.data.t1$l.peb
  police.imp.1.l$l.district_pca = police.data.t1$l.district_pca
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + l.ssc + l.dgp_tenure + l.o_tenure + l.invest_law + l.peb + gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:10, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:10, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 9, ncol = 3)

for (i in 1:9) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.p.dir.cl.c = result.p.1

model.p.dir.cl.c[2:10, 1] = result.t3[, 1]
model.p.dir.cl.c[2:10, 2] = result.t3[, 2]
model.p.dir.cl.c[2:10, 4] = result.t3[, 3]


## Dir.i1
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add directives
  police.imp.1.l$l.ssc = police.data.t1$l.ssc
  police.imp.1.l$l.dgp_tenure = police.data.t1$l.dgp_tenure
  police.imp.1.l$l.o_tenure = police.data.t1$l.o_tenure
  police.imp.1.l$l.invest_law = police.data.t1$l.invest_law
  police.imp.1.l$l.peb = police.data.t1$l.peb
  police.imp.1.l$l.district_pca = police.data.t1$l.district_pca
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + l.ssc + l.dgp_tenure + l.o_tenure + l.invest_law + l.peb + 
                  l.state_pca * l.invest_law +gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[c(2:10, 61), 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[c(2:10, 61), 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 10, ncol = 3)

for (i in 1:10) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.p.dir.i1.cl.c = result.p.1

model.p.dir.i1.cl.c[c(2:10, 61), 1] = result.t3[, 1]
model.p.dir.i1.cl.c[c(2:10, 61), 2] = result.t3[, 2]
model.p.dir.i1.cl.c[c(2:10, 61), 4] = result.t3[, 3]

## Dir.i2
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add directives
  police.imp.1.l$l.ssc = police.data.t1$l.ssc
  police.imp.1.l$l.dgp_tenure = police.data.t1$l.dgp_tenure
  police.imp.1.l$l.o_tenure = police.data.t1$l.o_tenure
  police.imp.1.l$l.invest_law = police.data.t1$l.invest_law
  police.imp.1.l$l.peb = police.data.t1$l.peb
  police.imp.1.l$l.district_pca = police.data.t1$l.district_pca
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + l.state_pca + l.district_pca + l.ssc + l.dgp_tenure + l.o_tenure + l.invest_law + l.peb + 
                  l.state_pca * l.o_tenure +gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[c(2:10, 61), 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[c(2:10, 61), 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 10, ncol = 3)

for (i in 1:10) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.p.dir.i2.cl.c = result.p.1

model.p.dir.i2.cl.c[c(2:10, 61), 1] = result.t3[, 1]
model.p.dir.i2.cl.c[c(2:10, 61), 2] = result.t3[, 2]
model.p.dir.i2.cl.c[c(2:10, 61), 4] = result.t3[, 3]

# ERROR 
# stargazer(model.ols.dis.cl.c, model.ols.dir.cl.c, model.ols.dir.i1.cl.c,  model.ols.dir.i2.cl.c,
#                     model.p.dis.cl.c,  model.p.dir.cl.c,  model.p.dir.i1.cl.c,  model.p.dir.i2.cl.c)

###############
###Table A18###
###############
## Add pca_bind
police.data.t1$pca_bind = police.data$pca_bind

## Lag pca_ind
police.data.t1 = ddply(police.data.t1, .(state_ut), transform, l.pca_bind = c(NA, pca_bind[-length(pca_bind)]))

## fill NA with 0
police.data.t1$l.pca_bind = ifelse(is.na(police.data.t1$l.pca_bind), 0, police.data.t1$l.pca_bind)

## Poisson no control binding
model.poisson.bind = glm(death_not_remanded ~ 1 + l.state_pca + l.pca_bind + as.factor(state_ut) + as.factor(year), data = police.data.t1, family = "poisson")

model.poisson.bind.cl = cl(police.data.t1, model.poisson.bind, police.data.t1$state_ut)

stargazer(model.poisson.bind.cl)

## categorical variable
police.data.t1$bindinglvl = ifelse(police.data.t1$l.pca_bind == 1 & police.data.t1$l.state_pca == 1, "Binding", 
                                    ifelse(police.data.t1$l.pca_bind == 0 & police.data.t1$l.state_pca == 1, "Regular", "No PCA"))
police.data.t1$bindinglvl = as.factor(police.data.t1$bindinglvl)
police.data.t1$bindinglvl = relevel(police.data.t1$bindinglvl, ref = "Regular")
levels(police.data.t1$bindinglvl)

## Poisson no control binding categorical
model.poisson.bind.ca = glm(death_not_remanded ~ 1 + bindinglvl + as.factor(state_ut) + as.factor(year), data = police.data.t1, family = "poisson")

model.poisson.bind.ca.cl = cl(police.data.t1, model.poisson.bind.ca, police.data.t1$state_ut)

stargazer(model.poisson.bind.ca.cl)


##############
###Table A19##
##############
## Add media_women
police.data.t1$media_women = police.data$media_women_0

police.data.media = police.data.t1[-which(is.na(police.data.t1$media_women)), ]

police.data.media$state_ut = as.factor(as.character(police.data.media$state_ut))

levels(police.data.media$state_ut)

## OLS Model with media_women
model.ols.media = lm(death_not_remanded ~ 1 + l.state_pca + media_women + state_ut + as.factor(year), data = police.data.media)

model.ols.media.cl = cl(police.data.media, model.ols.media, police.data.media$state_ut)

## Poisson Model with media_women
model.p.media = glm(death_not_remanded ~ 1 + l.state_pca + media_women + state_ut + as.factor(year), data = police.data.media, family = "poisson")

model.p.media.cl = cl(police.data.media, model.p.media, police.data.media$state_ut)

## OLS Model with media
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add media women
  police.imp.1.l$media_women = police.data.t1$media_women
  
  police.imp.1.l = police.imp.1.l[-which(is.na(police.imp.1.l$media_women)), ]
  
  police.imp.1.l$state_ut = as.factor(as.character(police.imp.1.l$state_ut))
  
  levels(police.imp.1.l)
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded ~ 1 + l.state_pca + media_women + gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.ols.media.cl.c = result.p.1

model.ols.media.cl.c[2:5, 1] = result.t3[, 1]
model.ols.media.cl.c[2:5, 2] = result.t3[, 2]
model.ols.media.cl.c[2:5, 4] = result.t3[, 3]

## Poisson Model with religion and controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add media women
  police.imp.1.l$media_women = police.data.t1$media_women
  
  police.imp.1.l = police.imp.1.l[-which(is.na(police.imp.1.l$media_women)), ]
  
  police.imp.1.l$state_ut = as.factor(as.character(police.imp.1.l$state_ut))
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + l.state_pca + media_women + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.p.media.cl.c = result.p.1

model.p.media.cl.c[2:5, 1] = result.t3[, 1]
model.p.media.cl.c[2:5, 2] = result.t3[, 2]
model.p.media.cl.c[2:5, 4] = result.t3[, 3]

# ERROR 
#stargazer(model.ols.media.cl, model.ols.media.cl.c, model.p.media.cl, model.p.media.cl.c)


##############
###Table A20##
##############
## Add literacy
police.data.t1$literacy = police.data$literacy

police.data.liter = police.data.t1[-which(is.na(police.data.t1$literacy)), ]

police.data.liter$state_ut = as.factor(as.character(police.data.liter$state_ut))

levels(police.data.liter$state_ut)

## OLS Model with literacy
model.ols.liter = lm(death_not_remanded ~ 1 + l.state_pca + literacy + state_ut + as.factor(year), data = police.data.liter)

model.ols.liter.cl = cl(police.data.liter, model.ols.liter, police.data.liter$state_ut)

## Poisson Model with literacy
model.p.liter = glm(death_not_remanded ~ 1 + l.state_pca + literacy + state_ut + as.factor(year), data = police.data.liter, family = "poisson")

model.p.liter.cl = cl(police.data.liter, model.p.liter, police.data.liter$state_ut)

## OLS Model with literacy
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add literacy
  police.imp.1.l$literacy = police.data.t1$literacy
  
  police.imp.1.l = police.imp.1.l[-which(is.na(police.imp.1.l$literacy)), ]
  
  police.imp.1.l$state_ut = as.factor(as.character(police.imp.1.l$state_ut))
  
  levels(police.imp.1.l)
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded ~ 1 + l.state_pca + literacy + gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.ols.liter.cl.c = result.p.1

model.ols.liter.cl.c[2:5, 1] = result.t3[, 1]
model.ols.liter.cl.c[2:5, 2] = result.t3[, 2]
model.ols.liter.cl.c[2:5, 4] = result.t3[, 3]

## Poisson Model with literacy and controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add literacy
  police.imp.1.l$literacy = police.data.t1$literacy
  
  police.imp.1.l = police.imp.1.l[-which(is.na(police.imp.1.l$literacy)), ]
  
  police.imp.1.l$state_ut = as.factor(as.character(police.imp.1.l$state_ut))
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + l.state_pca + literacy + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.p.liter.cl.c = result.p.1

model.p.liter.cl.c[2:5, 1] = result.t3[, 1]
model.p.liter.cl.c[2:5, 2] = result.t3[, 2]
model.p.liter.cl.c[2:5, 4] = result.t3[, 3]

# ERROR 
#stargazer(model.ols.liter.cl, model.ols.liter.cl.c, model.p.liter.cl, model.p.liter.cl.c)


##############
###Table A21##
##############

## Total Death
police.data.t1$total_death = police.data.t1$death_not_remanded + police.data.t1$death_remanded

## OLS Model with logged
police.data.t1$total_death_ln = log(police.data.t1$total_death + 1)

model.ols.total.l = lm(total_death_ln ~ 1 + l.state_pca + state_ut + as.factor(year), data = police.data.t1)

model.ols.total.l.cl = cl(police.data.t1, model.ols.total.l, police.data.t1$state_ut)

## OLS Model with SHRC
model.ols.total = lm(total_death ~ 1 + l.state_pca + state_ut + as.factor(year), data = police.data.t1)

model.ols.total.cl = cl(police.data.t1, model.ols.total, police.data.t1$state_ut)


## Poisson Model with SHRC
model.p.total = glm(total_death ~ 1 + l.state_pca + state_ut + as.factor(year), data = police.data.t1, family = "poisson")

model.p.total.cl = cl(police.data.t1, model.p.total, police.data.t1$state_ut)

## OLS Model with media
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Total Death
  police.imp.1.l$total_death = police.imp.1.l$death_not_remanded + police.imp.1.l$death_remanded
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(total_death ~ 1 + l.state_pca + gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.ols.total.cl.c = result.p.1

model.ols.total.cl.c[2:5, 1] = result.t3[, 1]
model.ols.total.cl.c[2:5, 2] = result.t3[, 2]
model.ols.total.cl.c[2:5, 4] = result.t3[, 3]

## Poisson Model with religion and controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Total Death
  police.imp.1.l$total_death = police.imp.1.l$death_not_remanded + police.imp.1.l$death_remanded
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(total_death ~ 1 + l.state_pca + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.p.total.cl.c = result.p.1

model.p.total.cl.c[2:5, 1] = result.t3[, 1]
model.p.total.cl.c[2:5, 2] = result.t3[, 2]
model.p.total.cl.c[2:5, 4] = result.t3[, 3]

# ERROR 
#stargazer(model.ols.total.cl, model.ols.total.cl.c, model.p.total.cl, model.p.total.cl.c)

##############
###Table A22##
##############
## OLS Model
model.ols.remanded = lm(death_remanded ~ 1 + l.state_pca + state_ut + as.factor(year), data = police.data.t1)

model.ols.remanded.cl = cl(police.data.t1, model.ols.remanded, police.data.t1$state_ut)


## Poisson Model
model.p.remanded = glm(death_remanded ~ 1 + l.state_pca + state_ut + as.factor(year), data = police.data.t1, family = "poisson")

model.p.remanded.cl = cl(police.data.t1, model.p.remanded, police.data.t1$state_ut)

## OLS Model
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_remanded ~ 1 + l.state_pca + gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.ols.remanded.cl.c = result.p.1

model.ols.remanded.cl.c[2:5, 1] = result.t3[, 1]
model.ols.remanded.cl.c[2:5, 2] = result.t3[, 2]
model.ols.remanded.cl.c[2:5, 4] = result.t3[, 3]

## Poisson Model with controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_remanded ~ 1 + l.state_pca + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))
  
}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])

}

result.t3

## Replace results to model result
model.p.remanded.cl.c = result.p.1

model.p.remanded.cl.c[2:5, 1] = result.t3[, 1]
model.p.remanded.cl.c[2:5, 2] = result.t3[, 2]
model.p.remanded.cl.c[2:5, 4] = result.t3[, 3]

# ERROR 
#stargazer(model.ols.remanded.cl, model.ols.remanded.cl.c, model.p.remanded.cl, model.p.remanded.cl.c)

###############
###Table A23###
###############

## Add pca_ind
police.data.t1$pca_ind = police.data$pca_ind

## Lag pca_ind
police.data.t1 = ddply(police.data.t1, .(state_ut), transform, l.pca_ind = c(NA, pca_ind[-length(pca_ind)]))

## fill NA with 0
police.data.t1$l.pca_ind = ifelse(is.na(police.data.t1$l.pca_ind), 0, police.data.t1$l.pca_ind)

## categorical variable
police.data.t1$indlvl = ifelse(police.data.t1$l.pca_ind == 1 & police.data.t1$l.state_pca == 1, "Ind", 
                                ifelse(police.data.t1$l.pca_ind == 0 & police.data.t1$l.state_pca == 1, "Regular", "No PCA"))
police.data.t1$indlvl = as.factor(police.data.t1$indlvl)
police.data.t1$indlvl = relevel(police.data.t1$indlvl, ref = "Regular")
levels(police.data.t1$indlvl)

## Poisson no control binding categorical
model.poisson.ind.ca = glm(death_not_remanded ~ 1 + indlvl + as.factor(state_ut) + as.factor(year), data = police.data.t1, family = "poisson")

model.poisson.ind.ca.cl = cl(police.data.t1, model.poisson.ind.ca, police.data.t1$state_ut)

stargazer(model.poisson.ind.ca.cl)

##############
###Table A24##
##############
## Add ngo
police.data.t1$ngo = police.data$ngo * 100

police.data.ngo = police.data.t1[-which(is.na(police.data.t1$ngo)), ]

police.data.ngo$state_ut = as.factor(as.character(police.data.ngo$state_ut))

levels(police.data.ngo$state_ut)

## OLS Model with ngo
model.ols.ngo = lm(death_not_remanded ~ 1 + l.state_pca + ngo + state_ut + as.factor(year), data = police.data.ngo)

model.ols.ngo.cl = cl(police.data.ngo, model.ols.ngo, police.data.ngo$state_ut)

## Poisson Model with ngo
model.p.ngo = glm(death_not_remanded ~ 1 + l.state_pca + ngo + state_ut + as.factor(year), data = police.data.ngo, family = "poisson")

model.p.ngo.cl = cl(police.data.ngo, model.p.ngo, police.data.ngo$state_ut)

## OLS Model with ngo
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add ngo
  police.imp.1.l$ngo = police.data$ngo * 100
  
  police.imp.1.l = police.imp.1.l[-which(is.na(police.imp.1.l$ngo)), ]
  
  police.imp.1.l$state_ut = as.factor(as.character(police.imp.1.l$state_ut))
  
  levels(police.imp.1.l)
  
  ## Poisson with outdata1.csv
  imp.1.p = lm(death_not_remanded ~ 1 + l.state_pca + ngo + gdp + 
                  head_trans + state_ut + 
                  as.factor(year), data = police.imp.1.l)
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))

}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])
  
}

result.t3

## Replace results to model result
model.ols.ngo.cl.c = result.p.1

model.ols.ngo.cl.c[2:5, 1] = result.t3[, 1]
model.ols.ngo.cl.c[2:5, 2] = result.t3[, 2]
model.ols.ngo.cl.c[2:5, 4] = result.t3[, 3]

## Poisson Model with ngo and controls
## Loop models for 5 imputation datasets
for (i in c(1:5)) {
  
  filename = paste("outdata", i, sep = "")
  filename.csv = paste(filename, "csv", sep = ".")
  police.imp.1 = read.csv(filename.csv)
  
  ## Lagged state_pca
  police.imp.1.l = ddply(police.imp.1, .(state_ut), transform, l.state_pca = c(NA, state_pca[-length(state_pca)]))
  
  ## fill NA with 0
  police.imp.1.l$l.state_pca = ifelse(is.na(police.imp.1.l$l.state_pca), 0, police.imp.1.l$l.state_pca)
  
  ## delete DAMAN & DIU 2001
  police.imp.1.l = police.imp.1.l[-500, ]
  
  ## Rescale GDP
  police.imp.1.l$gdp = police.imp.1.l$gdp / 1000000
  
  ## Add ngo
  police.imp.1.l$ngo = police.data$ngo * 100
  
  police.imp.1.l = police.imp.1.l[-which(is.na(police.imp.1.l$ngo)), ]
  
  police.imp.1.l$state_ut = as.factor(as.character(police.imp.1.l$state_ut))
  
  ## Poisson with outdata1.csv
  imp.1.p = glm(death_not_remanded ~ 1 + l.state_pca + ngo + gdp + 
                   head_trans + state_ut + 
                   as.factor(year), data = police.imp.1.l, family = "poisson")
  
  result.p.1 = cl(police.imp.1.l, imp.1.p, police.imp.1.l$state_ut)
  
  nam.e = paste("e", i, sep = "")
  assign(nam.e, result.p.1[2:5, 1])
  
  nam.se = paste("se", i, sep = "")
  assign(nam.se, result.p.1[2:5, 2])
  
}

beta.t = cbind(e1, e2, e3, e4, e5)

beta.se = cbind(se1, se2, se3, se4, se5)

## Calculate imputed beta and SEs
se_calc = function(q, se) {
  
  part1 = sum((se)^2) / length(se)
  part2 = sum((q - mean(q))^2) / (length(q) - 1) * (1 + 1 / length(q))
  se.imp = sqrt(part1 + part2)
  q.imp = mean(q)
  p.value = 2 * pnorm(abs(q.imp / se.imp), lower.tail = FALSE)
  
  return(c(q.imp, se.imp, p.value))

}

## Print poisson results
result.t3 = matrix(NA, nrow = 4, ncol = 3)

for (i in 1:4) {
  
  result.t3[i, ] = se_calc(q = beta.t[i, ], se = beta.se[i, ])

}

result.t3

## Replace results to model result
model.p.ngo.cl.c = result.p.1

model.p.ngo.cl.c[2:5, 1] = result.t3[, 1]
model.p.ngo.cl.c[2:5, 2] = result.t3[, 2]
model.p.ngo.cl.c[2:5, 4] = result.t3[, 3]

stargazer(model.ols.ngo.cl, model.ols.ngo.cl.c, model.p.ngo.cl, model.p.ngo.cl.c)

