# set working directory
setwd('E:/GitHub/APC-of-Another-3-Kinds-of-Models')

# load packages
library(tidyverse); library(lubridate); library(RColorBrewer)
library(mgcv); library(splines); library(dlnm)

# load data
load('Data/dat.RData')

# linear model with splines -----------------------------------------------

df.day <- 21 # degrees of freedom for 'day'
df.x5 <- 3 # degrees of freedom for 'x6'

m1 <- glm(y ~ x1 + as.factor(x6) + ns(day, df = df.day) + ns(x5, df = df.x5), family = poisson, data = dat)
m2 <- gam(y ~ x1 + as.factor(x6) + s(day, bs='cr', fx=TRUE, k=df.day+1) + s(x5, bs='cr', fx=TRUE, k=df.x5+1), family = poisson, data = dat)
m3 <- gam(y ~ x1 + as.factor(x6) + s(day, bs='cs', fx=TRUE, k=df.day+1) + s(x5, bs='cs', fx=TRUE, k=df.x5+1), family = poisson, data = dat)

# relative risk of y per 1 unit change in x1
RR1 <- exp(coef(m1)[2] * 1)

# APC 
## 1. get pairs
u <- 'x1'
v <- c('x5', 'x6', 'day')
XA <- dat %>% select(u, v) %>% mutate(No.A = 1:nrow(.)) # 'No.A' = 'OriginalRowNumber' in GetPairs() of predcomps package
XB <- dat %>% select(u, v) %>% mutate(No.B = 1:nrow(.)) # 'No.B' = 'OriginalRowNumber.B'

## 2. calculate distance
vMatrixA <- as.matrix(XA[v])
vMatrixB <- as.matrix(XB[v])

covV <- cov(vMatrixB)
distMatrix <- apply(vMatrixA, 1, function(row) mahalanobis(vMatrixB, row, covV))

DistDat <- distMatrix %>% as.tibble %>% 
  mutate(No.B = XB$No.B) %>% 
  gather(key = 'No.A', value = 'MahalanobisDistance', -No.B) %>% 
  mutate(No.A = as.integer(str_sub(No.A, 2, -1)))

## 3. calculate weights
mahalanobisConstantTerm <- 1
pairs <- DistDat %>% 
  left_join(., XA, by = 'No.A') %>% 
  left_join(., XB, by = 'No.B', suffix = c(".A", ".B")) %>% 
  mutate(WeightInPair.raw = 1/(mahalanobisConstantTerm + MahalanobisDistance)) %>% 
  group_by(No.A) %>% 
  mutate(WeightInPair.normalized = WeightInPair.raw/sum(WeightInPair.raw)) %>% 
  ungroup

## 4. get fitted values
comp0 <- c(paste0(u, '.A'), paste0(v, '.A'))
comp1 <- c(paste0(u, '.B'), paste0(v, '.A'))

Xcomp0 <- pairs %>% select(comp0) %>% rename_all(., function(x) str_sub(x, 1, -3))
Xcomp1 <- pairs %>% select(comp1) %>% rename_all(., function(x) str_sub(x, 1, -3))

yhat0.m1 <- predict(m1, Xcomp0, type = 'response')
yhat1.m1 <- predict(m1, Xcomp1, type = 'response')

yhat0.m2 <- predict(m2, Xcomp0, type = 'response')
yhat1.m2 <- predict(m2, Xcomp1, type = 'response')

yhat0.m3 <- predict(m3, Xcomp0, type = 'response')
yhat1.m3 <- predict(m3, Xcomp1, type = 'response')

summary(yhat1.m1)
summary(yhat1.m2)
summary(yhat1.m3)

## 5. compute APC from pairs
uDiff <- pairs[comp1[1]] - pairs[comp0[1]]
yDiff <- yhat1.m1 - yhat0.m1
weight <- pairs$WeightInPair.normalized

numer <- sum(weight * yDiff * sign(uDiff), na.rm = TRUE)
denom <- sum(weight * uDiff * sign(uDiff), na.rm = TRUE)
Apc <-  numer / denom

### TO SOLVE
# How to caculate distance of vectors with binary element
# How to get fitted values from model with splines.

# Linear model with interaction -------------------------------------------

df.day <- 21 # degrees of freedom for 'day'
df.x5 <- 3 # degrees of freedom for 'x6'

m1 <- glm(y ~ x1 + as.factor(x6) + ns(day, df = df.day) + x1 : ns(x5, df = df.x5), family = poisson, data = dat)
m2 <- gam(y ~ x1 + as.factor(x6) + s(day, bs='cr', fx=TRUE, k=df.day+1) + s(x5, bs='cr', fx=TRUE, k=df.x5+1, by=x1), family = poisson, data = dat)
m3 <- gam(y ~ x1 + as.factor(x6) + s(day, bs='cs', fx=TRUE, k=df.day+1) + s(x5, bs='cs', fx=TRUE, k=df.x5+1, by=x1), family = poisson, data = dat)

m2 <- gam(y ~ x1 + as.factor(x6) + 
            s(day, bs='cr', fx=TRUE, k=df.day+1) +
            te(x5, bs='cr', fx=TRUE, k=df.x5+1, by=x1), family = poisson, data = dat)


# relative risk of y per 1 unit change in x1
RR1 <- exp(coef(m1)[2] * 1)

# APC 
## 1. get pairs
u <- 'x1'
v <- c('x5', 'x6', 'day')
XA <- dat %>% select(u, v) %>% mutate(No.A = 1:nrow(.)) # 'No.A' = 'OriginalRowNumber' in GetPairs() of predcomps package
XB <- dat %>% select(u, v) %>% mutate(No.B = 1:nrow(.)) # 'No.B' = 'OriginalRowNumber.B'

## 2. calculate distance
vMatrixA <- as.matrix(XA[v])
vMatrixB <- as.matrix(XB[v])

covV <- cov(vMatrixB)
distMatrix <- apply(vMatrixA, 1, function(row) mahalanobis(vMatrixB, row, covV))

DistDat <- distMatrix %>% as.tibble %>% 
  mutate(No.B = XB$No.B) %>% 
  gather(key = 'No.A', value = 'MahalanobisDistance', -No.B) %>% 
  mutate(No.A = as.integer(str_sub(No.A, 2, -1)))

## 3. calculate weights
mahalanobisConstantTerm <- 1
pairs <- DistDat %>% 
  left_join(., XA, by = 'No.A') %>% 
  left_join(., XB, by = 'No.B', suffix = c(".A", ".B")) %>% 
  mutate(WeightInPair.raw = 1/(mahalanobisConstantTerm + MahalanobisDistance)) %>% 
  group_by(No.A) %>% 
  mutate(WeightInPair.normalized = WeightInPair.raw/sum(WeightInPair.raw)) %>% 
  ungroup

## 4. get fitted values
comp0 <- c(paste0(u, '.A'), paste0(v, '.A'))
comp1 <- c(paste0(u, '.B'), paste0(v, '.A'))

Xcomp0 <- pairs %>% select(comp0) %>% rename_all(., function(x) str_sub(x, 1, -3))
Xcomp1 <- pairs %>% select(comp1) %>% rename_all(., function(x) str_sub(x, 1, -3))

yhat0.m1 <- predict(m1, Xcomp0, type = 'response')
yhat1.m1 <- predict(m1, Xcomp1, type = 'response')

yhat0.m2 <- predict(m2, Xcomp0, type = 'response')
yhat1.m2 <- predict(m2, Xcomp1, type = 'response')

yhat0.m3 <- predict(m3, Xcomp0, type = 'response')
yhat1.m3 <- predict(m3, Xcomp1, type = 'response')

summary(yhat1.m1)
summary(yhat1.m2)
summary(yhat1.m3)

## 5. compute APC from pairs
uDiff <- pairs[comp1[1]] - pairs[comp0[1]]
yDiff <- yhat1.m1 - yhat0.m1
weight <- pairs$WeightInPair.normalized

numer <- sum(weight * yDiff * sign(uDiff), na.rm = TRUE)
denom <- sum(weight * uDiff * sign(uDiff), na.rm = TRUE)
Apc <-  numer / denom
