# Tutorial on parametric survival analysis in R ----
## Using the Cox mixed effects model has revealed that a lot of the miRNA survival data is not within the proportional hazard assumptions. A paper from Kutzer, Kurtz and Armitage (2019) used the parametric regression model to support analyses performed using the Cox mixed effects model.
## Tutorial here is using the Weibull regression in the survreg() package and adapted from  Zhang (2016) and Incerti (2019).

# Load libraries ----
library(survival)

## CoxPH model for comparison
cox.lung <- coxph(Surv(time, status) ~ ph.ecog + sex + age, data = lung)
summary(cox.lung)
cox.zph(cox.lung)

## Weibull distribution
wei.lung <- survreg(Surv(time, status) ~ ph.ecog + sex + age, data = lung, dist = 'weibull')
summary(wei.lung)

## Exponential distribution
expo.lung <- survreg(Surv(time, status) ~ ph.ecog + sex + age, data = lung, dist = 'exponential')
summary(expo.lung)

anova(wei.lung, expo.lung)
