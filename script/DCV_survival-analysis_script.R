# Load packages ----

library(dplyr)
library(coxme)
library(survminer)

# Rename the data files, replace placeholder with csv file name ----

DCV <- placeholder

# Filter out the PBS samples to stop them from confounding the analysis ----

DCV_noPBS <- filter(DCV, Virus == 1)

# Test for proportional hazards assumption ----
# coxph assumes that hazard ratios remain the same over time
# if the chi-square p-value < 0.05, then assumption is invalid and coxph cannot be used
# kaplan-meier might be a suitable alternative (it has no assumptions about the data)

fitCoxphDCV <- coxph(Surv(Days, Survival)~ Mutant+Virus, data = DCV_noPBS)
ftestDCV <- cox.zph(fitCoxphDCV)
ftestDCV
ggcoxzph(ftestDCV)

# Fit the Cox mixed effects model to the survival data ----
# Where the fixed effects Virus and the random effects 
# are those of the technical and biological replications
# Tests if the survival curves of miRNA KO mutants are 
# significantly different  to the wild-type. 
# Only uses the data of virus-infected flies. 

fitDCV <- coxme(Surv(Days, Survival)~ Mutant + (1|BiolRep/TechRep), 
                data = DCV_noPBS)
fitDCV

# Export coxme output ----
sink("miR-_Virus_coxme.txt") # add miRNA name and virus | creates the txt file in WD
print(fitDCV) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

# Fit survival curves ----
miRNAcurves <- survfit(Surv(Days, Survival) ~ Mutant, data=DCV_noPBS)
ggsurvplot(miRNAcurves,
           data = DCV_noPBS,               #change to the virus under analysis
           axes.offset = FALSE,            #changes axes to start at zero
           line = c(1,1),
           pval = placeholder,             #change to coxme pvalue
           xlab = "Time in days",
           ylab = "Proportional survival",
           break.time.by = 1, 
           break.y.by = 0.2,
           ggtheme = theme_classic(),
           legend.title = "Genotype",      #change legend title to suit your analysis
           legend.labs = 
             c("w1118", "miR- KO")         #change the legends to the miRNA
)