# Template for survival analysis with Cox mixed effects model.
# Setting up the analysis ----
install.packages("dplyr") #Installs the dplyr program for data cleanup
install.packages("coxme") #Installs the Cox Mixed Effect Model
install.packages("survminer") #Installs the survminer program for ggplot survival curves plotting

library(dplyr)
library(coxme)
library(survminer)

#Renames the data files, replace placeholder with csv file name 

DCV <- placeholder

FHV <- placeholder

#Filters out the PBS samples to stop them from confounding the analysis

DCV_noPBS <- filter(DCV, Virus == 1)

FHV_noPBS <- filter(FHV, Virus == 1)

#test for proportional hazards assumption
#coxph assumes that hazard ratios remain the same over time
#if the chi-square p-value < 0.05, then assumption is invalid and coxph cannot be used
#kaplan-meier might be a suitable alternative (it has no assumptions about the data)

fitCoxphDCV <- coxph(Surv(Days, Survival)~ Mutant+Virus, data = DCV_noPBS)
ftestDCV <- cox.zph(fitCoxphDCV)
ftestDCV

fitCoxphFHV <- coxph(Surv(Days, Survival)~ Mutant+Virus, data = FHV_noPBS)
ftestFHV <- cox.zph(fitCoxphFHV)
ftestFHV

#Fits the coxme to the survival data. 
#Where the fixed effects Virus and the random effects 
#are those of the technical and biological replications
#Tests if the survival curves of miRNA KO mutants are 
#significantly different  to the wild-type. 
#Only uses the data of virus-infected flies. 

fitDCV <- coxme(Surv(Days, Survival)~ Mutant + (1|BiolRep/TechRep), 
                data = DCV_noPBS)
fitDCV

fitFHV <- coxme(Surv(Days, Survival)~ Mutant + (1|BiolRep/TechRep), 
                data = FHV_noPBS)
fitFHV

#Effect of mutant and virus infection on mortality 
#(simpler model, no interactions between mutant and virus)
#fitmiRNA_Virus1 <- coxme(Surv(Days, Survival)~ Mutant+Virus + (1|BiolRep/TechRep), data=FHV)
#fitmiRNA_Virus1

#complex model
#fitmiRNA_Virus2 <- coxme(Surv(Days, Survival)~ Mutant*Virus + (1|BiolRep/TechRep), data=FHV)
#fitmiRNA_Virus2

#Fits survival curves based on the survival data
#Wild-type vs mutants
miRNAcurves <- survfit(Surv(Days, Survival) ~ Mutant, data=DCV_noPBS)
ggsurvplot(miRNAcurves,
           data = DCV_noPBS,               #change to the virus under analysis
           axes.offset = FALSE,            #changes axes to start at zero
           line = c(1,1),
           pval = placeholder,                   #change to coxme pvalue
           xlab = "Time in days",
           ylab = "Proportional survival",
           break.time.by = 1, 
           break.y.by = 0.2,
           ggtheme = theme_classic(),
           legend.title = "Genotype",      #change legend title to suit your analysis
           legend.labs = 
             c("w1118", "miR- KO")         #change the legends to the miRNA
)

# Export coxme output
sink("miR-_Virus_coxme.txt") # add miRNA name and virus | creates the txt file in WD
print(fitDCV) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

#save the plot as miRNA_Virus in the folder of your choice
#use the export function to save as pdf
#pdf -> device size -> preview
#preview before saving to make sure it is correct.
