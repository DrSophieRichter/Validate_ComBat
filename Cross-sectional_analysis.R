knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)

library(dplyr)
library(tidyr)
library(janitor)
library(ggplot2)
library(lsr)
library(knitr)
library(kableExtra)
library(stringr)
library(neuroCombat)
library(lme4)
library(lmerTest)
library(ggpubr)
library(rstatix)
library(flextable)
library(magrittr)
library(officer)
library(reticulate)
source("Functions/My_table_funs.R")
source("Functions/My_table_funs_corthick.R")
source("Functions/Select_struc_funs.R")
source("Functions/Select_dti_funs.R")
source("Functions/longCombat.R")
source("Functions/Harmonize_funs.R")
source("Functions/Simulate_funs.R")
## #load required python packages

## from neuroHarmonize import harmonizationLearn

## import pandas as pd

## import numpy as np


#Load data
scans <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Legacy_data/Most_uptodate_curated_data/Sophies_scan_database_20220121.csv", na.strings = c("", " ", "NA"))
qc1 <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Analysis_Aug2021/Current_data/Check_these_scans_1.csv", na.strings = c("", " ", "NA")) %>%
  select(Scan_ID, Decision)
qc2 <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Analysis_Aug2021/Current_data/Check_these_scans_2.csv", na.strings = c("", " ", "NA")) %>%
  select(Scan_ID, Decision)
qc3 <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Analysis_Aug2021/Current_data/Check_these_scans_3.csv", na.strings = c("", " ", "NA")) %>%
  select(Scan_ID, Decision)
qc_bad <- rbind(qc1, qc2, qc3) %>% as.data.frame() %>% filter(!Decision == "keep")

vol <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Raw_data/Sophies_processed_data/malpem138_volume_sophie.csv")
mymalp <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Analysis_Aug2021/Current_data/MALPEM_parcellation_Sophie.csv") %>% remove_empty()

corthick <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Raw_data/Sophies_processed_data/malpem021_corthick_sophie.csv", na.strings = c("", "NA"))
add <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Raw_data/Sophies_processed_data/malpem005_corthick_sophie.csv", na.strings = c("", " "))

dti <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Raw_data/Sophies_processed_data/malpem007_dti_calc_Sophie.csv")

clinical <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Legacy_data/Most_uptodate_curated_data/Sophies_clinical_database_20220121.csv", na.strings = c("", " ", "NA"))

#Add ICV as a covariate to the scans spreadsheet
icv <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Raw_data/Sophies_processed_data/malpem138_volume_sophie.csv", na.strings = c("", " ", "NA"))
icv <- icv %>% select(Scan_ID, ICV = TotalBrain)
scans <- merge(scans, icv, by = "Scan_ID", all.x = TRUE, all.y = FALSE)

#Prepare structural data

#collapse metrics into desired ROIs
vol         <- collapse_vol(vol, mymalp)
corthick    <- collapse_corthick(corthick, add)

#Select eligible scans for structural analysis
struc_within    <- select_struc_within(scans, qc_bad) #same scanner rescan
struc_across    <- select_struc_across(scans, qc_bad) #different scanner rescan

#Add volume
slong_vol   <- add_vol_within(struc_within, vol)
sshort_vol  <- calculate_CoV(slong_vol)#calculate CoV for each subject and ROI

tlong_vol   <- add_vol_across(struc_across, vol)
tshort_vol  <- calculate_CoV(tlong_vol) #make a short version of the table with perc difference

#Add cortical thickness
slong_corthick   <- add_corthick_within(struc_within, corthick)
sshort_corthick  <- calculate_CoV_corthick(slong_corthick)#make a short version of the table with perc difference

tlong_corthick   <- add_corthick_across(struc_across, corthick)
tshort_corthick  <- calculate_CoV_corthick(tlong_corthick) #make a short version of the table with perc difference

#combine all scans for harmonisation and add clinical covariates
corthick_for_combat <- all_struc_scans(corthick, slong_corthick, tlong_corthick, clinical)
vol_for_combat      <- all_struc_scans(vol, slong_vol, tlong_vol, clinical)

#Prepare DTI data

#extract fa or md and DTI series name
md <- extract_dti_fresh(dti,metric = "md")
fa <- extract_dti_fresh(dti,metric = "fa")

#Select eligible scans and combine with md
slong_md  <- select_dti_within(md, scans, qc_bad)
sshort_md <- calculate_CoV(slong_md) #make a short version of the table with perc difference

tlong_md  <- select_dti_across(md, scans, qc_bad)
tshort_md <- calculate_CoV(tlong_md) #make a short version of the table with perc difference

#Select eligible scans and select with fa
slong_fa  <- select_dti_within(fa, scans, qc_bad)
sshort_fa <- calculate_CoV(slong_fa) #make a short version of the table with perc difference

tlong_fa  <- select_dti_across(fa, scans, qc_bad)
tshort_fa <- calculate_CoV(tlong_fa) #make a short version of the table with perc difference

#combine all scans for harmonisation and add clinical covariates
md_for_combat <- all_dti_scans(md, slong_md, tlong_md, clinical)
fa_for_combat <- all_dti_scans(fa, slong_fa, tlong_fa, clinical)

#Make summary table of cohort characteristics

####for within-scanner cohort####
#structural
within_cohort_struc <- make_cohort_table_1(sshort_vol, slong_vol, scans, clinical)
within_cohort_struc

#dti
within_cohort_dti <- make_cohort_table_1(sshort_fa, slong_fa, scans, clinical)
within_cohort_dti

####for across-scanner cohort####
#structural
across_cohort_struc <- make_cohort_table_1(tshort_vol, tlong_vol, scans, clinical)
across_cohort_struc

#dti
across_cohort_dti <- make_cohort_table_1(tshort_fa, tlong_fa, scans, clinical)
across_cohort_dti

all_cohorts <- cbind(within_cohort_struc, within_cohort_dti, across_cohort_struc, across_cohort_dti)
colnames(all_cohorts) <- c("within_cohort_struc", "within_cohort_dti", "across_cohort_struc", "across_cohort_dti")
all_cohorts

write.table(all_cohorts, file = "Figures/Cohorts.txt", sep = ",", quote = FALSE, row.names = T)
#summarise scanners
scanners_struc1 <- slong_vol %>% ungroup %>% select(Scanner)
scanners_struc2 <- tlong_vol %>% ungroup %>% select(Scanner)
scanners_struc <- rbind(scanners_struc1, scanners_struc2) %>% as.data.frame()
summary(as.factor(scanners_struc$Scanner))

scanners_struc1 <- slong_fa %>% ungroup %>% select(Scanner)
scanners_struc2 <- tlong_fa %>% ungroup %>% select(Scanner)
scanners_struc <- rbind(scanners_struc1, scanners_struc2) %>% as.data.frame()
summary(as.factor(scanners_struc$Scanner))

#check dates of acquisition
temp <- scans %>% filter(Scan_ID %in% slong_vol$Scan_ID |
                           Scan_ID %in% tlong_vol$Scan_ID |
                           Scan_ID %in% slong_fa$Scan_ID |
                           Scan_ID %in% tlong_fa$Scan_ID)
temp <- temp %>% select(Scan_ID) %>% filter(!str_detect(Scan_ID, "Con"))
temp <- temp %>% separate(Scan_ID, into = c("wbic", "date", "uid"), sep = "_")


#Harmonise

#Standard crosssectional combat on all available scans
#followed by selecting only scans from the across-scan cohort

crosscombat_vol      <- 
  get_crosscombat(vol_for_combat, from = "Ventricles", to = "Brainstem") %>%
  extract_across(., tlong_vol)

crosscombat_corthick <- 
  get_crosscombat(corthick_for_combat, from = "Frontal", to = "WholeCortex") %>%
  extract_across(., tlong_corthick, is_corthick = "yes")

crosscombat_md       <- 
  get_crosscombat(md_for_combat, from = "Ventricles", to = "Brainstem") %>%
  extract_across(., tlong_md)

crosscombat_fa       <- 
  get_crosscombat(fa_for_combat, from = "Ventricles", to = "Brainstem") %>%
  extract_across(., tlong_fa)

#Longitudinal combat - on all available scans
#followed by selecting only scans from the across-scan cohort

longcombat_vol       <- 
  get_longcombat(vol_for_combat) %>%
  extract_across(., tlong_vol)

longcombat_corthick  <- 
  get_longcombat(corthick_for_combat, is_corthick = "yes") %>%
  extract_across(., tlong_corthick, is_corthick = "yes")

longcombat_md        <- 
  get_longcombat(md_for_combat) %>%
  extract_across(., tlong_md)

longcombat_fa        <- 
  get_longcombat(fa_for_combat) %>%
  extract_across(., tlong_fa)


#Combat-GAM - on all available scans
#followed by selecting only scans from the across-scan cohort
gamcombat_vol      <- 
  get_gamcombat(vol_for_combat, from = "Ventricles", to = "Brainstem") %>%
  extract_across(., tlong_vol)

gamcombat_corthick <- 
  get_gamcombat(corthick_for_combat, from = "Frontal", to = "WholeCortex") %>%
  extract_across(., tlong_corthick, is_corthick = "yes")

gamcombat_md       <- 
  get_gamcombat(md_for_combat, from = "Ventricles", to = "Brainstem") %>%
  extract_across(., tlong_md)

gamcombat_fa       <- 
  get_gamcombat(fa_for_combat, from = "Ventricles", to = "Brainstem") %>%
  extract_across(., tlong_fa)



#make a summary table for harmonisation results

#volume results
t1a <- make_summary_table(sshort_vol, sshort_vol)
t1b <- make_summary_table(tshort_vol, sshort_vol)
t1c <- make_summary_table(crosscombat_vol, sshort_vol)
t1d <- make_summary_table(longcombat_vol, sshort_vol)
t1e <- make_summary_table(gamcombat_vol, sshort_vol)
t1  <- cbind(t1a[,1:4], t1b[,-c(1, 8)], t1c[,-c(1, 2, 8)], t1d[,-c(1, 2, 8)], t1e[,-c(1, 2, 8)])

#cortical thickness
t2a <- make_summary_table_corthick(sshort_corthick, sshort_corthick)
t2b <- make_summary_table_corthick(tshort_corthick, sshort_corthick)
t2c <- make_summary_table_corthick(crosscombat_corthick, sshort_corthick)
t2d <- make_summary_table_corthick(longcombat_corthick, sshort_corthick)
t2e <- make_summary_table_corthick(gamcombat_corthick, sshort_corthick)
t2  <- cbind(t2a[,1:4], t2b[,-c(1, 8)], t2c[,-c(1, 2, 8)], t2d[,-c(1, 2, 8)], t2e[,-c(1, 2, 8)])

#md
t3a <- make_summary_table(sshort_md, sshort_md)
t3b <- make_summary_table(tshort_md, sshort_md)
t3c <- make_summary_table(crosscombat_md, sshort_md)
t3d <- make_summary_table(longcombat_md, sshort_md)
t3e <- make_summary_table(gamcombat_md, sshort_md)
t3  <- cbind(t3a[,1:4], t3b[,-c(1, 8)], t3c[,-c(1, 2, 8)], t3d[,-c(1, 2, 8)], t3e[,-c(1, 2, 8)])

#fa
t4a <- make_summary_table(sshort_fa, sshort_fa)
t4b <- make_summary_table(tshort_fa, sshort_fa)
t4c <- make_summary_table(crosscombat_fa, sshort_fa)
t4d <- make_summary_table(longcombat_fa, sshort_fa)
t4e <- make_summary_table(gamcombat_fa, sshort_fa)
t4  <- cbind(t4a[,1:4], t4b[,-c(1, 8)], t4c[,-c(1, 2, 8)], t4d[,-c(1, 2, 8)], t4e[,-c(1, 2, 8)])

tall <- rbind(t1, t2, t3, t4)

colnames(tall) <- c("ROI" , "Scan_pairs_within", "ICC_within", "Perc_change_within",
                    "Scan_pairs_across",
                    "ICC_unharm", "Perc_change_unharm", "Delta_unharm", "Praw_unharm", "Padj_unharm", "Size_unharm",
                    "ICC_cross", "Perc_change_cross", "Delta_cross","Praw_cross", "Padj_cross", "Size_cross",
                    "ICC_long", "Perc_change_long", "Delta_long","Praw_long", "Padj_long", "Size_long",
                     "ICC_gam", "Perc_change_gam", "Delta_gam","Praw_gam", "Padj_gam", "Size_gam")

#format table for word
ft <- make_flextable(tall)

#save table to word
sect_properties <- prop_section(
  page_size = page_size(orient = "landscape"),
  type = "continuous",
  page_margins = page_mar()
)

save_as_docx(ft, path = "Figures/Cross-sec-CoV_table.docx", pr_section = sect_properties)



#As a sensitivity analysis for whether the 
# parametric prior of CrossCombat is appropriate
# I will replicat the above table but using
# a crosscombat_nonpara and crosscombat_nonbays
# as the two harmonization methods

####Harmonise

#Non-parametric prior
crosscombat_nonpara_vol      <- 
  get_crosscombat(vol_for_combat, from = "Ventricles", to = "Brainstem", prior = "non-parametric") %>%
  extract_across(., tlong_vol)

crosscombat_nonpara_corthick <- 
  get_crosscombat(corthick_for_combat, from = "Frontal", to = "WholeCortex", prior = "non-parametric") %>%
  extract_across(., tlong_corthick, is_corthick = "yes")

crosscombat_nonpara_md       <- 
  get_crosscombat(md_for_combat, from = "Ventricles", to = "Brainstem", prior = "non-parametric") %>%
  extract_across(., tlong_md)

crosscombat_nonpara_fa       <- 
  get_crosscombat(fa_for_combat, from = "Ventricles", to = "Brainstem", prior = "non-parametric") %>%
  extract_across(., tlong_fa)

#No prior (location shift model)
crosscombat_nonbays_vol      <- 
  get_crosscombat(vol_for_combat, from = "Ventricles", to = "Brainstem", prior = "none") %>%
  extract_across(., tlong_vol)

crosscombat_nonbays_corthick <- 
  get_crosscombat(corthick_for_combat, from = "Frontal", to = "WholeCortex", prior = "none") %>%
  extract_across(., tlong_corthick, is_corthick = "yes")

crosscombat_nonbays_md       <- 
  get_crosscombat(md_for_combat, from = "Ventricles", to = "Brainstem", prior = "none") %>%
  extract_across(., tlong_md)

crosscombat_nonbays_fa       <- 
  get_crosscombat(fa_for_combat, from = "Ventricles", to = "Brainstem", prior = "none") %>%
  extract_across(., tlong_fa)

#make a summary table for harmonisation results

#volume results
t1f <- make_summary_table(crosscombat_nonpara_vol, sshort_vol)
t1g <- make_summary_table(crosscombat_nonbays_vol, sshort_vol)
s1  <- cbind(t1a[,1:4], t1b[,-c(1, 8)], t1f[,-c(1, 2, 8)], t1g[,-c(1, 2, 8)])

#cortical thickness
t2f <- make_summary_table_corthick(crosscombat_nonpara_corthick, sshort_corthick)
t2g <- make_summary_table_corthick(crosscombat_nonbays_corthick, sshort_corthick)
s2  <- cbind(t2a[,1:4], t2b[,-c(1, 8)], t2f[,-c(1, 2, 8)], t2g[,-c(1, 2, 8)])

#md
t3f <- make_summary_table(crosscombat_nonpara_md, sshort_md)
t3g <- make_summary_table(crosscombat_nonbays_md, sshort_md)
s3  <- cbind(t3a[,1:4], t3b[,-c(1, 8)], t3f[,-c(1, 2, 8)], t3g[,-c(1, 2, 8)])

#fa
t4f <- make_summary_table(crosscombat_nonpara_fa, sshort_fa)
t4g <- make_summary_table(crosscombat_nonbays_fa, sshort_fa)
s4  <- cbind(t4a[,1:4], t4b[,-c(1, 8)], t4f[,-c(1, 2, 8)], t4g[,-c(1, 2, 8)])

s_tall <- rbind(s1, s2, s3, s4)

colnames(s_tall) <- c("ROI" , "Scan_pairs_within", "ICC_within", "Perc_change_within",
                    "Scan_pairs_across",
                    "ICC_unharm", "Perc_change_unharm", "Delta_unharm", "Praw_unharm", "Padj_unharm", "Size_unharm",
                    "ICC_nonpara", "Perc_change_nonpara", "Delta_nonpara","Praw_nonpara", "Padj_nonpara", "Size_nonpara",
                    "ICC_nonbays", "Perc_change_nonbays", "Delta_nonbays","Praw_nonbays", "Padj_nonbays", "Size_nonbays")
#format table for word
ft <- make_flextable_sens(s_tall)

#save table to word
sect_properties <- prop_section(
  page_size = page_size(orient = "landscape"),
  type = "continuous",
  page_margins = page_mar()
)

save_as_docx(ft, path = "Figures/Sensitivity_Cross-sec-CoV_table.docx", pr_section = sect_properties)



# Simulate testing for difference in intercept

mymax = 1000
#simulate volume
fp_unharm_vol      <- sim_unharm(mymax = mymax, vol_for_combat)
fp_crosscombat_vol <- sim_crosscombat(mymax = mymax, vol_for_combat)
fp_longcombat_vol  <- sim_longcombat(mymax = mymax, vol_for_combat)
fp_gamcombat_vol   <- sim_gamcombat(mymax = mymax, vol_for_combat)
fp_vol             <- rbind(fp_unharm_vol, 
                            fp_crosscombat_vol, 
                            fp_longcombat_vol,
                            fp_gamcombat_vol)
write.csv(fp_vol, "Data/fp_vol.csv")

#simulate corthick
fp_unharm_corthick      <- sim_unharm(mymax = mymax, corthick_for_combat, 
                             is_corthick = "yes")
fp_crosscombat_corthick <- sim_crosscombat(mymax = mymax, corthick_for_combat, 
                                  from = "Frontal", to = "WholeCortex", is_corthick = "yes")
fp_longcombat_corthick  <- sim_longcombat(mymax = mymax, corthick_for_combat,
                                 is_corthick = "yes")
fp_gamcombat_corthick   <- sim_gamcombat(mymax = mymax, corthick_for_combat,, 
                                  from = "Frontal", to = "WholeCortex",
                                 is_corthick = "yes")
fp_corthick             <- rbind(fp_unharm_corthick, 
                                 fp_crosscombat_corthick, 
                                 fp_longcombat_corthick,
                                 fp_gamcombat_corthick)
write.csv(fp_corthick, "Data/fp_corthick.csv")

#simulate md
fp_unharm_md      <- sim_unharm(mymax = mymax, md_for_combat)
fp_crosscombat_md <- sim_crosscombat(mymax = mymax, md_for_combat)
fp_longcombat_md  <- sim_longcombat(mymax = mymax, md_for_combat)
fp_gamcombat_md   <- sim_gamcombat(mymax = mymax, md_for_combat)
fp_md             <- rbind(fp_unharm_md, 
                           fp_crosscombat_md, 
                           fp_longcombat_md,
                           fp_gamcombat_md)
write.csv(fp_md, "Data/fp_md.csv")

#simulate fa
fp_unharm_fa      <- sim_unharm(mymax = mymax, fa_for_combat)
fp_crosscombat_fa <- sim_crosscombat(mymax = mymax, fa_for_combat)
fp_longcombat_fa  <- sim_longcombat(mymax = mymax, fa_for_combat)
fp_gamcombat_fa   <- sim_gamcombat(mymax = mymax, fa_for_combat)
fp_fa             <- rbind(fp_unharm_fa, 
                           fp_crosscombat_fa, 
                           fp_longcombat_fa,
                           fp_gamcombat_fa)
write.csv(fp_fa, "Data/fp_fa.csv")
