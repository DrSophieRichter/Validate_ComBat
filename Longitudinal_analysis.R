## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
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
library(forestplot)
library(ggpubr)
library(rstatix)
library(magrittr)
library(effects)
library(sjPlot)
source("Functions/My_table_funs.R")
source("Functions/My_table_funs_corthick.R")
source("Functions/Select_struc_funs.R")
source("Functions/Select_dti_funs.R")
source("Functions/longCombat.R")
source("Functions/Harmonize_funs.R")
source("Functions/Simulate_funs.R")
source("Functions/Followup_funs.R")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
#Load data
scans <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Legacy_data/Most_uptodate_curated_data/Sophies_scan_database_20220121.csv", na.strings = c("", " ", "NA"))
qc1 <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Analysis_Aug2021/Current_data/Check_these_scans_1.csv", na.strings = c("", " ", "NA")) %>%
  select(Scan_ID, Decision)
qc2 <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Analysis_Aug2021/Current_data/Check_these_scans_2.csv", na.strings = c("", " ", "NA")) %>%
  select(Scan_ID, Decision)
qc3 <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Analysis_Aug2021/Current_data/Check_these_scans_3.csv", na.strings = c("", " ", "NA")) %>%
  select(Scan_ID, Decision)
qc4 <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Analysis_Aug2021/CompareCombats/Data/Combat_DTI_QC_filled.csv", na.strings = c("", " ", "NA")) %>%
  select(Scan_ID, Series, Decision = Decision_final)
qc_bad <- rbind(qc1, qc2, qc3) %>% as.data.frame() %>% filter(!Decision == "keep")

vol <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Raw_data/Sophies_processed_data/malpem138_volume_sophie.csv")
mymalp <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Analysis_Aug2021/Current_data/MALPEM_parcellation_Sophie.csv") %>% remove_empty()

corthick <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Raw_data/Sophies_processed_data/malpem021_corthick_sophie.csv", na.strings = c("", "NA"))
add <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Raw_data/Sophies_processed_data/malpem005_corthick_sophie.csv", na.strings = c("", " "))

dti <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Raw_data/Sophies_processed_data/malpem007_dti_calc_Sophie.csv")

clinical <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Legacy_data/Most_uptodate_curated_data/Sophies_clinical_database_20220121.csv", na.strings = c("", " ", "NA"))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
#Prepare structural data
#collapse metrics into desired ROIs
vol         <- collapse_vol(vol, mymalp)
corthick    <- collapse_corthick(corthick, add)

#select structural images for longitudinal cohort
struc_wide <- select_struc_fu(scans, qc_bad, clinical)
struc_long <- wide_to_long_struc(struc_wide, scans, clinical)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
#Now prepare diffusion data
#extract fa or md and DTI series name
md <- extract_dti_fresh(dti,metric = "md") # fresh means I do not require the same scans to have been used for volume
fa <- extract_dti_fresh(dti,metric = "fa")

#select dti images for longitudinal cohort
dti_wide <- select_dti_fu(scans, qc_bad, clinical)
dti_long <- wide_to_long_dti(dti_wide, scans, clinical)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
#Prepare summary table of cohort characteristics
cohort_table_struc <- make_cohort_table_2(struc_wide, struc_long)
cohort_table_struc
cohort_table_dti   <- make_cohort_table_2(dti_wide, dti_long)
cohort_table_dti

## -------------------------------------------------------------------------------------------------------------------------------------------------------------
#check dates of acquisition
temp <- scans %>% filter(Scan_ID %in% struc_long$Scan_ID |
                           Scan_ID %in% dti_long$Scan_ID )
temp <- temp %>% separate(Scan_ID, into = c("wbic", "date", "uid"), sep = "_")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
###############################
#### VOLUME         ###########
###############################
#Create one dataframe for each harmonisation method: unharmonised, crosscombat and longcombat
vol_fu_unharm <- merge(struc_long, vol, by = "Scan_ID", all.x = TRUE, all.y = FALSE)
vol_fu_crossharm <- get_crosscombat(vol_fu_unharm)
vol_fu_crossharm <- cbind(vol_fu_crossharm, vol_fu_unharm[, c("Age_at_this_scan", "Sex", "Label")])
vol_fu_longharm <- get_longcombat(vol_fu_unharm)
vol_fu_longharm <- cbind(vol_fu_longharm, vol_fu_unharm[, c("Age_at_this_scan", "Sex", "Label")])
vol_fu_longharm_plusslope <- get_longcombat_plusslope(vol_fu_unharm)
vol_fu_longharm_plusslope <- cbind(vol_fu_longharm_plusslope, vol_fu_unharm[, c("Age_at_this_scan", "Sex", "Label")])

#bring these tables into the wide format, so that there is one row per subject per roi
#and one column per scan
vol_fu_unharm <- calculate_scanner_induced_difference(vol_fu_unharm, struc_wide)
vol_fu_crossharm <- calculate_scanner_induced_difference(vol_fu_crossharm, struc_wide)
vol_fu_longharm <- calculate_scanner_induced_difference(vol_fu_longharm, struc_wide)
vol_fu_longharm_plusslope <- calculate_scanner_induced_difference(vol_fu_longharm_plusslope, struc_wide)

#Combine the three tables to compare results
vol_fu_unharm$Method <- "Unharmonised"
vol_fu_crossharm$Method <- "CrossCombat"
vol_fu_longharm$Method <- "LongCombat_i"
vol_fu_longharm_plusslope$Method <- "LongCombat_i+s"
vol_fu_all <- rbind(vol_fu_unharm,
                       vol_fu_crossharm,
                       vol_fu_longharm,
                       vol_fu_longharm_plusslope)

#order methods in the order I want them displayed in the graph
vol_fu_all$Method <- factor(vol_fu_all$Method, levels = c("Unharmonised", "CrossCombat", "LongCombat_i", "LongCombat_i+s"))


ggplot(vol_fu_all, aes(x=ROI, CoV, fill=Method)) + 
  geom_boxplot() +ggtitle("Volume - CoV in longitudinal scans") + 
  ylab("Difference in volume (%)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1))



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
###############################
#### CORTICAL THICKNESS #######
###############################
#Create one dataframe for each harmonisation method: unharmonised, crosscombat and longcombat
corthick_fu_unharm <- merge(struc_long, corthick, by = "Scan_ID", all.x = TRUE, all.y = FALSE)
corthick_fu_crossharm <- get_crosscombat(corthick_fu_unharm, from = "Frontal", to = "WholeCortex")
corthick_fu_crossharm <- cbind(corthick_fu_crossharm, corthick_fu_unharm[, c("Age_at_this_scan", "Sex", "Label")])
corthick_fu_longharm <- get_longcombat(corthick_fu_unharm, is_corthick = "yes")
corthick_fu_longharm <- cbind(corthick_fu_longharm, corthick_fu_unharm[, c("Age_at_this_scan", "Sex", "Label")])
corthick_fu_longharm_plusslope <- get_longcombat_plusslope(corthick_fu_unharm, is_corthick = "yes")
corthick_fu_longharm_plusslope <- cbind(corthick_fu_longharm_plusslope, corthick_fu_unharm[, c("Age_at_this_scan", "Sex", "Label")])

#bring these tables into the wide format, so that there is one row per subject per roi
#and one column per scan
corthick_fu_unharm <- calculate_scanner_induced_difference(corthick_fu_unharm, struc_wide, from = "Frontal", to = "WholeCortex")
corthick_fu_crossharm <- calculate_scanner_induced_difference(corthick_fu_crossharm, struc_wide, from = "Frontal", to = "WholeCortex")
corthick_fu_longharm <- calculate_scanner_induced_difference(corthick_fu_longharm, struc_wide, from = "Frontal", to = "WholeCortex")
corthick_fu_longharm_plusslope <- calculate_scanner_induced_difference(corthick_fu_longharm_plusslope, struc_wide, from = "Frontal", to = "WholeCortex")

#Combine the three tables to compare results
corthick_fu_unharm$Method <- "Unharmonised"
corthick_fu_crossharm$Method <- "CrossCombat"
corthick_fu_longharm$Method <- "LongCombat_i"
corthick_fu_longharm_plusslope$Method <- "LongCombat_i+s"
corthick_fu_all <- rbind(corthick_fu_unharm,
                       corthick_fu_crossharm,
                       corthick_fu_longharm,
                       corthick_fu_longharm_plusslope)

#order methods in the order I want them displayed in the graph
corthick_fu_all$Method <- factor(corthick_fu_all$Method, levels = c("Unharmonised", "CrossCombat", "LongCombat_i", "LongCombat_i+s"))

ggplot(corthick_fu_all, aes(x=ROI, CoV, fill=Method)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +ggtitle("Cortical thickness - CoV in longitudinal scans") + 
  ylab("Difference in cortical thickness (%)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
###############################
#### MD             ###########
###############################
#Create one dataframe for each harmonisation method: unharmonised, crosscombat and longcombat
md <- md %>% select(Merge_ID, Ventricles:Brainstem) %>% distinct(Merge_ID, .keep_all = TRUE)
md_fu_unharm <- merge(dti_long, md, by = "Merge_ID", all.x = TRUE, all.y = FALSE)
md_fu_crossharm <- get_crosscombat(md_fu_unharm, is_dti = "yes")
md_fu_crossharm <- cbind(md_fu_crossharm, md_fu_unharm[, c("Age_at_this_scan", "Sex", "Label")])
md_fu_longharm <- get_longcombat(md_fu_unharm, is_dti = "yes")
md_fu_longharm <- cbind(md_fu_longharm, md_fu_unharm[, c("Age_at_this_scan", "Sex", "Label")])
md_fu_longharm_plusslope <- get_longcombat(md_fu_unharm, is_dti = "yes")
md_fu_longharm_plusslope <- cbind(md_fu_longharm_plusslope, md_fu_unharm[, c("Age_at_this_scan", "Sex", "Label")])

#bring these tables into the wide format, so that there is one row per subject per roi
#and one column per scan
md_fu_unharm <- calculate_scanner_induced_difference(md_fu_unharm, struc_wide)
md_fu_crossharm <- calculate_scanner_induced_difference(md_fu_crossharm, struc_wide)
md_fu_longharm <- calculate_scanner_induced_difference(md_fu_longharm, struc_wide)
md_fu_longharm_plusslope <- calculate_scanner_induced_difference(md_fu_longharm_plusslope, struc_wide)

#Combine the three tables to compare results
md_fu_unharm$Method <- "Unharmonised"
md_fu_crossharm$Method <- "CrossCombat"
md_fu_longharm$Method <- "LongCombat_i"
md_fu_longharm_plusslope$Method <- "LongCombat_i+s"
md_fu_all <- rbind(md_fu_unharm,
                       md_fu_crossharm,
                       md_fu_longharm,
                       md_fu_longharm_plusslope)

#order methods in the order I want them displayed in the graph
md_fu_all$Method <- factor(md_fu_all$Method, levels = c("Unharmonised", "CrossCombat", "LongCombat_i", "LongCombat_i+s"))

ggplot(md_fu_all, aes(x=ROI, CoV, fill=Method)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +ggtitle("MD - CoV in longitudinal scans") + 
  ylab("Difference in md (%)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1))

## -------------------------------------------------------------------------------------------------------------------------------------------------------------
###############################
#### FA             ###########
###############################
#Create one dataframe for each harmonisation method: unharmonised, crosscombat and longcombat
fa <- fa %>% select(Merge_ID, Ventricles:Brainstem) %>% distinct(Merge_ID, .keep_all = TRUE)
fa_fu_unharm <- merge(dti_long, fa, by = "Merge_ID", all.x = TRUE, all.y = FALSE)
fa_fu_crossharm <- get_crosscombat(fa_fu_unharm, is_dti = "yes")
fa_fu_crossharm <- cbind(fa_fu_crossharm, fa_fu_unharm[, c("Age_at_this_scan", "Sex", "Label")])
fa_fu_longharm <- get_longcombat(fa_fu_unharm, is_dti = "yes")
fa_fu_longharm <- cbind(fa_fu_longharm, fa_fu_unharm[, c("Age_at_this_scan", "Sex", "Label")])
fa_fu_longharm_plusslope <- get_longcombat_plusslope(fa_fu_unharm, is_dti = "yes")
fa_fu_longharm_plusslope <- cbind(fa_fu_longharm_plusslope, fa_fu_unharm[, c("Age_at_this_scan", "Sex", "Label")])

#bring these tables into the wide format, so that there is one row per subject per roi
#and one column per scan
fa_fu_unharm <- calculate_scanner_induced_difference(fa_fu_unharm, struc_wide)
fa_fu_crossharm <- calculate_scanner_induced_difference(fa_fu_crossharm, struc_wide)
fa_fu_longharm <- calculate_scanner_induced_difference(fa_fu_longharm, struc_wide)
fa_fu_longharm_plusslope <- calculate_scanner_induced_difference(fa_fu_longharm_plusslope, struc_wide)

#Combine the three tables to compare results
fa_fu_unharm$Method <- "Unharmonised"
fa_fu_crossharm$Method <- "CrossCombat"
fa_fu_longharm$Method <- "LongCombat_i"
fa_fu_longharm_plusslope$Method <- "LongCombat_i+s"
fa_fu_all <- rbind(fa_fu_unharm,
                       fa_fu_crossharm,
                       fa_fu_longharm,
                       fa_fu_longharm_plusslope)

#order methods in the order I want them displayed in the graph
fa_fu_all$Method <- factor(fa_fu_all$Method, levels = c("Unharmonised", "CrossCombat", "LongCombat_i", "LongCombat_i+s"))

ggplot(fa_fu_all, aes(x=ROI, CoV, fill=Method)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +ggtitle("FA - CoV in longitudinal scans") + 
  ylab("Difference in FA (%)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
#Add the metric name to each dataframe so that they can be combined
vol_fu_all$Metric       <- "volume"
corthick_fu_all$Metric  <- "corthick"
md_fu_all$Metric        <- "md"
fa_fu_all$Metric        <- "fa"
all_fu                  <- rbind(vol_fu_all,
                                 corthick_fu_all, 
                                 md_fu_all, 
                                 fa_fu_all) %>%
                            as.data.frame()

all_fu$Metric           <- as.factor(all_fu$Metric)



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test if the Cov after harmonization differs significantly from the CoV before harmonization
# by fitting a mixed model for each metric and ROI
# and then store the results in a list called res

data <- all_fu
res <- list()
metrics <- levels(as.factor(data$Metric))
  
for (j in 1:length(metrics)){
    
    metric.data <- 
      data %>% 
      filter(Metric == metrics[j])
    
    rois <- levels(as.factor(metric.data$ROI))
    
    res[[j]] <- as.data.frame(matrix(,ncol=6,nrow=3*length(rois)))
    colnames(res[[j]]) <- c("Method", "Metric", "ROI", "Coef", "LLCI", "ULCI")
    
    for (i in 1:length(rois)){
      
      roi.data <- 
        metric.data %>%
        filter(ROI == rois[i]) 
      
      fit <- lmer(CoV ~ Method +  (1| Master_subject_ID), data = roi.data)
      sum <- summary(fit)
      
      res[[j]][(i-1)*3+1, "Coef"]   <- sum$coefficients["MethodCrossCombat","Estimate"]
      res[[j]][(i-1)*3+2, "Coef"]   <- sum$coefficients["MethodLongCombat_i","Estimate"]
      res[[j]][(i-1)*3+3, "Coef"]   <- sum$coefficients["MethodLongCombat_i+s","Estimate"]
      res[[j]][c((i-1)*3+1,(i-1)*3+2, (i-1)*3+3), "ROI"]    <- rois[i]
      res[[j]][c((i-1)*3+1,(i-1)*3+2, (i-1)*3+3), "Metric"] <- metrics[j]
      res[[j]][(i-1)*3+1, "Method"] <- "neuroCombat"
      res[[j]][(i-1)*3+2, "Method"] <- "longCombat_i"
      res[[j]][(i-1)*3+3, "Method"] <- "longCombat_i+s"
      res[[j]][(i-1)*3+1, "LLCI"]   <- confint(fit)["MethodCrossCombat","2.5 %"]
      res[[j]][(i-1)*3+2, "LLCI"]   <- confint(fit)["MethodLongCombat_i","2.5 %"]
      res[[j]][(i-1)*3+3, "LLCI"]   <- confint(fit)["MethodLongCombat_i+s","2.5 %"]
      res[[j]][(i-1)*3+1, "ULCI"]   <- confint(fit)["MethodCrossCombat","97.5 %"]
      res[[j]][(i-1)*3+2, "ULCI"]   <- confint(fit)["MethodLongCombat_i","97.5 %"]
      res[[j]][(i-1)*3+3, "ULCI"]   <- confint(fit)["MethodLongCombat_i+s","97.5 %"]
    }
    
  }

#combine all results from the list into a single table
res_methods <- bind_rows(res)

#I want to illustrate these effects in a forest plot

#Re-order metrics prior to plotting
res_methods$Metric <- factor(res_methods$Metric,
                         levels = c("volume", "corthick", "md", "fa"),
                         labels = c("Volume", "Cortical thickness", "Mean diffusivity", "Fractional anisotropy"))

data <- res_methods
data$ROI <-  factor(data$ROI)
  data$Method <- factor(data$Method, 
                        levels = c("longCombat_i+s", "longCombat_i", "neuroCombat"),
                        labels = c("longCombat i+s", "longCombat i", "neuroCombat"))
  
#define colours for dots and bars
dotCOLS = c("hotpink4", "lightskyblue3", "steelblue3")
barCOLS = c("hotpink3", "lightskyblue1", "steelblue4")
  
  
p <- ggplot(data, aes(x=ROI, y=Coef, ymin=LLCI, ymax=ULCI, col=Method, fill=Method)) + 
     geom_hline(yintercept=0, size = 1, colour = "grey80") +
     #specify position here
     geom_linerange(size=1.8,position=position_dodge(width = 0.8)) +
     #specify position here too
     geom_point(size=1.8, shape=21, colour="white", stroke = 0.5 ,position=position_dodge(width = 0.8)) +
     scale_fill_manual(values=dotCOLS, breaks = c("neuroCombat", "longCombat i", "longCombat i+s"))+
     scale_color_manual(values=barCOLS, breaks = c("neuroCombat", "longCombat i", "longCombat i+s"))+
     scale_x_discrete(name="Regions of interest") +
     scale_y_continuous(name="Scanner effect relative to unharmonised data",  
                       limits = c(-18, 18),
                       breaks = seq(-18, 18, by = 2)) +
     coord_flip() +
     theme_bw() +
     facet_grid(rows = vars(Metric), space = "free", scales = "free") +
     theme(legend.position = "bottom",
          strip.background = element_rect(fill="grey95"))
 
ggsave(plot = p, 
       filename = "Figures/forrest_plot_CoV.tiff",
       width = 16, height = 24, dpi = 300,
       units = "cm")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Assessing false positive rate

# Get dataframes of unharmonised data in long format
vol_for_combat <- merge(struc_long, vol, by = "Scan_ID", all.x = TRUE, all.y = FALSE)
corthick_for_combat <- merge(struc_long, corthick, by = "Scan_ID", all.x = TRUE, all.y = FALSE)
md_for_combat <- merge(dti_long, md, by = "Merge_ID", all.x = TRUE, all.y = FALSE)
fa_for_combat <- merge(dti_long, fa, by = "Merge_ID", all.x = TRUE, all.y = FALSE)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Simulate and test for difference in INTERCEPT between groups

mymax = 1000
#simulate volume
fp_unharm_vol      <- sim_unharm(mymax = mymax, vol_for_combat)
fp_crosscombat_vol <- sim_crosscombat(mymax = mymax, vol_for_combat)
fp_longcombat_vol  <- sim_longcombat(mymax = mymax, vol_for_combat)
fp_vol <- rbind(fp_unharm_vol, fp_crosscombat_vol, fp_longcombat_vol)
write.csv(fp_vol, "Data/fu_fp_vol.csv")

#simulate corthick
fp_unharm_corthick      <- sim_unharm(mymax = mymax, corthick_for_combat, 
                             is_corthick = "yes")
fp_crosscombat_corthick <- sim_crosscombat(mymax = mymax, corthick_for_combat, 
                                  from = "Frontal", to = "WholeCortex", is_corthick = "yes")
fp_longcombat_corthick  <- sim_longcombat(mymax = mymax, corthick_for_combat,
                                 is_corthick = "yes")
fp_corthick    <- rbind(fp_unharm_corthick, fp_crosscombat_corthick, fp_longcombat_corthick)
write.csv(fp_corthick, "Data/fu_fp_corthick.csv")

#simulate md
fp_unharm_md      <- sim_unharm(mymax = mymax, md_for_combat)
fp_crosscombat_md <- sim_crosscombat(mymax = mymax, md_for_combat, is_fudti = "yes")
fp_longcombat_md  <- sim_longcombat(mymax = mymax, md_for_combat)
fp_md             <- rbind(fp_unharm_md, fp_crosscombat_md, fp_longcombat_md)
write.csv(fp_md, "Data/fu_fp_md.csv")

#simulate fa
fp_unharm_fa      <- sim_unharm(mymax = mymax, fa_for_combat)
fp_crosscombat_fa <- sim_crosscombat(mymax = mymax, fa_for_combat, is_fudti = "yes")
fp_longcombat_fa  <- sim_longcombat(mymax = mymax, fa_for_combat)
fp_fa             <- rbind(fp_unharm_fa, fp_crosscombat_fa, fp_longcombat_fa)
write.csv(fp_fa, "Data/fu_fp_fa.csv")


# Simulate and test for difference in SLOPE between groups

mymax = 1000
#simulate volume
fp_unharm_vol      <- sim_unharm(mymax = mymax, vol_for_combat, intercept_only = "no")
fp_crosscombat_vol <- sim_crosscombat(mymax = mymax, vol_for_combat, intercept_only = "no")
fp_longcombat_vol  <- sim_longcombat(mymax = mymax, vol_for_combat, intercept_only = "no")
fp_vol <- rbind(fp_unharm_vol, fp_crosscombat_vol, fp_longcombat_vol)
write.csv(fp_vol, "Data/fu_fp_vol_slope.csv")

#simulate corthick
fp_unharm_corthick      <- sim_unharm(mymax = mymax, corthick_for_combat,
                             is_corthick = "yes", intercept_only = "no")
fp_crosscombat_corthick <- sim_crosscombat(mymax = mymax, corthick_for_combat, 
                                  from = "Frontal", to = "WholeCortex", is_corthick = "yes", , intercept_only = "no")
fp_longcombat_corthick  <- sim_longcombat(mymax = mymax, corthick_for_combat,
                                 is_corthick = "yes", intercept_only = "no")
fp_corthick    <- rbind(fp_unharm_corthick, fp_crosscombat_corthick, fp_longcombat_corthick)
write.csv(fp_corthick, "Data/fu_fp_corthick_slope.csv")

#simulate md
fp_unharm_md      <- sim_unharm(mymax = mymax, md_for_combat, intercept_only = "no")
fp_crosscombat_md <- sim_crosscombat(mymax = mymax, md_for_combat, is_fudti = "yes", intercept_only = "no")
fp_longcombat_md  <- sim_longcombat(mymax = mymax, md_for_combat, intercept_only = "no")
fp_md             <- rbind(fp_unharm_md, fp_crosscombat_md, fp_longcombat_md)
write.csv(fp_md, "Data/fu_fp_md_slope.csv")

#simulate fa
fp_unharm_fa      <- sim_unharm(mymax = mymax, fa_for_combat, intercept_only = "no")
fp_crosscombat_fa <- sim_crosscombat(mymax = mymax, fa_for_combat, is_fudti = "yes", intercept_only = "no")
fp_longcombat_fa  <- sim_longcombat(mymax = mymax, fa_for_combat, intercept_only = "no")
fp_fa             <- rbind(fp_unharm_fa, fp_crosscombat_fa, fp_longcombat_fa)
write.csv(fp_fa, "Data/fu_fp_fa_slope.csv")


