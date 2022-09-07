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
library(forestplot)
library(ggpubr)
library(rstatix)
library(magrittr)
library(effects)
library(sjPlot)
library(flextable)
library(grid)
library(gridExtra)
library(reticulate)
source("Functions/My_table_funs.R")
source("Functions/My_table_funs_corthick.R")
source("Functions/Select_struc_funs.R")
source("Functions/Select_dti_funs.R")
source("Functions/longCombat.R")
source("Functions/Harmonize_funs.R")
source("Functions/Simulate_funs.R")
source("Functions/Followup_funs.R")
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
qc4 <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Analysis_Aug2021/CompareCombats/Data/Combat_DTI_QC_filled.csv", na.strings = c("", " ", "NA")) %>%
  select(Scan_ID, Series, Decision = Decision_final)
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

#select structural images for longitudinal cohort
struc_wide <- select_struc_fu(scans, qc_bad, clinical)
struc_long <- wide_to_long_struc(struc_wide, scans, clinical)

#Now prepare diffusion data
#extract fa or md and DTI series name
md <- extract_dti_fresh(dti,metric = "md") # fresh means I do not require the same scans to have been used for volume
fa <- extract_dti_fresh(dti,metric = "fa")

#select dti images for longitudinal cohort
dti_wide <- select_dti_fu(scans, qc_bad, clinical)
dti_long <- wide_to_long_dti(dti_wide, scans, clinical)

#Prepare summary table of cohort characteristics
cohort_table_struc <- make_cohort_table_2(struc_wide, struc_long)
cohort_table_struc
cohort_table_dti   <- make_cohort_table_2(dti_wide, dti_long)
cohort_table_dti
#check dates of acquisition
temp <- scans %>% filter(Scan_ID %in% struc_long$Scan_ID |
                           Scan_ID %in% dti_long$Scan_ID )
temp <- temp %>% separate(Scan_ID, into = c("wbic", "date", "uid"), sep = "_")

###############################
#### VOLUME         ###########
###############################
#Create one dataframe for unharmonised data
vol_fu_unharm             <- 
  merge(struc_long, vol, by = "Scan_ID", all.x = TRUE, all.y = FALSE) 

#create identical dataframes with differently harmnonised data
vol_fu_crossharm          <- 
  get_crosscombat(vol_fu_unharm) %>% #harmonize data
  cbind(., vol_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>% #covariates
  calculate_scanner_induced_difference(., struc_wide) %>% #bring into wide format and calculate differences btw scans
  tibble::add_column(Method = "CrossCombat") #indicate harmonisation method

vol_fu_crossharm_nonpara  <- 
  get_crosscombat(vol_fu_unharm, prior = "non-parametric") %>% #harmonize data
  cbind(., vol_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>% #covariates
  calculate_scanner_induced_difference(., struc_wide) %>% #bring into wide format and calculate differences btw scans
  tibble::add_column(Method = "CrossCombat_nonpara") #indicate harmonisation method

vol_fu_crossharm_nonbays  <- 
  get_crosscombat(vol_fu_unharm, prior = "none") %>% #harmonize data
  cbind(., vol_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>% #covariates
  calculate_scanner_induced_difference(., struc_wide) %>% #bring into wide format and calculate differences btw scans
  tibble::add_column(Method = "CrossCombat_nonbays") #indicate harmonisation method

vol_fu_longharm           <- 
  get_longcombat(vol_fu_unharm) %>% #harmonize data
  cbind(., vol_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>% #covariates
  calculate_scanner_induced_difference(., struc_wide) %>% #bring into wide format and calculate differences btw scans
  tibble::add_column(Method = "LongCombat_i") #indicate harmonisation method

vol_fu_longharm_plusslope <- 
  get_longcombat_plusslope(vol_fu_unharm) %>% #harmonize data
  cbind(., vol_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>% #covariates
  calculate_scanner_induced_difference(., struc_wide) %>% #bring into wide format and calculate differences btw scans
  tibble::add_column(Method = "LongCombat_i+s") #indicate harmonisation method

vol_fu_gamharm          <- 
  get_gamcombat(vol_fu_unharm) %>% #harmonize data
  cbind(., vol_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>% #covariates
  calculate_scanner_induced_difference(., struc_wide) %>% #bring into wide format and calculate differences btw scans
  tibble::add_column(Method = "GamCombat") #indicate harmonisation method

vol_fu_unharm <-
  vol_fu_unharm %>%
  calculate_scanner_induced_difference(., struc_wide)  %>% #bring into wide format and calculate differences btw scans
  tibble::add_column(Method = "Unharmonised") #indicate harmonisation method

#Combine the three tables to compare results
vol_fu_all <- rbind(vol_fu_unharm,
                    vol_fu_crossharm,
                    vol_fu_crossharm_nonpara,
                    vol_fu_crossharm_nonbays,
                    vol_fu_longharm,
                    vol_fu_longharm_plusslope,
                    vol_fu_gamharm)

#order methods in the order I want them displayed in the graph
vol_fu_all$Method <- 
  factor(vol_fu_all$Method, 
         levels = c("Unharmonised", 
                    "CrossCombat", 
                    "CrossCombat_nonpara",
                    "CrossCombat_nonbays",
                    "LongCombat_i", 
                    "LongCombat_i+s",
                    "GamCombat"))

#colorblind friendly palette
cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

vol_plot <-
  ggplot(vol_fu_all, aes(x=ROI, CoV, fill=Method)) + 
  geom_boxplot() +ggtitle("Volume - CoV in longitudinal scans") + 
  ylab("Coefficient of variation (%)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1),
        legend.position = "none") +
  scale_fill_manual(values=cbPalette)
vol_plot

###############################
#### CORTICAL THICKNESS #######
###############################
#Create one dataframe for unharmonised data
corthick_fu_unharm             <- 
  merge(struc_long, corthick, by = "Scan_ID", all.x = TRUE, all.y = FALSE) 

#create identical dataframes with differently harmnonised data
corthick_fu_crossharm          <- 
  get_crosscombat(corthick_fu_unharm, from = "Frontal", to = "WholeCortex") %>% #harmonize data
  cbind(., corthick_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>% #covariates
  calculate_scanner_induced_difference(., 
                                       struc_wide, 
                                       from = "Frontal", to = "WholeCortex") %>% #bring into wide format and calculate differences btw scans
  tibble::add_column(Method = "CrossCombat") #indicate harmonisation method

corthick_fu_crossharm_nonpara  <- 
  get_crosscombat(corthick_fu_unharm, prior = "non-parametric", from = "Frontal", to = "WholeCortex") %>% #harmonize data
  cbind(., corthick_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>% #covariates
  calculate_scanner_induced_difference(., 
                                       struc_wide, 
                                       from = "Frontal", to = "WholeCortex") %>% #bring into wide format and calculate differences btw scans
  tibble::add_column(Method = "CrossCombat_nonpara") #indicate harmonisation method

corthick_fu_crossharm_nonbays  <- 
  get_crosscombat(corthick_fu_unharm, prior = "none", from = "Frontal", to = "WholeCortex") %>% #harmonize data
  cbind(., corthick_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>% #covariates
  calculate_scanner_induced_difference(., 
                                       struc_wide, 
                                       from = "Frontal", to = "WholeCortex") %>% #bring into wide format and calculate differences btw scans
  tibble::add_column(Method = "CrossCombat_nonbays") #indicate harmonisation method

corthick_fu_longharm           <- 
  get_longcombat(corthick_fu_unharm, is_corthick = "yes") %>% #harmonize data
  cbind(., corthick_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>% #covariates
  calculate_scanner_induced_difference(., 
                                       struc_wide, 
                                       from = "Frontal", to = "WholeCortex") %>% #bring into wide format and calculate differences btw scans
  tibble::add_column(Method = "LongCombat_i") #indicate harmonisation method

corthick_fu_longharm_plusslope <- 
  get_longcombat_plusslope(corthick_fu_unharm, is_corthick = "yes") %>% #harmonize data
  cbind(., corthick_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>% #covariates
  calculate_scanner_induced_difference(., 
                                       struc_wide, 
                                       from = "Frontal", to = "WholeCortex") %>% #bring into wide format and calculate differences btw scans
  tibble::add_column(Method = "LongCombat_i+s") #indicate harmonisation method

corthick_fu_gamharm          <- 
  get_gamcombat(corthick_fu_unharm, from = "Frontal", to = "WholeCortex") %>% #harmonize data
  cbind(., corthick_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>% #covariates
  calculate_scanner_induced_difference(., 
                                       struc_wide, 
                                       from = "Frontal", to = "WholeCortex") %>% #bring into wide format and calculate differences btw scans
  tibble::add_column(Method = "GamCombat") #indicate harmonisation method


corthick_fu_unharm <-
  corthick_fu_unharm %>%
  calculate_scanner_induced_difference(., 
                                       struc_wide, 
                                       from = "Frontal", to = "WholeCortex")  %>% #bring into wide format and calculate differences btw scans
  tibble::add_column(Method = "Unharmonised") #indicate harmonisation method

#Combine the three tables to compare results
corthick_fu_all <- rbind(corthick_fu_unharm,
                    corthick_fu_crossharm,
                    corthick_fu_crossharm_nonpara,
                    corthick_fu_crossharm_nonbays,
                    corthick_fu_longharm,
                    corthick_fu_longharm_plusslope,
                    corthick_fu_gamharm)

#order methods in the order I want them displayed in the graph
corthick_fu_all$Method <- 
  factor(corthick_fu_all$Method, 
         levels = c("Unharmonised", 
                    "CrossCombat", 
                    "CrossCombat_nonpara",
                    "CrossCombat_nonbays",
                    "LongCombat_i", 
                    "LongCombat_i+s",
                    "GamCombat"))


corthick_plot <-
  ggplot(corthick_fu_all, aes(x=ROI, CoV, fill=Method)) + 
  geom_boxplot() +ggtitle("Cortical thickness - CoV in longitudinal scans") + 
  ylab("Coefficient of variation (%)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1),
        legend.position = "none") +
  scale_fill_manual(values=cbPalette)

corthick_plot

###############################
#### MD             ###########
###############################
#Create one dataframe for unharmonised data
md <- 
  md %>% 
  select(Merge_ID, Ventricles:Brainstem) %>% 
  distinct(Merge_ID, .keep_all = TRUE)
md_fu_unharm <- 
  merge(dti_long, md, by = "Merge_ID", all.x = TRUE, all.y = FALSE)

#create identical dataframes with differently harmnonised data
md_fu_crossharm <- 
  get_crosscombat(md_fu_unharm, is_dti = "yes") %>%
   cbind(., md_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>%
  calculate_scanner_induced_difference(., dti_wide) %>%
  tibble::add_column(Method = "CrossCombat")

md_fu_crossharm_nonpara <- 
  get_crosscombat(md_fu_unharm, is_dti = "yes", prior = "non-parametric") %>%
  cbind(., md_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>%
  calculate_scanner_induced_difference(., dti_wide) %>%
  tibble::add_column(Method = "CrossCombat_nonpara")

md_fu_crossharm_nonbays <- 
  get_crosscombat(md_fu_unharm, is_dti = "yes", prior = "none") %>%
  cbind(., md_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>%
  calculate_scanner_induced_difference(., dti_wide) %>%
  tibble::add_column(Method = "CrossCombat_nonbays")
  
md_fu_longharm <- 
  get_longcombat(md_fu_unharm, is_dti = "yes")%>%
  cbind(., md_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>%
  calculate_scanner_induced_difference(., dti_wide) %>%
  tibble::add_column(Method = "LongCombat_i")

md_fu_longharm_plusslope <- 
  get_longcombat(md_fu_unharm, is_dti = "yes")%>%
  cbind(., md_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>%
  calculate_scanner_induced_difference(., dti_wide) %>%
  tibble::add_column(Method = "LongCombat_i+s")

md_fu_gamharm <- 
  get_gamcombat(md_fu_unharm, is_dti = "yes") %>%
   cbind(., md_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>%
  calculate_scanner_induced_difference(., dti_wide) %>%
  tibble::add_column(Method = "GamCombat")

md_fu_unharm <-
  md_fu_unharm %>%
  calculate_scanner_induced_difference(., dti_wide)%>%
  tibble::add_column(Method = "Unharmonised")

#Combine the three tables to compare results
md_fu_all <- rbind(md_fu_unharm,
                   md_fu_crossharm,
                   md_fu_crossharm_nonpara,
                   md_fu_crossharm_nonbays,
                   md_fu_longharm,
                   md_fu_longharm_plusslope,
                   md_fu_gamharm)

#order methods in the order I want them displayed in the graph
md_fu_all$Method <- 
  factor(md_fu_all$Method, 
         levels = c("Unharmonised", 
                    "CrossCombat", 
                    "CrossCombat_nonpara",
                    "CrossCombat_nonbays",
                    "LongCombat_i", 
                    "LongCombat_i+s",
                    "GamCombat"))

md_plot <- 
  ggplot(md_fu_all, aes(x=ROI, CoV, fill=Method)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +ggtitle("MD - CoV in longitudinal scans") + 
  ylab("Coefficient of variation (%)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1),
        legend.position = "none") +
  scale_fill_manual(values=cbPalette)

md_plot

###############################
#### FA             ###########
###############################
#Create one dataframe for unharmonised data
fa <- 
  fa %>% 
  select(Merge_ID, Ventricles:Brainstem) %>% 
  distinct(Merge_ID, .keep_all = TRUE)
fa_fu_unharm <- 
  merge(dti_long, fa, by = "Merge_ID", all.x = TRUE, all.y = FALSE)

#create identical dataframes with differently harmnonised data
fa_fu_crossharm <- 
  get_crosscombat(fa_fu_unharm, is_dti = "yes") %>%
   cbind(., fa_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>%
  calculate_scanner_induced_difference(., dti_wide) %>%
  tibble::add_column(Method = "CrossCombat")

fa_fu_crossharm_nonpara <- 
  get_crosscombat(fa_fu_unharm, is_dti = "yes", prior = "non-parametric") %>%
  cbind(., fa_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>%
  calculate_scanner_induced_difference(., dti_wide) %>%
  tibble::add_column(Method = "CrossCombat_nonpara")

fa_fu_crossharm_nonbays <- 
  get_crosscombat(fa_fu_unharm, is_dti = "yes", prior = "none") %>%
  cbind(., fa_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>%
  calculate_scanner_induced_difference(., dti_wide) %>%
  tibble::add_column(Method = "CrossCombat_nonbays")
  
fa_fu_longharm <- 
  get_longcombat(fa_fu_unharm, is_dti = "yes")%>%
  cbind(., fa_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>%
  calculate_scanner_induced_difference(., dti_wide) %>%
  tibble::add_column(Method = "LongCombat_i")

fa_fu_longharm_plusslope <- 
  get_longcombat(fa_fu_unharm, is_dti = "yes")%>%
  cbind(., fa_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>%
  calculate_scanner_induced_difference(., dti_wide) %>%
  tibble::add_column(Method = "LongCombat_i+s")

fa_fu_gamharm <- 
  get_gamcombat(fa_fu_unharm, is_dti = "yes") %>%
   cbind(., fa_fu_unharm[, c("Age_at_this_scan", "Sex","ICV","Label")]) %>%
  calculate_scanner_induced_difference(., dti_wide) %>%
  tibble::add_column(Method = "GamCombat")

fa_fu_unharm <-
  fa_fu_unharm %>%
  calculate_scanner_induced_difference(., dti_wide)%>%
  tibble::add_column(Method = "Unharmonised")

#Combine the three tables to compare results
fa_fu_all <- rbind(fa_fu_unharm,
                   fa_fu_crossharm,
                   fa_fu_crossharm_nonpara,
                   fa_fu_crossharm_nonbays,
                   fa_fu_longharm,
                   fa_fu_longharm_plusslope,
                   fa_fu_gamharm)

#order methods in the order I want them displayed in the graph
fa_fu_all$Method <- 
  factor(fa_fu_all$Method, 
         levels = c("Unharmonised", 
                    "CrossCombat", 
                    "CrossCombat_nonpara",
                    "CrossCombat_nonbays",
                    "LongCombat_i", 
                    "LongCombat_i+s",
                    "GamCombat"))

fa_plot <- 
  ggplot(fa_fu_all, aes(x=ROI, CoV, fill=Method)) + 
  geom_boxplot() +ggtitle("FA - CoV in longitudinal scans") + 
  ylab("Coefficient of variation (%)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1),
        legend.position = "bottom") +
  scale_fill_manual(values=cbPalette)

fa_plot
get_legend  <- function(myggplot)
  {
  tmp       <- ggplot_gtable(ggplot_build(myggplot))
  leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend    <- tmp$grobs[[leg]]
  return(legend)
  }

legend      <- get_legend(fa_plot)

fa_plot     <- fa_plot + theme(legend.position="none")

myplot <- 
  grid.arrange(
  vol_plot, 
  corthick_plot, 
  md_plot,
  fa_plot,
  legend, 
  ncol = 1,  
  layout_matrix = rbind(1, 1, 1, 1, 1, 1, 1, 1,
                        2, 2, 2, 2, 2, 2, 2, 2,
                        3, 3, 3, 3, 3, 3, 3, 3,
                        4, 4, 4, 4, 4, 4, 4, 4,
                        5))

ggsave(plot = myplot, 
       filename = "Figures/CoV_plot.tiff",
       width = 16, height = 28, dpi = 300,
       units = "cm")

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

all_fu$Metric           <- factor(all_fu$Metric,
                                  levels = c("volume", "corthick", "md", "fa"))

levels(all_fu$Method)


#Calculate the ICC for Rate_within vs Rate_across
#for each method and each metrics and each roi
#as a measure of true positives after harmonisation
metrics <- levels(all_fu$Metric)
methods <- levels(all_fu$Method)
cols <- c("Metric", "Method", "ROI", "Rate_within", "Rate_across", "ICCmean", "ICCsd")
plusminus <-"\u00b1" 
Encoding(plusminus)<-"UTF-8"

icc_data <- 
  all_fu %>%
  ungroup() %>%
  group_by(Method) %>%
  group_split()
names(icc_data) <- levels(all_fu$Method)

res_list<-lapply(my.list<-vector(mode = 'list',length(methods)),
                 function(x) x<-vector(mode='list',length(metrics)))
names(res_list) <- methods

for (i in 1:length(methods)){
  for (j in 1:length(metrics)){
    sub2 <-
  icc_data[[i]] %>%
  filter(Metric == metrics[j])

  rois <- rev(levels(as.factor(sub2$ROI)))

  res <- as.data.frame(matrix(,
                            ncol = length(cols),
                            nrow = 7))
  colnames(res) <- cols
  
for (k in 1:length(rois)){
  
  sub3 <-
    sub2 %>%
    filter(ROI == rois[k])
  
  res[k, "ROI"]    <- rois[k]
  res[k, "Method"] <- methods[i]
  res[k, "Metric"] <- metrics[j]
  res[k, "Rate_within"] <- mean(sub3$Rate_within)*1000 #converting cubic cm/year to cubic mm/year
  res[k, "Rate_across"] <- mean(sub3$Rate_across)*1000 #converting cubic cm/year to cubic mm/year
  res[k, "Rate_within"] <- format(round(as.numeric(res[k, "Rate_within"]), 2), nsmall = 2)
  res[k, "Rate_across"] <- format(round(as.numeric(res[k, "Rate_across"]), 2), nsmall = 2)
  
  #I can calculate ICC for rates
  #but I cannot calculate CoV for rates as
  #the result would be affected by the units chosen
  tab  <- psych::ICC(sub3[,c("Rate_within", "Rate_across")])
  icc  <- tab$results["Average_raters_absolute", "ICC"]
  icc  <- format(round(as.numeric(icc), 2), nsmall = 2)
  sd   <- get_icc_sd(tab)
  sd   <- format(round(as.numeric(sd), 2), nsmall = 2)
  res[k, "ICCmean"] <- icc
  res[k, "ICCsd"]   <-sd
  
}
  res_list[[i]][[j]] <- res
  }
}


all_icc <- bind_cols(bind_rows(res_list[[1]])[c("Metric", "ROI","ICCmean", "ICCsd")], 
                     bind_rows(res_list[[2]])[c("ICCmean", "ICCsd")], 
                     bind_rows(res_list[[3]])[c("ICCmean", "ICCsd")], 
                     bind_rows(res_list[[4]])[c("ICCmean", "ICCsd")],
                     bind_rows(res_list[[5]])[c("ICCmean", "ICCsd")], 
                     bind_rows(res_list[[6]])[c("ICCmean", "ICCsd")],
                     bind_rows(res_list[[7]])[c("ICCmean", "ICCsd")])

colnames(all_icc) <- c("Metric", 
                       "ROI", 
                       paste(levels(all_fu$Method)[1], "mean"),
                       paste(levels(all_fu$Method)[1], "sd"),
                       paste(levels(all_fu$Method)[2], "mean"),
                       paste(levels(all_fu$Method)[2], "sd"),
                       paste(levels(all_fu$Method)[3], "mean"),
                       paste(levels(all_fu$Method)[3], "sd"),
                       paste(levels(all_fu$Method)[4], "mean"),
                       paste(levels(all_fu$Method)[4], "sd"),
                       paste(levels(all_fu$Method)[5], "mean"),
                       paste(levels(all_fu$Method)[5], "sd"),
                       paste(levels(all_fu$Method)[6], "mean"),
                       paste(levels(all_fu$Method)[6], "sd"),
                       paste(levels(all_fu$Method)[7], "mean"),
                       paste(levels(all_fu$Method)[7], "sd"))

colnames(all_icc)

#True positives - check correlation between the ground truth rate and the across-scanner rate

colourer <- scales::col_numeric(
  palette = "PiYG",
  domain = c(0, 1))


all_icc[,-c(1:2)] <- 
  all_icc[,-c(1:2)] %>%
  mutate_all(as.numeric)



ft <- 
    all_icc %>% 
    regulartable() %>%
    set_header_labels("Unharmonised mean"        = "Mean",
                      "CrossCombat mean"         = "Mean",
                      "CrossCombat_nonpara mean" = "Mean",
                      "CrossCombat_nonbays mean" = "Mean",
                      "LongCombat_i mean"        = "Mean",
                      "LongCombat_i+s mean"      = "Mean",
                      "GamCombat mean"           = "Mean",
                      "Unharmonised sd"        = "SD",
                      "CrossCombat sd"         = "SD",
                      "CrossCombat_nonpara sd" = "SD",
                      "CrossCombat_nonbays sd" = "SD",
                      "LongCombat_i sd"        = "SD",
                      "LongCombat_i+s sd"      = "SD",
                      "GamCombat sd"           = "SD",
                      "ROI"                      = "Region of interest") %>%
  bg(
    bg = colourer,
    j = c("Unharmonised mean", 
          "CrossCombat mean", 
          "CrossCombat_nonpara mean", 
          "CrossCombat_nonbays mean", 
          "LongCombat_i mean", 
          "LongCombat_i+s mean", 
          "GamCombat mean"),
    part = "body") %>%
    autofit() %>%
    add_header_row(values = c(" ", "Un-\nharmonized", "neuro\nCombat\n(para)", "neuro\nCombat\n(non-para)", "neuro\nCombat\n(non-bays)", "long\nCombat\n(i only)", "long\nCombat\n(i+s)", "gam\nCombat\n(age)"),
                   colwidths = c(2, 2, 2, 2, 2, 2, 2, 2))%>%
    add_header_row(values = c(" ", "Preserved true biological effect\nmeasured as ICC with ground truth"),
                   colwidths = c(2, 14))%>%
    vline_left() %>%
    vline_right() %>%
    fontsize(size = 8, part = "all") %>%
    set_table_properties(layout = "autofit")%>% 
    padding(padding = 0, part = "all", padding.right = 1.5, padding.left = 1.5) %>%
    align(align = "right", part = "all")%>%
    align(align = "left", j = c(1), part = "all") %>%
    align(align = "center", part = "header", i = c(1, 2)) %>%
  vline(j = c(2, 4, 6, 8, 10, 12, 14, 16), border = fp_border_default()) %>%
  hline(i = c(7, 14, 21, 28), border = fp_border_default(width = 2)) %>%
  rotate(j = 1, rotation = "lrtb", align = "center",part = "all") %>%
  merge_v(j = ~Metric)
  

save_as_docx(ft, path = "Figures/True_positive_table.docx")
levels(all_fu$Method)

# Test if the Cov after harmonization differs significantly from the CoV before harmonization
# by fitting a mixed model for each metric and ROI
# and then store the results in a list called res

data <- 
  all_fu %>% 
  filter(Method %in% c("Unharmonised", "CrossCombat", "LongCombat_i", "GamCombat"))
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
      res[[j]][(i-1)*3+3, "Coef"]   <- sum$coefficients["MethodGamCombat","Estimate"]
      res[[j]][c((i-1)*3+1,(i-1)*3+2, (i-1)*3+3), "ROI"]    <- rois[i]
      res[[j]][c((i-1)*3+1,(i-1)*3+2, (i-1)*3+3), "Metric"] <- metrics[j]
      res[[j]][(i-1)*3+1, "Method"] <- "neuroCombat"
      res[[j]][(i-1)*3+2, "Method"] <- "longCombat"
      res[[j]][(i-1)*3+3, "Method"] <- "gamCombat"
      res[[j]][(i-1)*3+1, "LLCI"]   <- confint(fit)["MethodCrossCombat","2.5 %"]
      res[[j]][(i-1)*3+2, "LLCI"]   <- confint(fit)["MethodLongCombat_i","2.5 %"]
      res[[j]][(i-1)*3+3, "LLCI"]   <- confint(fit)["MethodGamCombat","2.5 %"]
      res[[j]][(i-1)*3+1, "ULCI"]   <- confint(fit)["MethodCrossCombat","97.5 %"]
      res[[j]][(i-1)*3+2, "ULCI"]   <- confint(fit)["MethodLongCombat_i","97.5 %"]
      res[[j]][(i-1)*3+3, "ULCI"]   <- confint(fit)["MethodGamCombat","97.5 %"]
    }
    
  }

#combine all results from the list into a single table
res_methods <- bind_rows(res)

#I want to illustrate these effects in a forest plot

#Re-order metrics prior to plotting
res_methods$Metric <- factor(res_methods$Metric,
                         levels = c("volume", "corthick", "md", "fa"),
                         labels = c("Volume", "Cortical thickness", "Mean diffusivity", "Fractional anisotropy"))

data        <- res_methods
data$ROI    <- factor(data$ROI)
data$Method <- factor(data$Method, 
                        levels = c("gamCombat", "longCombat", "neuroCombat"),
                        labels = c("gamCombat", "longCombat", "neuroCombat"))
  
#define colours for dots and bars
dotCOLS = c("darkgreen", "burlywood4", "cadetblue4")
barCOLS = c("chartreuse4", "burlywood4", "cadetblue3")
  
  
p <- ggplot(data, aes(x=ROI, y=Coef, ymin=LLCI, ymax=ULCI, col=Method, fill=Method)) + 
     geom_hline(yintercept=0, size = 1, colour = "grey80") +
     #specify position here
     geom_linerange(size=1.8,position=position_dodge(width = 0.8)) +
     #specify position here too
     geom_point(size=1.8, shape=21, colour="white", stroke = 0.5 ,position=position_dodge(width = 0.8)) +
     scale_fill_manual(values=dotCOLS, breaks = c("neuroCombat", "longCombat", "gamCombat"))+
     scale_color_manual(values=barCOLS, breaks = c("neuroCombat", "longCombat", "gamCombat"))+
     scale_x_discrete(name="Regions of interest") +
     scale_y_continuous(name="Scanner effect relative to unharmonised data",  
                       limits = c(-10.5, 10.5),
                       breaks = seq(-10, 10, by = 2)) +
     coord_flip() +
     theme_bw() +
     facet_grid(rows = vars(Metric), space = "free", scales = "free") +
     theme(legend.position = "bottom",
          strip.background = element_rect(fill="grey95"))
 
ggsave(plot = p, 
       filename = "Figures/forrest_plot_CoV.tiff",
       width = 16, height = 22, dpi = 300,
       units = "cm")

# Repeat for remaining methods as supplementary material

data <- 
  all_fu %>% 
  filter(Method %in% c("Unharmonised", "CrossCombat_nonpara", "CrossCombat_nonbays", "LongCombat_i+s"))
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
      
      res[[j]][(i-1)*3+1, "Coef"]   <- sum$coefficients["MethodCrossCombat_nonpara","Estimate"]
      res[[j]][(i-1)*3+2, "Coef"]   <- sum$coefficients["MethodCrossCombat_nonbays","Estimate"]
      res[[j]][(i-1)*3+3, "Coef"]   <- sum$coefficients["MethodLongCombat_i+s","Estimate"]
      res[[j]][c((i-1)*3+1,(i-1)*3+2, (i-1)*3+3), "ROI"]    <- rois[i]
      res[[j]][c((i-1)*3+1,(i-1)*3+2, (i-1)*3+3), "Metric"] <- metrics[j]
      res[[j]][(i-1)*3+1, "Method"] <- "CrossCombat_nonpara"
      res[[j]][(i-1)*3+2, "Method"] <- "CrossCombat_nonbays"
      res[[j]][(i-1)*3+3, "Method"] <- "LongCombat_i+s"
      res[[j]][(i-1)*3+1, "LLCI"]   <- confint(fit)["MethodCrossCombat_nonpara","2.5 %"]
      res[[j]][(i-1)*3+2, "LLCI"]   <- confint(fit)["MethodCrossCombat_nonbays","2.5 %"]
      res[[j]][(i-1)*3+3, "LLCI"]   <- confint(fit)["MethodLongCombat_i+s","2.5 %"]
      res[[j]][(i-1)*3+1, "ULCI"]   <- confint(fit)["MethodCrossCombat_nonpara","97.5 %"]
      res[[j]][(i-1)*3+2, "ULCI"]   <- confint(fit)["MethodCrossCombat_nonbays","97.5 %"]
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

data        <- res_methods
data$ROI    <- factor(data$ROI)
data$Method <- factor(data$Method, 
                        levels = c("CrossCombat_nonbays", "CrossCombat_nonpara", "LongCombat_i+s"),
                        labels = c("neuroCombat\n(non-bays)", "neuroCombat\n(non-para)", "LongCombat\n(i+s)"))
  
#define colours for dots and bars
dotCOLS = c("coral3", "darkslategray", "bisque4")
barCOLS = c("coral2", "darkslategray4", "bisque3")
  
  
p <- ggplot(data, aes(x=ROI, y=Coef, ymin=LLCI, ymax=ULCI, col=Method, fill=Method)) + 
     geom_hline(yintercept=0, size = 1, colour = "grey80") +
     #specify position here
     geom_linerange(size=1.8,position=position_dodge(width = 0.8)) +
     #specify position here too
     geom_point(size=1.8, shape=21, colour="white", stroke = 0.5 ,position=position_dodge(width = 0.8)) +
     scale_fill_manual(values=dotCOLS, breaks = c("LongCombat\n(i+s)", "neuroCombat\n(non-para)", "neuroCombat\n(non-bays)" ))+
     scale_color_manual(values=barCOLS, breaks = c("LongCombat\n(i+s)", "neuroCombat\n(non-para)", "neuroCombat\n(non-bays)"))+
     scale_x_discrete(name="Regions of interest") +
     scale_y_continuous(name="Scanner effect relative to unharmonised data",  
                       limits = c(-10.5, 10.5),
                       breaks = seq(-10, 10, by = 2)) +
     coord_flip() +
     theme_bw() +
     facet_grid(rows = vars(Metric), space = "free", scales = "free") +
     theme(legend.position = "bottom",
          strip.background = element_rect(fill="grey95"))
 
ggsave(plot = p, 
       filename = "Figures/forrest_plot_CoV_supplement.tiff",
       width = 16, height = 22, dpi = 300,
       units = "cm")

# Assessing false positive rate

# Get dataframes of unharmonised data in long format
vol_for_combat      <- merge(struc_long, vol, by = "Scan_ID", all.x = TRUE, all.y = FALSE)
corthick_for_combat <- merge(struc_long, corthick, by = "Scan_ID", all.x = TRUE, all.y = FALSE)
md_for_combat       <- merge(dti_long, md, by = "Merge_ID", all.x = TRUE, all.y = FALSE)
fa_for_combat       <- merge(dti_long, fa, by = "Merge_ID", all.x = TRUE, all.y = FALSE)

# Simulate and test for difference in INTERCEPT between groups

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
write.csv(fp_vol, "Data/fu_fp_vol.csv")

#simulate corthick
fp_unharm_corthick      <- sim_unharm(mymax = mymax, corthick_for_combat, 
                             is_corthick = "yes")
fp_crosscombat_corthick <- sim_crosscombat(mymax = mymax, corthick_for_combat, 
                                  from = "Frontal", to = "WholeCortex", 
                                  is_corthick = "yes")
fp_longcombat_corthick  <- sim_longcombat(mymax = mymax, corthick_for_combat,
                                 is_corthick = "yes")
fp_gamcombat_corthick   <- sim_gamcombat(mymax = mymax, corthick_for_combat, 
                                  from = "Frontal", to = "WholeCortex", 
                                  is_corthick = "yes")
fp_corthick             <- rbind(fp_unharm_corthick, 
                                 fp_crosscombat_corthick, 
                                 fp_longcombat_corthick,
                                 fp_gamcombat_corthick)
write.csv(fp_corthick, "Data/fu_fp_corthick.csv")

#simulate md
fp_unharm_md      <- sim_unharm(mymax = mymax, md_for_combat)
fp_crosscombat_md <- sim_crosscombat(mymax = mymax, md_for_combat, 
                                     is_fudti = "yes")
fp_longcombat_md  <- sim_longcombat(mymax = mymax, md_for_combat)
fp_gamcombat_md   <- sim_gamcombat(mymax = mymax, md_for_combat, 
                                   is_fudti = "yes")
fp_md             <- rbind(fp_unharm_md, 
                           fp_crosscombat_md, 
                           fp_longcombat_md,
                           fp_gamcombat_md)
write.csv(fp_md, "Data/fu_fp_md.csv")

#simulate fa
fp_unharm_fa      <- sim_unharm(mymax = mymax, fa_for_combat)
fp_crosscombat_fa <- sim_crosscombat(mymax = mymax, fa_for_combat, 
                                     is_fudti = "yes")
fp_longcombat_fa  <- sim_longcombat(mymax = mymax, fa_for_combat)
fp_gamcombat_fa   <- sim_gamcombat(mymax = mymax, fa_for_combat, 
                                   is_fudti = "yes")
fp_fa             <- rbind(fp_unharm_fa, 
                           fp_crosscombat_fa, 
                           fp_longcombat_fa,
                           fp_gamcombat_fa)
write.csv(fp_fa, "Data/fu_fp_fa.csv")


# Simulate and test for difference in SLOPE between groups

#simulate volume
fp_unharm_vol      <- sim_unharm(mymax = mymax, vol_for_combat, intercept_only = "no")
fp_crosscombat_vol <- sim_crosscombat(mymax = mymax, vol_for_combat, intercept_only = "no")
fp_longcombat_vol  <- sim_longcombat(mymax = mymax, vol_for_combat, intercept_only = "no")
fp_gamcombat_vol   <- sim_gamcombat(mymax = mymax, vol_for_combat, intercept_only = "no")
fp_vol             <- rbind(fp_unharm_vol, 
                            fp_crosscombat_vol, 
                            fp_longcombat_vol,
                            fp_gamcombat_vol)
write.csv(fp_vol, "Data/fu_fp_vol_slope.csv")

#simulate corthick
fp_unharm_corthick      <- sim_unharm(mymax = mymax, corthick_for_combat,
                             is_corthick = "yes", 
                             intercept_only = "no")
fp_crosscombat_corthick <- sim_crosscombat(mymax = mymax, corthick_for_combat, 
                                  from = "Frontal", to = "WholeCortex", 
                                  is_corthick = "yes", 
                                  intercept_only = "no")
fp_longcombat_corthick  <- sim_longcombat(mymax = mymax, corthick_for_combat,
                                 is_corthick = "yes", 
                                 intercept_only = "no")
fp_gamcombat_corthick   <- sim_gamcombat(mymax = mymax, corthick_for_combat, 
                                  from = "Frontal", to = "WholeCortex", 
                                  is_corthick = "yes", 
                                  intercept_only = "no")
fp_corthick             <- rbind(fp_unharm_corthick, 
                                 fp_crosscombat_corthick, 
                                 fp_longcombat_corthick,
                                 fp_gamcombat_corthick)

write.csv(fp_corthick, "Data/fu_fp_corthick_slope.csv")

#simulate md
fp_unharm_md      <- sim_unharm(mymax = mymax, md_for_combat, 
                                intercept_only = "no")
fp_crosscombat_md <- sim_crosscombat(mymax = mymax, md_for_combat, 
                                     is_fudti = "yes", 
                                     intercept_only = "no")
fp_longcombat_md  <- sim_longcombat(mymax = mymax, md_for_combat, 
                                    intercept_only = "no")
fp_gamcombat_md   <- sim_gamcombat(mymax = mymax, md_for_combat, 
                                   is_fudti = "yes", 
                                   intercept_only = "no")
fp_md             <- rbind(fp_unharm_md, 
                           fp_crosscombat_md, 
                           fp_longcombat_md,
                           fp_gamcombat_md)
write.csv(fp_md, "Data/fu_fp_md_slope.csv")

#simulate fa
fp_unharm_fa      <- sim_unharm(mymax = mymax, fa_for_combat, 
                                intercept_only = "no")
fp_crosscombat_fa <- sim_crosscombat(mymax = mymax, fa_for_combat, 
                                     is_fudti = "yes", 
                                     intercept_only = "no")
fp_longcombat_fa  <- sim_longcombat(mymax = mymax, fa_for_combat, 
                                    intercept_only = "no")
fp_gamcombat_fa   <- sim_gamcombat(mymax = mymax, fa_for_combat, 
                                   is_fudti = "yes", 
                                   intercept_only = "no")
fp_fa             <- rbind(fp_unharm_fa, 
                           fp_crosscombat_fa, 
                           fp_longcombat_fa,
                           fp_gamcombat_fa)
write.csv(fp_fa, "Data/fu_fp_fa_slope.csv")

