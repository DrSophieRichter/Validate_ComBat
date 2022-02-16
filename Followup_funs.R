#Functions for the longitudinal analysis

####################################################
#                                                  #
# Select structural scans for longitudinal cohort  #
#                                                  #
####################################################


select_struc_fu <- function(scans, qc_bad, clinical){
 
  df <- scans %>% 
    filter(Category %in% c("HEALTHY", "ORTHO") & 
             Cohort != "CamCAN" & 
             Malpem138_volume == "available" &
             !Scan_ID %in% qc_bad$Scan_ID) 
  
  #define Scanner by Site and manufacturer model
  #acknowledge that Site Site-06-a72b20 is Cambridge
  df$Site <- ifelse(df$Site %in% c("Site-06-a72b20", "detect/legacy"), "Cambridge", df$Site)
  df$Scanner <- paste0(df$Site, "_", df$Model)
  
  #keep only subjects with at least 3 scans of which 2 scans on the same scanner
  df <- 
    df %>% 
    select(Master_subject_ID, 
           Scan_ID, 
           Days_since_first_scan, 
           Scan_date, 
           Scanner)
  temp1 <- 
    df %>% 
    group_by(Master_subject_ID, Scanner) %>% 
    summarise(num = n()) %>% 
    filter(num > 1)
  temp2 <- 
    df %>% 
    group_by(Master_subject_ID) %>% 
    summarise(num = n()) %>% filter(num > 2)
  df <- 
    df %>% 
    filter(Master_subject_ID %in% temp1$Master_subject_ID & Master_subject_ID %in% temp2$Master_subject_ID)
  
  #to temp1, add the time difference between the two scans done on the same scanner
  #select subjects with 2 scans done on the same scanner
  temp1 <- 
    temp1 %>% 
    filter(Master_subject_ID %in% df$Master_subject_ID)
  #make "within" listing all possible same-scanner pairs per eligible subject
  temp1$Subject_Scanner <- paste0(temp1$Master_subject_ID, "_", temp1$Scanner)
  df$Subject_Scanner    <- paste0(df$Master_subject_ID, "_", df$Scanner)
  within <- df %>% 
    filter(Subject_Scanner %in% temp1$Subject_Scanner)
 
   #out of all scans per subject with the same Subject_Scanner (i.e. all done on the same scanner), pick the two furthest apart
  withinmax <- 
    within %>% 
    group_by(Subject_Scanner) %>% 
    slice(which.max(Days_since_first_scan))
  withinmin <- 
    within %>% 
    group_by(Subject_Scanner) %>% 
    slice(which.min(Days_since_first_scan))
  withinmin <- 
    withinmin %>% 
    ungroup %>% 
    select(Master_subject_ID, 
           Scan_ID_withinmin = Scan_ID, 
           Days_withinmin = Days_since_first_scan, 
           Scanner_within = Scanner,
           Subject_Scanner_within = Subject_Scanner)
  withinmax <- 
    withinmax %>% 
    ungroup %>% 
    select(Master_subject_ID, 
           Scan_ID_withinmax = Scan_ID, 
           Days_withinmax = Days_since_first_scan, 
           Scanner_within = Scanner,
           Subject_Scanner_within = Subject_Scanner)
  
  
  within <- merge(withinmax, withinmin, 
                  by = c("Master_subject_ID", "Scanner_within", "Subject_Scanner_within"), 
                  all.x = TRUE, all.y = TRUE)
  within$DeltaDays_within <- within$Days_withinmax - within$Days_withinmin
  within$Pair_ID          <- paste0(within$Scan_ID_withinmin, "_", within$Scan_ID_withinmax)
  
  #pick all subjects that have a difference of more than 365 days between the two scans on the same scanner
  within <- 
    within %>% 
    filter(DeltaDays_within > 365)
  n_distinct(within$Master_subject_ID) # 21 distinct subjects
  
  #pick the scan pair with the biggest time gap per subject (but keep alternatives)
  within_chosen <- 
    within %>% 
    group_by(Master_subject_ID) %>% 
    slice(which.max(DeltaDays_within))
  alternative   <- 
    within %>% 
    filter(!Pair_ID %in% within_chosen$Pair_ID)
  
  #Now look at the third scan
  across <- 
    df %>% 
    filter(Master_subject_ID %in% within_chosen$Master_subject_ID & !Subject_Scanner %in% within_chosen$Subject_Scanner_within)
  n_distinct(across$Master_subject_ID)#17
  across <- 
    across %>% 
    select(Master_subject_ID, 
           Scan_ID_across = Scan_ID, 
           Days_since_first_scan_across = Days_since_first_scan,
           Scanner_across = Scanner,
           Subject_Scanner_across = Subject_Scanner)
  across <- merge(across, within_chosen, by = "Master_subject_ID", all.x = TRUE, all.y = FALSE)
  
  #The third scan must have happened either 365+ days after the first within-scan (if the first within-scan is the reference scan)
  #or 365+ days before the second within-scans (if the second within-scan is the reference)
  #As I want both follow up scans to have either happened both after the refrence scan or both before the reference scan
  #as to maximise the overlap in the time window of  a patients life they cover
  across$DeltaDays_after_withinmin  <- across$Days_since_first_scan_across - across$Days_withinmin
  across$DeltaDays_before_withinmax <- across$Days_withinmax - across$Days_since_first_scan_across
  across <- 
    across %>% 
    filter(DeltaDays_after_withinmin > 365 | DeltaDays_before_withinmax > 365) # do not use abs() here as I don't want days in the negative direction
  #all across scans are at least 365 days from one of the within scans, so I do not have to exclude any
  #So I do not have to go back to my alternative scan pairs
  
  #For each subject, pick the across follow up scan that is closes to the within follow up scan
  across <- 
    across %>% 
    rowwise %>% 
    mutate(DeltaDays_across = max(DeltaDays_after_withinmin, DeltaDays_before_withinmax))
  across$DeltaDelta <- abs(across$DeltaDays_within - across$DeltaDays_across)
  across <- 
    across %>% 
    group_by(Master_subject_ID) %>% 
    slice(which.min(DeltaDelta))
  
  #identify the reference scan for each subject
  across$Ref               <- ifelse(across$DeltaDays_across == across$DeltaDays_after_withinmin, "within_min", "within_max")
  across$Scan_ID_ref       <- ifelse(across$Ref == "within_min", across$Scan_ID_withinmin, across$Scan_ID_withinmax)
  across$Scan_ID_fu_within <- ifelse(across$Ref == "within_min", across$Scan_ID_withinmax, across$Scan_ID_withinmin)
  across$Scan_ID_fu_across <- across$Scan_ID_across
  #correct time signs based on the reference scan
  across$Time_fu_within <- ifelse(across$Ref == "within_min", across$DeltaDays_within, -1*across$DeltaDays_within)
  across$Time_fu_across <- ifelse(across$Ref == "within_min", across$DeltaDays_across, -1*across$DeltaDays_across)
  
  struc_wide <- 
    across %>% 
    select(Master_subject_ID, 
           Scan_ID_ref, 
           Scan_ID_fu_within, 
           Scan_ID_fu_across, 
           Time_fu_within, 
           Time_fu_across, 
           Scanner_within, 
           Scanner_across)
  
  
  return(struc_wide)
}

####################################################
#                                                  #
# Wide to long format for structural data          #
#                                                  #
####################################################

wide_to_long_struc <- function(struc_wide, scans, clinical){
  
  struc_long <- 
    struc_wide %>% 
    gather(key = "Label", value = "Scan_ID", contains("Scan_ID"))
  
  struc_long$Label <- gsub("Scan_ID_", "", struc_long$Label)
  
  struc_long$Time <- ifelse(struc_long$Label == "ref", 0,
                            ifelse(struc_long$Label == "fu_within", struc_long$Time_fu_within, struc_long$Time_fu_across))
  
  struc_long$Scanner <- ifelse(struc_long$Label == "ref", struc_long$Scanner_within,
                               ifelse(struc_long$Label == "fu_within", struc_long$Scanner_within, struc_long$Scanner_across))
  
  struc_long <- 
    struc_long %>% 
    select(Master_subject_ID, 
           Label, Scan_ID, 
           Scanner, 
           Time)
  
  #check if each scanner has at least 2 scans
  summary(as.factor(struc_long$Scanner)) #good, at least 5 scans per Scanner
  
  #Now add age and sex
  age <- 
    scans %>% 
    select(Scan_ID, 
           Age_at_this_scan)
  
  sex <- 
    clinical %>% 
    select(Master_subject_ID, 
           Sex)
  sex$Sex <- ifelse(sex$Sex == "female", 0, 1)
  
  struc_long <- merge(age, struc_long, by = "Scan_ID", all.x = FALSE, all.y = TRUE)
  
  struc_long <- merge(sex, struc_long, by = "Master_subject_ID", all.x = FALSE, all.y = TRUE)
  
  return(struc_long)
}

####################################################
#                                                  #
# Select diffusion scans for longitudinal cohort   #
#                                                  #
####################################################


select_dti_fu <- function(scans, qc_bad, clinical){
  #Only healthy controls with available structural imaging also
  df <- scans %>% 
    filter(Category %in% c("HEALTHY", "ORTHO") & 
             Cohort != "CamCAN" & 
             Malpem138_volume == "available" &
             !Scan_ID %in% qc_bad$Scan_ID)
  
  #select only scans with diffusion and create one entry per series
  #as sometimes more than one series was acquired per scan
  df <- merge(df, md, by = "Scan_ID", all.x = FALSE, all.y = TRUE)
  df$Scanner <- paste0(df$Model,"_", df$Series)
  
  #Make sure that for the same day and series I only have one copy
  df <- df %>% group_by(Scan_ID, Scanner) %>% slice(1)
  
  #keep only subjects with at least 3 scans of which 2 scans on the same scanner
  df <- df %>% select(Master_subject_ID, Scan_ID, Merge_ID, Days_since_first_scan, Scan_date, Scanner)
  temp1 <- df %>% group_by(Master_subject_ID, Scanner) %>% summarise(num = n()) %>% filter(num > 1)
  temp2 <- df %>% group_by(Master_subject_ID) %>% summarise(num = n()) %>% filter(num > 2)
  df <- df %>% filter(Master_subject_ID %in% temp1$Master_subject_ID & Master_subject_ID %in% temp2$Master_subject_ID)
  
  #to temp1, add the time difference between the two scans done on the same scanner
  #select subjects with 2 scans done on the same scanner
  temp1 <- temp1 %>% filter(Master_subject_ID %in% df$Master_subject_ID)
  #make "within" listing all possible same-scanner pairs per eligible subject
  temp1$Subject_Scanner <- paste0(temp1$Master_subject_ID, "_", temp1$Scanner)
  df$Subject_Scanner <- paste0(df$Master_subject_ID, "_", df$Scanner)
  within <- df %>% filter(Subject_Scanner %in% temp1$Subject_Scanner)
  #out of all scans per subject with the same Subject_Scanner (i.e. all done on the same scanner), pick the two furthest apart
  withinmax <- within %>% group_by(Subject_Scanner) %>% slice(which.max(Days_since_first_scan))
  withinmin <- within %>% group_by(Subject_Scanner) %>% slice(which.min(Days_since_first_scan))
  withinmax <- 
    withinmax %>% 
    ungroup %>% 
    select(Master_subject_ID, 
           Merge_ID_withinmax = Merge_ID,
           Days_withinmax = Days_since_first_scan, 
           Scanner_within = Scanner,
           Subject_Scanner_within = Subject_Scanner)
  withinmin <- 
    withinmin %>% 
    ungroup %>% 
    select(Master_subject_ID, 
           Merge_ID_withinmin = Merge_ID, 
           Days_withinmin = Days_since_first_scan, 
           Scanner_within = Scanner,
           Subject_Scanner_within = Subject_Scanner)
  within <- merge(withinmax, withinmin, 
                  by = c("Master_subject_ID", "Scanner_within", "Subject_Scanner_within"), 
                  all.x = TRUE, all.y = TRUE)
  within$DeltaDays_within <- within$Days_withinmax - within$Days_withinmin
  within$Pair_ID <- paste0(within$Merge_ID_withinmin, "_", within$Merge_ID_withinmax)
  
  #pick all subjects that have a difference of more than 365 days between the two scans on the same scanner
  within <- within %>% filter(DeltaDays_within > 365)
  n_distinct(within$Master_subject_ID) # 14 distinct subjects
  
  #pick the scan pair with the biggest time gap per subject (but keep alternatives)
  within_chosen <- within %>% group_by(Master_subject_ID) %>% slice(which.max(DeltaDays_within))
  alternative <- within %>% filter(!Pair_ID %in% within_chosen$Pair_ID)
  
  #Now look at the third scan
  across <- df %>% filter(Master_subject_ID %in% within_chosen$Master_subject_ID & !Subject_Scanner %in% within_chosen$Subject_Scanner_within)
  n_distinct(across$Master_subject_ID)#14
  across <- across %>% select(Master_subject_ID, 
                              Merge_ID_across = Merge_ID, 
                              Days_since_first_scan_across = Days_since_first_scan,
                              Scanner_across = Scanner,
                              Subject_Scanner_across = Subject_Scanner)
  across <- merge(across, within_chosen, by = "Master_subject_ID", all.x = TRUE, all.y = FALSE)
  
  
  #The third scan must have happened either 365+ days after the first within-scan (if the first within-scan is the reference scan)
  #or 365+ days before the second within-scans (if the second within-scan is the reference)
  #As I want both follow up scans to have either happened both after the refrence scan or both before the reference scan
  #as to maximise the overlap in the time window of  a patients life they cover
  across$DeltaDays_after_withinmin <- across$Days_since_first_scan_across - across$Days_withinmin
  across$DeltaDays_before_withinmax <- across$Days_withinmax - across$Days_since_first_scan_across
  across <- across %>% filter(DeltaDays_after_withinmin > 365 | DeltaDays_before_withinmax > 365) # do not use abs() here as I don't want days in the negative direction
  #all across scans are at least 365 days from one of the within scans, so I do not have to exclude any
  #So I do not have to go back to my alternative scan pairs
  
  across <- across %>% rowwise %>% mutate(DeltaDays_across = max(DeltaDays_after_withinmin, DeltaDays_before_withinmax))
  
  
  #For each subject, pick the across follow up scan that is closes to the within follow up scan
  #but keep alternative options for later
  across$DeltaDelta <- abs(across$DeltaDays_within - across$DeltaDays_across)
  alternative_across <- across
  across <- across %>% group_by(Master_subject_ID) %>% slice(which.min(DeltaDelta))
  alternative_across <- alternative_across %>% filter(!Merge_ID_across %in% across$Merge_ID_across)
  #identify the reference scan for each subject
  across$Ref <- ifelse(across$DeltaDays_across == across$DeltaDays_after_withinmin, "within_min", "within_max")
  across$Merge_ID_ref <- ifelse(across$Ref == "within_min", across$Merge_ID_withinmin, across$Merge_ID_withinmax)
  across$Merge_ID_fu_within <- ifelse(across$Ref == "within_min", across$Merge_ID_withinmax, across$Merge_ID_withinmin)
  across$Merge_ID_fu_across <- across$Merge_ID_across
  #correct time signs based on the reference scan
  across$Time_fu_within <- ifelse(across$Ref == "within_min", across$DeltaDays_within, -1*across$DeltaDays_within)
  across$Time_fu_across <- ifelse(across$Ref == "within_min", across$DeltaDays_across, -1*across$DeltaDays_across)
  
  dti_wide <- 
    across %>% 
    select(Master_subject_ID, 
           Merge_ID_ref, 
           Merge_ID_fu_within, 
           Merge_ID_fu_across, 
           Time_fu_within, 
           Time_fu_across, 
           Scanner_within, 
           Scanner_across)
  
  dti_long <- dti_wide %>% gather(key = "Label", value = "Merge_ID", contains("Merge_ID"))
  dti_long$Label <- gsub("Merge_ID_", "", dti_long$Label)
  dti_long$Time <- ifelse(dti_long$Label == "ref", 0,
                          ifelse(dti_long$Label == "fu_within", dti_long$Time_fu_within, dti_long$Time_fu_across))
  dti_long$Scanner <- ifelse(dti_long$Label == "ref", dti_long$Scanner_within,
                             ifelse(dti_long$Label == "fu_within", dti_long$Scanner_within, dti_long$Scanner_across))
  dti_long <- dti_long %>% select(Master_subject_ID, Label, Merge_ID, Scanner, Time)
  
  #check if each scanner has at least 2 scans
  summary(as.factor(dti_long$Scanner)) #oh no, problems!
  #Single scan per scanner in 4 instances!
  scanners <- summary(as.factor(dti_long$Scanner)) %>% as.data.frame()
  colnames(scanners) <- "count"
  single_scanners <- scanners %>% filter(count == 1)
  single_scanners
  problem <- dti_long %>% filter(Scanner %in% row.names(single_scanners)) #for all 4 subjects its the fu_across scan that is the problem
  #first lets check if they have a different option for the fu_across scan
  alternative_across <- alternative_across %>% filter(Master_subject_ID %in% problem$Master_subject_ID)
  #any subject who doesn't have an alternative across_fu?
  problem %>% filter(!Master_subject_ID %in% alternative_across$Master_subject_ID)
  #can I avoid the Prisma cbu multishell?
  temp <- scans %>% filter(Master_subject_ID == "21091") %>% select(Scan_ID, contains("DTIPATH"))
  #nope, this subject has only 3 scans
  #can I make one of the other fu_across Prisma cbu multishell?
  #yes could make 1000 and 10179 CBU multishell (albeit deltadelta being very large, esp for 1000)
  #that would sort 3/4
  #fourth is 18545 Verio free 12
  #I could use Prisma_proc_set1_nobzero instead as this already exists in my set of scanners
  
  #I will replace 3 fu_across scan and keep the 4th one (which now has other scans done on the same scanner)
  replacement <- 
    alternative_across %>% 
    filter(Subject_Scanner_across %in% c("10000_Prisma_cbu_multi-shell_sequence",
                                         "10179_Prisma_cbu_multi-shell_sequence",
                                         "18545_Prisma_proc_set1_nobzero"))
  #all 3 replacement scans are measured with within_min as the reference scan
  #I will consider this when renaming the columns to match dti_wide
  replacement <- replacement %>% select(Master_subject_ID, 
                                        Merge_ID_ref = Merge_ID_withinmin,
                                        Merge_ID_fu_within = Merge_ID_withinmax,
                                        Merge_ID_fu_across = Merge_ID_across,
                                        Time_fu_within = DeltaDays_within,
                                        Time_fu_across = DeltaDays_across,
                                        Scanner_within,
                                        Scanner_across)
  #remove the row to be replaced from dti_wide
  dti_wide <- dti_wide %>% filter(!Master_subject_ID %in% replacement$Master_subject_ID)
  dti_wide <- rbind(dti_wide, replacement) %>% as.data.frame()
  
  
  return(dti_wide)
}

####################################################
#                                                  #
# Wide to long format for dti data                 #
#                                                  #
####################################################

wide_to_long_dti <- function(dti_wide, scans, clinical){
  
  #now re-arrange dti_wide into dti_long
  dti_long <- 
    dti_wide %>% 
    gather(key = "Label", value = "Merge_ID", contains("Merge_ID"))
  dti_long$Label <- gsub("Merge_ID_", "", dti_long$Label)
  dti_long$Time <- ifelse(dti_long$Label == "ref", 0,
                          ifelse(dti_long$Label == "fu_within", dti_long$Time_fu_within, dti_long$Time_fu_across))
  dti_long$Scanner <- ifelse(dti_long$Label == "ref", dti_long$Scanner_within,
                             ifelse(dti_long$Label == "fu_within", dti_long$Scanner_within, dti_long$Scanner_across))
  dti_long <- dti_long %>% select(Master_subject_ID, Label, Merge_ID, Scanner, Time)
  
  #check if each scanner has at least 2 scans
  summary(as.factor(dti_long$Scanner)) #ok, at least 3 scans per scanner
  
  #Now add age and sex
  #all Merge_IDs are in the Cambridge format
  dti_long <- 
    dti_long %>% 
    separate(Merge_ID, into = c("wbic", "date", "uid"), remove = FALSE, extra = "drop", sep = "_")
  dti_long <- 
    dti_long %>% 
    unite("Scan_ID", wbic:uid)
  
  age <- 
    scans %>% 
    select(Scan_ID, Age_at_this_scan)
  sex <- 
    clinical %>% 
    select(Master_subject_ID, Sex)
  sex$Sex <- ifelse(sex$Sex == "female", 0, 1)
  
  dti_long <- merge(age, dti_long, by = "Scan_ID", all.x = FALSE, all.y = TRUE)
  dti_long <- merge(sex, dti_long, by = "Master_subject_ID", all.x = FALSE, all.y = TRUE)
  
  return(dti_long)
}

####################################################
#                                                  #
# Make summary table for cohort characteristics    #
#                                                  #
####################################################


make_cohort_table_2 <- function(struc_wide, struc_long){
  
  df <- df <- as.data.frame(matrix(,ncol=1,nrow=1))
  colnames(df) <- "Subjects"
  
  df$Subjects <- n_distinct(struc_long$Master_subject_ID)
  df$Scans <- nrow(struc_long)
  df$Scanners <- n_distinct(struc_long$Scanner)
  
  df$Time_within_median <- format(round(median(abs(struc_wide$Time_fu_within))/365, 1), nsmall = 1) 
  df$Time_within_min <- format(round(min(abs(struc_wide$Time_fu_within))/365, 1), nsmall = 1)
  df$Time_within_max <- format(round(max(abs(struc_wide$Time_fu_within))/365, 1), nsmall = 1)
  df$Time_fu_within <- paste0(df$Time_within_median, " (", df$Time_within_min, " - ", df$Time_within_max, ") years")
  
  df$Time_across_median <- format(round(median(abs(struc_wide$Time_fu_across))/365, 1), nsmall = 1) 
  df$Time_across_min <- format(round(min(abs(struc_wide$Time_fu_across))/365, 1), nsmall = 1)
  df$Time_across_max <- format(round(max(abs(struc_wide$Time_fu_across))/365, 1), nsmall = 1)
  df$Time_fu_across <- paste0(df$Time_across_median, " (", df$Time_across_min, " - ", df$Time_across_max, ") years")
  
  df <- df %>% select(-c(Time_within_median, Time_within_min, Time_within_max, 
                         Time_across_median, Time_across_min, Time_across_max))
  
  df$Age <- paste0(median(struc_long$Age_at_this_scan) %>% round(0), 
                   " (",
                   min(struc_long$Age_at_this_scan),
                   " - ",
                   max(struc_long$Age_at_this_scan),
                   ")")
  
  df$Sex <- paste0(sum(struc_long$Sex==1)/3, 
                   " (", 
                   ((sum(struc_long$Sex==1)/3)/df$Subjects*100) %>% round(0), 
                   "%)")
  
  df <- t(df) %>% as.data.frame()
  
  return(df)
}

####################################################
#                                                  #
# Calculate scanner induced difference             #
#                                                  #
####################################################

calculate_scanner_induced_difference <- function(vol_fu_unharm, struc_wide, from = "Ventricles", to = "Brainstem"){
  
  #bring these data into the wide format, so that there is one row per subject per roi
  #and one column per scan label (ref, fu_within, fu_across)
  vol_fu_unharm <- vol_fu_unharm %>% gather(key = "ROI", value = "Metric", from:to)
  vol_fu_unharm <- vol_fu_unharm %>% select(Master_subject_ID, Label, ROI, Metric)
  vol_fu_unharm <- vol_fu_unharm %>% spread(key = Label, value = Metric)
  vol_fu_unharm[,c(3:5)] <- sapply(vol_fu_unharm[,c(3:5)], as.numeric)
  
  #Add time differences (in years)
  time <- 
    struc_wide %>% 
    select(Master_subject_ID, Time_fu_within, Time_fu_across)
  time$Time_fu_across <- time$Time_fu_across/365
  time$Time_fu_within <- time$Time_fu_within/365
  vol_fu_unharm <- merge(vol_fu_unharm, time, by = "Master_subject_ID", all.x = TRUE, all.y = FALSE)
  
  #calculate the change from the initial to the follow up scan PER year
  #careful with signs as sometimes the reference scan is AFTER the follow up scan
  vol_fu_unharm$Rate_within <- (vol_fu_unharm$fu_within - vol_fu_unharm$ref)/vol_fu_unharm$Time_fu_within
  vol_fu_unharm$Rate_across <- (vol_fu_unharm$fu_across - vol_fu_unharm$ref)/vol_fu_unharm$Time_fu_across
  
  #the scanner induced discrepancy between within and across scan pairs needs to be
  #expressed as a percentage CoV 
  #to adjust for the fact that there is a time difference between the within and across
  #follow-up scans I will calculate what the across follow-up scan
  #would have been if it had been done exactly the same day as the within scan
  #based on the rate calculate for the across scan
  vol_fu_unharm$fu_across_adj <- vol_fu_unharm$ref + vol_fu_unharm$Rate_across*vol_fu_unharm$Time_fu_within
  fun <- function(x){
    (sd(x)/mean(x))*100
  }
  vol_fu_unharm$CoV <- apply(vol_fu_unharm[,c("fu_within", "fu_across_adj")], 1, fun)
  
  return(vol_fu_unharm)
}

