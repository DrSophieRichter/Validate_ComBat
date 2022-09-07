# Functions to select scans for volumetric analysis


##########################################
#                                        #
# Collapse ROI volumes                   #
#                                        #
##########################################

collapse_vol <- function(vol, mymalp){
  
  #Make my combined ROIs
  mymalp$ROI.Index <- paste0("ROI_", mymalp$ROI.Index)
  vol$Ventricles <- vol$Ventricle
  vol$SupraCortex <- vol$Cortical
  vol$SupraWM <- vol$ROI_16 + vol$ROI_17
  dgm <- mymalp %>% filter(Sophie == "SupDGM")
  vol$SupraDGM <- rowSums(vol[, which(colnames(vol) %in% dgm$ROI.Index)])
  clumGM <- mymalp %>% filter(Sophie == "CerebellumGM")
  vol$CerebellumGM <- rowSums(vol[, which(colnames(vol) %in% clumGM$ROI.Index)])
  clumWM <- mymalp %>% filter(Sophie == "CerebellumWM")
  vol$CerebellumWM <- rowSums(vol[, which(colnames(vol) %in% clumWM$ROI.Index)])
  vol$Brainstem <- vol$ROI_7
  
  vol <- vol %>% select(Scan_ID,  
                        Ventricles, 
                        SupraCortex, 
                        SupraWM, 
                        SupraDGM,
                        CerebellumWM,
                        CerebellumGM,
                        Brainstem)
  
  return(vol)
}



##########################################
#                                        #
# Collapse cortical thickness ROIs       #
#                                        #
##########################################

collapse_corthick <- function(corthick, add){
  corthick$Frontal      <- (corthick$ROI_6  + corthick$ROI_15)/2
  corthick$Insular      <- (corthick$ROI_7  + corthick$ROI_16)/2
  corthick$Occipital    <- (corthick$ROI_8  + corthick$ROI_17)/2
  corthick$Parietal     <- (corthick$ROI_9  + corthick$ROI_18)/2
  corthick$Hippocampal  <- (corthick$ROI_10 + corthick$ROI_19)/2
  corthick$Temporal     <- (corthick$ROI_11 + corthick$ROI_20)/2
  
  corthick <- corthick %>% select(Scan_ID, 
                                  Frontal,
                                  Insular,
                                  Occipital,
                                  Parietal,
                                  Hippocampal,
                                  Temporal)
  
  add <- add %>% select(Scan_ID, WholeCortex = Cortical)
  add <- add %>% distinct(Scan_ID, .keep_all = TRUE)
  
  corthick <- merge(corthick, add, by = "Scan_ID", all.x = TRUE, all.y = FALSE)
  
  corthick <- corthick %>% select(Scan_ID,  
                                  Frontal, 
                                  Parietal, 
                                  Insular, 
                                  Occipital,
                                  Hippocampal,
                                  Temporal,
                                  WholeCortex)
  
  return(corthick)
}


########################################
#                                      #
# IDENTIFY THE WITHIN-SCANNER COHORT   #
#                                      #
########################################

select_struc_within <- function(scans, qc_bad){
  
  #select only non-head injured subjects with available structural imaging
  df <- scans %>% 
    filter(Category %in% c("HEALTHY", "ORTHO") & 
             Cohort != "CamCAN" & 
             Malpem138_volume == "available" &
             !Scan_ID %in% qc_bad$Scan_ID &
             Scan_ID != "25041_20170519_U-ID38019") #otherwise I end up with only one scan from this scanner
  
  #define Scanner by Site and manufacturer model
  #acknowledge that Site Site-06-a72b20 is Cambridge
  df$Site <- ifelse(df$Site %in% c("Site-06-a72b20", "detect/legacy"), "Cambridge", df$Site)
  df$Scanner <- paste0(df$Site, "_", df$Model)
  
  #identify subjects with more than one scan per scanner
  temp <- df %>% 
    group_by(Master_subject_ID, Scanner) %>% summarise(num = n()) %>% filter(num >1)
  temp$Series_ID <- paste0(temp$Master_subject_ID, "_", temp$Scanner)
  df$Series_ID <- paste0(df$Master_subject_ID, "_", df$Scanner)
  rescan_ever <- df %>% filter(Series_ID %in% temp$Series_ID)
  
  #Calculate the time interval between scans of the same series
  rescan_ever <- rescan_ever %>% 
    group_by(Series_ID) %>%
    arrange(Days_since_first_scan, .by_group = T)%>%
    mutate(Days_since_previous_scan = Days_since_first_scan - lag(Days_since_first_scan)) 
  rescan_ever <- rescan_ever %>% select(Master_subject_ID, Scan_ID, Category, Days_since_first_scan, Days_since_previous_scan, Scanner)
  
  #rank the scans in chronological order
  rescan_ever <- rescan_ever %>%
    group_by(Series_ID) %>%
    mutate(Chron_rank = order(Days_since_first_scan, decreasing=FALSE))
  rescan_ever$Rank_ID <- paste0(rescan_ever$Series_ID, "_", rescan_ever$Chron_rank)
  
  #Identify all scans that have been done within 6 months of the previous scan
  rescan_6m <- rescan_ever %>% filter(Days_since_previous_scan <= 180)
  rescan_6m$Prev_scan <- paste0(rescan_6m$Series_ID, "_", (rescan_6m$Chron_rank -1))
  
  #Keep only one scan-pair per subject (the one with the shortest interval)
  keep <-  
    rescan_6m %>% 
    group_by(Master_subject_ID) %>% 
    slice(which.min(Days_since_previous_scan))
  
  #Add Scan and Rank IDs for second scan
  temp <- rescan_ever %>% select(Scan_ID2 = Scan_ID, Rank_ID2 = Rank_ID, Series_ID2 = Series_ID)
  keep <- merge(temp, keep, by.x = "Rank_ID2", by.y = "Prev_scan", all.x = FALSE, all.y = TRUE)
  
  return(keep)

}


########################################
#                                      #
# ADD VOLUME TO WITHIN-SCANNER COHORT  #
#                                      #
########################################

add_vol_within <- function(keep, vol){
  
  #Add volume for first scan
  s1 <- vol %>% filter(Scan_ID %in% keep$Scan_ID)
  colnames(s1)[-1] <- paste0("S1_", colnames(s1)[-1]) 
  keep <- merge(keep, s1, by = "Scan_ID", all.x = TRUE, all.y = FALSE)
  
  #Add volume for second scan
  s2 <- vol %>% filter(Scan_ID %in% keep$Scan_ID2)
  colnames(s2)[-1] <- paste0("S2_", colnames(s2)[-1]) 
  keep <- merge(keep, s2, by.x = "Scan_ID2", by.y = "Scan_ID", all.x = TRUE, all.y = FALSE)
  
  #isolate first scan
  s1 <- keep %>% select(Master_subject_ID, Scan_ID, Scanner, S1_Ventricles:S1_Brainstem)
  s1$Days <- 0
  s1$Scan_num <- 1
  s1 <- s1 %>% select(Master_subject_ID, Scanner, Scan_ID, Scan_num, Days, everything())
  
  #isolate second scan
  s2 <- keep %>% 
    select(Master_subject_ID, 
           Scanner, 
           Scan_ID = Scan_ID2, 
           Days = Days_since_previous_scan, 
           S2_Ventricles:S2_Brainstem)
  s2$Scan_num <- 2
  s2 <- s2 %>% select(Master_subject_ID, Scanner, Scan_ID, Scan_num, Days, everything())
  
  #combine into long format
  colnames(s1)[6:12] <- gsub("S1_","", colnames(s1)[6:12])
  colnames(s2)[6:12] <- gsub("S2_","", colnames(s2)[6:12])
  slong <- rbind(s1, s2) %>% as.data.frame()
  
  return(slong)
  
}

####################################################
#                                                  #
# ADD CORTICAL THICKNESS TO WITHIN-SCANNER COHORT  #
#                                                  #
####################################################

add_corthick_within <- function(keep, corthick){
  
  #Add corthick for first scan
  s1 <- corthick %>% filter(Scan_ID %in% keep$Scan_ID)
  colnames(s1)[-1] <- paste0("S1_", colnames(s1)[-1]) 
  keep <- merge(keep, s1, by = "Scan_ID", all.x = TRUE, all.y = FALSE)
  
  #Add corthick for second scan
  s2 <- corthick %>% filter(Scan_ID %in% keep$Scan_ID2)
  colnames(s2)[-1] <- paste0("S2_", colnames(s2)[-1]) 
  keep <- merge(keep, s2, by.x = "Scan_ID2", by.y = "Scan_ID", all.x = TRUE, all.y = FALSE)
  
  #isolate first scan
  s1 <- keep %>% select(Master_subject_ID, Scan_ID, Scanner, S1_Frontal:S1_WholeCortex)
  s1$Days <- 0
  s1$Scan_num <- 1
  s1 <- s1 %>% select(Master_subject_ID, Scanner, Scan_ID, Scan_num, Days, everything())
  
  #isolate second scan
  s2 <- keep %>% 
    select(Master_subject_ID, 
           Scanner, 
           Scan_ID = Scan_ID2, 
           Days = Days_since_previous_scan, 
           S2_Frontal:S2_WholeCortex)
  s2$Scan_num <- 2
  s2 <- s2 %>% select(Master_subject_ID, Scanner, Scan_ID, Scan_num, Days, everything())
  
  #combine into long format
  colnames(s1)[6:12] <- gsub("S1_","", colnames(s1)[6:12])
  colnames(s2)[6:12] <- gsub("S2_","", colnames(s2)[6:12])
  slong <- rbind(s1, s2) %>% as.data.frame()
  
  return(slong)
  
}



  
  ##########################################
  #                                        #
  # IDENTIFY THE ACROSS-SCANNER COHORT     #
  #                                        #
  ##########################################

select_struc_across <- function(scans, qc_bad){
  
  #select only non-head injured subjects with available structural imaging
  df <- scans %>% 
    filter(Category %in% c("HEALTHY", "ORTHO") & 
             Cohort != "CamCAN" & 
             Malpem138_volume == "available" &
             !Scan_ID %in% qc_bad$Scan_ID &
             Scan_ID != "25041_20170519_U-ID38019") #otherwise I end up with only one scan from this scanner
  
  #define Scanner by Site and manufacturer model
  #acknowledge that Site Site-06-a72b20 is Cambridge
  df$Site <- ifelse(df$Site %in% c("Site-06-a72b20", "detect/legacy"), "Cambridge", df$Site)
  df$Scanner <- paste0(df$Site, "_", df$Model)
  
  #select only subjects with at least 2 scans
  temp <- df %>% group_by(Master_subject_ID) %>% summarise(num = n()) %>% filter(num>1)
  df <- df %>% filter(Master_subject_ID %in% temp$Master_subject_ID)
  
  #Look at the time difference between scans
  temp <- df %>% select(Scan_ID, Master_subject_ID,Model,Cohort,Protocols, Days_since_first_scan)
  temp <- temp %>% 
    group_by(Master_subject_ID) %>%
    arrange(Days_since_first_scan, .by_group = T)%>%
    mutate(Days_since_previous_scan = Days_since_first_scan - lag(Days_since_first_scan)) 
  
  #rank the scans in chronological order
  anyscanner <- temp %>%
    group_by(Master_subject_ID) %>%
    mutate(Chron_rank = order(Days_since_first_scan, decreasing=FALSE))
  anyscanner$Rank_ID <- paste0(anyscanner$Master_subject_ID, "_", anyscanner$Chron_rank)
  anyscanner$Scanner <- paste0(anyscanner$Cohort, "_", anyscanner$Model)
  
  #Identify all scans that have been done within 6 months of the previous scan
  any_6m <- anyscanner %>% filter(Days_since_previous_scan <= 180)
  n_distinct(any_6m$Master_subject_ID)
  any_6m$Prev_scan <- paste0(any_6m$Master_subject_ID, "_", (any_6m$Chron_rank -1))
  model2 <- anyscanner %>% select(Rank_ID, Scanner)
  model2 <- model2 %>% rename(Scanner2 = Scanner)
  any_6m <- merge(any_6m, model2, by.x = "Prev_scan", by.y = "Rank_ID", all.x = TRUE, all.y = FALSE)
  any_6m$Master_subject_ID <- any_6m$Master_subject_ID.x
  any_6m <- any_6m %>% select(-c(Master_subject_ID.x, Master_subject_ID.y))
  dif <- any_6m %>% filter(Scanner != Scanner2)
  
  #define and format Scanner variable
  dif$Scanner_pair <- paste0(dif$Scanner, "_VS_", dif$Scanner2)
  dif$Scanner_pair <- 
    ifelse(dif$Scanner_pair %in% c("DETECT_Verio_VS_LEGACY_Trio", "LEGACY_Trio_VS_DETECT_Verio"), "LEGACY_Trio_VS_DETECT_Verio",
           ifelse(dif$Scanner_pair %in% c("LEGACY_Trio_VS_LEGACY_Verio", "LEGACY_Verio_VS_LEGACY_Trio"), "LEGACY_Verio_VS_LEGACY_Trio", 
                  as.character(dif$Scanner_pair)))
  trav <-  dif %>% 
    group_by(Master_subject_ID) %>% 
    slice(which.min(Days_since_previous_scan))
  
  #get Scan_ID for the previous scan and tidy dataframe
  temp <- anyscanner %>% select(Scan_ID2 = Scan_ID, Prev_scan = Rank_ID)
  trav <- merge(temp, trav, by = "Prev_scan", all.x = FALSE, all.y = TRUE)
  trav <- trav %>% select(Master_subject_ID = Master_subject_ID.x,
                          Scan_ID1 = Scan_ID,
                          Scan_ID2,
                          Scanner1 = Scanner,
                          Scanner2,
                          Scanner_pair,
                          Days = Days_since_previous_scan)
  
  return(trav)
  
}


########################################
#                                      #
# ADD VOLUME TO ACROSS-SCANNER COHORT  #
#                                      #
########################################

add_vol_across <- function(trav, vol){
  
  #add volume
  v1 <- vol %>% filter(Scan_ID %in% trav$Scan_ID1)
  colnames(v1)[-1] <- paste0("T1_", colnames(v1)[-1])
  
  v2 <- vol %>% filter(Scan_ID %in% trav$Scan_ID2)
  colnames(v2)[-1] <- paste0("T2_", colnames(v2)[-1])
  
  trav <- merge(trav, v1, by.x = "Scan_ID1", by.y = "Scan_ID")
  trav <- merge(trav, v2, by.x = "Scan_ID2", by.y = "Scan_ID")
  
  #make into long format
  t1 <- trav %>% select(Master_subject_ID, Scan_ID1, Scanner1, T1_Ventricles:T1_Brainstem)
  t1$Scan_num <- 1
  t1$Days <- 0
  t1 <- t1 %>% select(Master_subject_ID,
                      Scan_ID = Scan_ID1,
                      Scanner = Scanner1,
                      Scan_num,
                      Days,
                      Ventricles = T1_Ventricles,
                      SupraWM = T1_SupraWM,
                      SupraCortex = T1_SupraCortex,
                      SupraDGM = T1_SupraDGM,
                      CerebellumGM = T1_CerebellumGM,
                      CerebellumWM = T1_CerebellumWM,
                      Brainstem = T1_Brainstem)
  
  t2 <- trav %>% select(Master_subject_ID, Scan_ID2, Scanner2, Days, T2_Ventricles:T2_Brainstem)
  t2$Scan_num <- 2
  t2 <- t2 %>% select(Master_subject_ID,
                      Scan_ID = Scan_ID2,
                      Scanner = Scanner2,
                      Scan_num,
                      Days,
                      Ventricles = T2_Ventricles,
                      SupraWM = T2_SupraWM,
                      SupraCortex = T2_SupraCortex,
                      SupraDGM = T2_SupraDGM,
                      CerebellumGM = T2_CerebellumGM,
                      CerebellumWM = T2_CerebellumWM,
                      Brainstem = T2_Brainstem)
  
  tlong <- rbind(t1, t2) %>% as.data.frame()
  
  return(tlong)
  
}

##########################################
#                                        #
# ADD CORTHICK TO ACROSS-SCANNER COHORT  #
#                                        #
##########################################

add_corthick_across <- function(trav, corthick){
  
  #add corthick
  v1 <- corthick %>% filter(Scan_ID %in% trav$Scan_ID1)
  colnames(v1)[-1] <- paste0("T1_", colnames(v1)[-1])
  
  v2 <- corthick %>% filter(Scan_ID %in% trav$Scan_ID2)
  colnames(v2)[-1] <- paste0("T2_", colnames(v2)[-1])
  
  trav <- merge(trav, v1, by.x = "Scan_ID1", by.y = "Scan_ID")
  trav <- merge(trav, v2, by.x = "Scan_ID2", by.y = "Scan_ID")
  
  #make into long format
  t1 <- trav %>% select(Master_subject_ID, Scan_ID1, Scanner1, T1_Frontal:T1_WholeCortex)
  t1$Scan_num <- 1
  t1$Days <- 0
  t1 <- t1 %>% select(Master_subject_ID,
                      Scan_ID = Scan_ID1,
                      Scanner = Scanner1,
                      Scan_num,
                      Days,
                      Frontal = T1_Frontal,
                      Insular = T1_Insular,
                      Parietal = T1_Parietal,
                      Occipital = T1_Occipital,
                      Temporal = T1_Temporal,
                      Hippocampal = T1_Hippocampal,
                      WholeCortex = T1_WholeCortex)
  
  t2 <- trav %>% select(Master_subject_ID, Scan_ID2, Scanner2, Days, T2_Frontal:T2_WholeCortex)
  t2$Scan_num <- 2
  t2 <- t2 %>% select(Master_subject_ID,
                      Scan_ID = Scan_ID2,
                      Scanner = Scanner2,
                      Scan_num,
                      Days,
                      Frontal = T2_Frontal,
                      Insular = T2_Insular,
                      Parietal = T2_Parietal,
                      Occipital = T2_Occipital,
                      Temporal = T2_Temporal,
                      Hippocampal = T2_Hippocampal,
                      WholeCortex = T2_WholeCortex)
  
  tlong <- rbind(t1, t2) %>% as.data.frame()
  
  return(tlong)
  
}

##########################################
#                                        #
# GATHER ALL SCANS FOR HARMONISATION     #
#                                        #
##########################################


all_struc_scans <- function(metric, slong, tlong, clinical) {
  
  #combine scans from within- and across-scanner cohorts
  forcombat <- metric %>% filter(Scan_ID %in% slong$Scan_ID | Scan_ID %in% tlong$Scan_ID)
  forcombat <- forcombat %>% group_by(Scan_ID) %>% slice(1) #no overlap between with and across scanners
  
  forcombat$Within <- ifelse(forcombat$Scan_ID %in% slong$Scan_ID, "yes", "no")
  forcombat$Across <- ifelse(forcombat$Scan_ID %in% tlong$Scan_ID, "yes", "no")
  
  #add clinical data
  clinical <- clinical %>% select(Master_subject_ID, Sex)
  scans <- scans %>% select(Master_subject_ID, Scan_ID, Age_at_this_scan, Model, Cohort, Days_since_first_scan, ICV)
  scans <- merge(scans, clinical, by = "Master_subject_ID", all.x = TRUE, all.y = TRUE)
  forcombat <-  merge(forcombat, scans, by = "Scan_ID", all.x = TRUE, all.y = FALSE)
  forcombat <- forcombat %>% unite("Scanner", c(Cohort, Model))
  forcombat$Scanner <- as.factor(forcombat$Scanner)
  
  
  #redo Days_since_first_scan relative to the first scan in my cohort
  forcombat <- forcombat %>% 
    group_by(Master_subject_ID) %>% 
    mutate(Time = Days_since_first_scan - min(Days_since_first_scan))
  forcombat <- forcombat %>% select(-Days_since_first_scan) 
  
  #binarise Sex for matrices
  forcombat$Sex <- ifelse(forcombat$Sex == "male", 1, 0)
  
  return(forcombat)
}