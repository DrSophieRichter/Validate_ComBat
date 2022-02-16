#Functions to select scans with DTI data


############################################
#                                          #
# Extract md or fa with DTI series name    #
# as a subset of those with structural data#
#                                          #
############################################

extract_dti <- function(dti, metric = "md", invol){
  
  #to avoid running into numerical problems I will multiply all values 
  #by 100,000 for md or
  #by 100 for fa
  
  if (metric == "md"){
    dti <- dti %>% filter(Metric == "md")
    dti <- dti %>% 
      mutate_at(vars(Ventricles:Brainstem),
                .funs = funs(. * 100000))
  } else {
    dti <- dti %>% filter(Metric == "fa")
    dti <- dti %>% 
      mutate_at(vars(Ventricles:Brainstem),
                .funs = funs(. * 100))
  }
  
  #rationalise series label
  dti_id <- dti %>% select(Scan_ID, DTI_ID)
  gsub_v <- Vectorize(gsub, c("pattern", "x"))
  dti_id$Series <- gsub_v(paste0(dti_id$Scan_ID,"-"), "", dti_id$DTI_ID)
  
  dti_id$Series <- ifelse(is.na(dti_id$Series)==TRUE, NA,
                          ifelse(str_detect(dti_id$Series, "HighRes")==TRUE, "HighRes",
                                 ifelse(str_detect(dti_id$Series, "FREE_12")==TRUE, "FREE_12", 
                                        ifelse(str_detect(dti_id$Series, "FREE_63")==TRUE, "FREE_63", dti_id$Series))))
  
  
  #Remove HighRes
  highres <- dti_id %>% filter(Series == "HighRes")
  dti <- dti %>% filter(!DTI_ID %in% highres$DTI_ID)

  #Add rationalised series label to the dataframe
  dti <- dti %>% select(Scan_ID,  DTI_ID,
                      Ventricles, 
                      SupraCortex, 
                      SupraWM, 
                      SupraDGM,
                      CerebellumGM, 
                      CerebellumWM,
                      Brainstem)
  dti <- merge(dti_id[-1], dti, by = "DTI_ID", all.x = FALSE, all.y = TRUE)
  dti$Merge_ID <- paste0(dti$Scan_ID, "_", dti$Series)
  
  #Remove scans not present in volume dataset for consistency of presentation
  dti <- dti %>% filter(Scan_ID %in% invol$Scan_ID)
  
  return(dti)
}

#########################################################
#                                                       #
# Extract md or fa with DTI series name                 #
# irrespective of whether they are in structural dataset#
#                                                       #
#########################################################

extract_dti_fresh <- function(dti, metric = "md"){
  
  #to avoid running into numerical problems I will multiply all values 
  #by 100,000 for md or
  #by 100 for fa
  
  if (metric == "md"){
    dti <- dti %>% filter(Metric == "md")
    dti <- dti %>% 
      mutate_at(vars(Ventricles:Brainstem),
                .funs = funs(. * 100000))
  } else {
    dti <- dti %>% filter(Metric == "fa")
    dti <- dti %>% 
      mutate_at(vars(Ventricles:Brainstem),
                .funs = funs(. * 100))
  }
  
  #rationalise series label
  dti_id <- dti %>% select(Scan_ID, DTI_ID)
  gsub_v <- Vectorize(gsub, c("pattern", "x"))
  dti_id$Series <- gsub_v(paste0(dti_id$Scan_ID,"-"), "", dti_id$DTI_ID)
  
  dti_id$Series <- ifelse(is.na(dti_id$Series)==TRUE, NA,
                          ifelse(str_detect(dti_id$Series, "HighRes")==TRUE, "HighRes",
                                 ifelse(str_detect(dti_id$Series, "FREE_12")==TRUE, "FREE_12", 
                                        ifelse(str_detect(dti_id$Series, "FREE_63")==TRUE, "FREE_63", dti_id$Series))))
  
  
  #Remove HighRes
  highres <- dti_id %>% filter(Series == "HighRes")
  dti <- dti %>% filter(!DTI_ID %in% highres$DTI_ID)
  
  #Add rationalised series label to the dataframe
  dti <- dti %>% select(Scan_ID,  DTI_ID,
                        Ventricles, 
                        SupraCortex, 
                        SupraWM, 
                        SupraDGM,
                        CerebellumGM, 
                        CerebellumWM,
                        Brainstem)
  dti <- merge(dti_id[-1], dti, by = "DTI_ID", all.x = FALSE, all.y = TRUE)
  dti$Merge_ID <- paste0(dti$Scan_ID, "_", dti$Series)
  
  return(dti)
}

##########################################
#                                        #
# SELECT WITHIN SCANNER COHORT           #
#                                        #
##########################################

select_dti_within <- function(md, scans, qc_bad){
  
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
  
  #Select patients scanned twice within 6 months
  temp <- df %>% 
    group_by(Master_subject_ID, Scanner) %>% summarise(num = n()) %>% filter(num >1)
  temp$Series_ID <- paste0(temp$Master_subject_ID, "_", temp$Scanner)
  df$Series_ID <- paste0(df$Master_subject_ID, "_", df$Scanner)
  rescan_ever <- df %>% filter(Series_ID %in% temp$Series_ID)
  
  #Time between scans of the same series
  rescan_ever <- rescan_ever %>% group_by(Series_ID) %>%
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
  
  #keep only one scan pair per subject
  keep <-  rescan_6m %>% group_by(Master_subject_ID) %>% slice(which.min(Days_since_previous_scan))
  
  #pick out the dti data for the identified serial scans
  s1 <- md %>% filter(Scan_ID %in% keep$Scan_ID)
  colnames(s1)[-2] <- paste0("S1_", colnames(s1)[-2]) 
  
  #identify subjects with more than one diffusion acquired per session
  temp <- s1 %>% group_by(S1_Scan_ID) %>% summarise(num = n()) %>% filter(num >1)
  s1 %>% filter(S1_Scan_ID %in% temp$S1_Scan_ID) %>% select(S1_Scan_ID, Series)
  #I will just keep the top diffusion for each session
  s1 <- s1 %>% group_by(S1_Scan_ID) %>% slice(1)
  
  temp <- rescan_ever %>% select(Scan_ID2 = Scan_ID, Rank_ID2 = Rank_ID, Series_ID2 = Series_ID)
  keep <- merge(temp, keep, by.x = "Rank_ID2", by.y = "Prev_scan", all.x = FALSE, all.y = TRUE)
  s2 <- md %>% filter(Scan_ID %in% keep$Scan_ID2)
  colnames(s2)[-2] <- paste0("S2_", colnames(s2)[-2]) 
  
  #identify subjects with more than one diffusion acquired per session
  temp <- s2 %>% group_by(S2_Scan_ID) %>% summarise(num = n()) %>% filter(num >1)
  s2 %>% filter(S2_Scan_ID %in% temp$S2_Scan_ID) %>%  select(S2_Scan_ID, Series)
  #I will just keep the top diffusion for each session
  s2 <- s2 %>% group_by(S2_Scan_ID) %>% slice(1)
  
  #Now bring them together
  keep <- keep %>% filter(Scan_ID %in% s1$S1_Scan_ID & Scan_ID2 %in% s2$S2_Scan_ID)
  keep <- merge(keep, s1, by.x = "Scan_ID", by.y = "S1_Scan_ID", all.x = TRUE, all.y = FALSE)
  keep <- merge(keep, s2, by.x = "Scan_ID2", by.y = "S2_Scan_ID", all.x = TRUE, all.y = FALSE)
  keep <- keep %>% rename(Series1 = Series.x,
                          Series2 = Series.y)
  
  #Ensuring that the repeat scans for each subject have been done using the same sequence
  keep <- keep %>% select(Series1, Series2, Scan_ID, Scan_ID2, everything())
  keep <- keep %>% filter(Series1 == Series2)
  
  #Make into long format AND add dti sequence
  s1 <- keep %>% select(Master_subject_ID, Scan_ID, Scanner, S1_Ventricles:S1_Brainstem)
  s1$Days <- 0
  s1$Scan_num <- 1
  s1 <- s1 %>% select(Master_subject_ID, Scanner, Scan_ID, Scan_num, Days, everything())
  
  s2 <- keep %>% select(Master_subject_ID, Scanner, Scan_ID = Scan_ID2, Days = Days_since_previous_scan, S2_Ventricles:S2_Brainstem)
  s2$Scan_num <- 2
  s2 <- s2 %>% select(Master_subject_ID, Scanner, Scan_ID, Scan_num, Days, everything())
  
  #long
  colnames(s1)[6:12] <- gsub("S1_","", colnames(s1)[6:12])
  colnames(s2)[6:12] <- gsub("S2_","", colnames(s2)[6:12])
  slong <- rbind(s1, s2) %>% as.data.frame()
  
  #Make a merge_ID or unique identifier per scan
  slong <- slong %>% separate(Scanner, remove = FALSE, into = c("Model", "Series"), extra = "merge", sep = "_")
  slong$Merge_ID <- paste0(slong$Scan_ID, "_", slong$Series)
  
  return(slong)
}


##########################################
#                                        #
# SELECT ACROSS-SCANNER COHORT           #
#                                        #
##########################################

select_dti_across <- function(md, scans, qc_bad){
  
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
  
  #select only subjects with at least 2 scans
  temp <- df %>% group_by(Master_subject_ID) %>% summarise(num = n()) %>% filter(num>1)
  df <- df %>% filter(Master_subject_ID %in% temp$Master_subject_ID)
  
  #Look at the time difference between scans
  temp <- df %>% select(Scan_ID, Master_subject_ID,Model,Cohort,Series, Days_since_first_scan)
  temp <- temp %>% group_by(Master_subject_ID) %>%
    mutate(Days_since_previous_scan = Days_since_first_scan - lag(Days_since_first_scan)) 
  
  #rank the scans in chronological order
  anyscanner <- temp %>%
    group_by(Master_subject_ID) %>%
    mutate(Chron_rank = order(Days_since_first_scan, decreasing=FALSE))
  anyscanner$Rank_ID <- paste0(anyscanner$Master_subject_ID, "_", anyscanner$Chron_rank)
  anyscanner$Scanner <- paste0(anyscanner$Model, "_", anyscanner$Series)
  
  #Identify all scans that have been done within 6 months of the previous scan
  any_6m <- anyscanner %>% filter(Days_since_previous_scan <= 180)
  n_distinct(any_6m$Master_subject_ID)
  any_6m$Prev_scan <- paste0(any_6m$Master_subject_ID, "_", (any_6m$Chron_rank -1))
  model2 <- anyscanner %>% select(Rank_ID, Scanner, Scan_ID)
  model2 <- model2 %>% rename(Scanner2 = Scanner, Scan_ID2 = Scan_ID)
  any_6m <- merge(any_6m, model2, by.x = "Prev_scan", by.y = "Rank_ID", all.x = TRUE, all.y = FALSE)
  any_6m$Master_subject_ID <- any_6m$Master_subject_ID.x
  any_6m <- any_6m %>% select(-c(Master_subject_ID.x, Master_subject_ID.y))
  dif <- any_6m %>% filter(Scanner != Scanner2)
  
  #Make two cohorts
  #a) same session (Scan_ID) different sequence
  #b) separate session (usually seperate days) with different sequence and/or model
  
  difa <- dif %>% filter(Scan_ID == Scan_ID2)
  difb <- dif %>% filter(Scan_ID != Scan_ID2)
  
  difb$Scanner_pair <- paste0(difb$Scanner, " vs ", difb$Scanner2)
  difb$Scanner_pair <- ifelse(difb$Scanner_pair == "Trio_FREE_63_VS_Trio_FREE_12", "Trio_FREE_12_VS_Trio_FREE_63",
                              ifelse(difb$Scanner_pair == "Trio_FREE_12_VS_Verio_FREE_63", "Verio_FREE_63_VS_Trio_FREE_12", as.character(difb$Scanner_pair)))
  
  #pick only one scan pair per subject
  trav <-  difb %>% group_by(Master_subject_ID) %>% slice(which.min(Days_since_previous_scan))
  
  trav <- trav %>% select(Master_subject_ID,
                          Scan_ID1 = Scan_ID,
                          Scan_ID2,
                          Scanner1 = Scanner,
                          Scanner2,
                          Scanner_pair,
                          Days = Days_since_previous_scan)
  
  #add md (or fa)
  #Add a merge ID to both dataframes (md and dataframe trav) consisting of Scan_ID and series name
  trav <- trav %>% separate(Scanner1, remove = FALSE, into = c("Model1", "Series1"), extra = "merge", sep = "_")
  trav$Merge_ID1 <- paste0(trav$Scan_ID1, "_", trav$Series1)
  trav <- trav %>% separate(Scanner2, remove = FALSE, into = c("Model2", "Series2"), extra = "merge", sep = "_")
  trav$Merge_ID2 <- paste0(trav$Scan_ID2, "_", trav$Series2)
  
  md1 <- md %>% filter(Merge_ID %in% trav$Merge_ID1)
  md1 <- md1 %>% group_by(Merge_ID)%>% slice(1)
  colnames(md1)[-2] <- paste0("T1_", colnames(md1)[-2])
  
  md2 <- md %>% filter(Merge_ID %in% trav$Merge_ID2)
  md2 <- md2 %>% group_by(Merge_ID)%>% slice(1)
  colnames(md2)[-2] <- paste0("T2_", colnames(md2)[-2])
  
  trav <- merge(trav, md1, by.x = "Merge_ID1", by.y = "T1_Merge_ID", all.x = TRUE, all.y = FALSE)
  trav <- merge(trav, md2, by.x = "Merge_ID2", by.y = "T2_Merge_ID", all.x = TRUE, all.y = FALSE)
  trav <- trav %>% select(-c(Series.x, Series.y))
  
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
  
  #make a merge_ID or unique identifier per scan
  tlong <- tlong %>% separate(Scanner, remove = FALSE, into = c("Model", "Series"), extra = "merge", sep = "_")
  tlong$Merge_ID <- paste0(tlong$Scan_ID, "_", tlong$Series)
  
  return(tlong)
}

##########################################
#                                        #
# GATHER ALL SCANS FOR HARMONISATION     #
#                                        #
##########################################


all_dti_scans <- function(metric, slong, tlong, clinical) {
  
  #combine scans from within- and across-scanner cohorts
  forcombat <- metric %>% filter(Merge_ID %in% slong$Merge_ID | Merge_ID %in% tlong$Merge_ID)
  forcombat <- forcombat %>% group_by(Merge_ID) %>% slice(1) #4 scans appear in both within and across scanner
  forcombat <- forcombat %>% select(Merge_ID, everything())
  
  forcombat$Within_fa <- ifelse(forcombat$Merge_ID %in% slong$Merge_ID, "yes", "no")
  forcombat$Across_fa <- ifelse(forcombat$Merge_ID %in% tlong$Merge_ID, "yes", "no")
  
  #add clinical data
  clinical <- clinical %>% select(Master_subject_ID, Sex)
  scans <- scans %>% select(Master_subject_ID, Scan_ID, Age_at_this_scan, Model, Cohort, Days_since_first_scan)
  scans <- merge(scans, clinical, by = "Master_subject_ID", all.x = TRUE, all.y = TRUE)
  forcombat <-  merge(forcombat, scans, by = "Scan_ID", all.x = TRUE, all.y = FALSE)
  forcombat <- forcombat %>% unite("Scanner", c(Model, Series))
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
