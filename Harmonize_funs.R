#Harmonisation functions

##############################
#                            #
# Cross-sectional Combat     #
#                            #
##############################

get_crosscombat <- function(forcombat, from = "Ventricles", to = "Brainstem", is_dti = "no", prior = "parametric"){
  
  #we need a matrix with Scan_ID as column names and ROI index as row names
  dat <- forcombat %>% ungroup %>% select(from:to)
  if (is_dti == "no"){
    rownames(dat) <- forcombat$Scan_ID
  } else {
    rownames(dat) <- forcombat$Merge_ID
  }
  
  dat <- t(dat)
  
  #batch covariate i.e. scanner
  batch <- forcombat$Scanner
  
  #making the model matrix
  #this does NOT include batch
  mod <- forcombat
  mod$Intercept <- 1
  mod <- 
    mod %>% 
    ungroup %>% 
    select(Intercept, Age_at_this_scan, Sex, ICV) %>% 
    as.matrix()
  
  #use neuroCombat
  if (prior == "parametric"){
    result <- neuroCombat(dat = dat, 
                          batch = batch, 
                          mod = mod, 
                          verbose = FALSE, 
                          mean.only = FALSE)
  } else if (prior == "non-parametric") {
    result <- neuroCombat(dat = dat, 
                          batch = batch, 
                          mod = mod, 
                          verbose = FALSE, 
                          mean.only = FALSE,
                          parametric = FALSE)
  } else if (prior == "none") {
    result <- neuroCombat(dat = dat, 
                         batch = batch, 
                         mod = mod, 
                         verbose = FALSE, 
                         mean.only = FALSE,
                         eb = FALSE)
  } else {
    print("Please select valid prior")
  }
  

  
  if (is_dti == "no"){
    crosscombat <- cbind(forcombat$Scan_ID,
                         as.character(forcombat$Scanner), 
                         forcombat$Master_subject_ID, 
                         t(result[[1]])) %>% 
      as.data.frame()
    colnames(crosscombat)[1:3] <- c("Scan_ID", "Scanner", "Master_subject_ID")
  } else {
    crosscombat <- cbind(forcombat$Merge_ID,
                         as.character(forcombat$Scanner), 
                         forcombat$Master_subject_ID, 
                         t(result[[1]])) %>% 
      as.data.frame()
    colnames(crosscombat)[1:3] <- c("Merge_ID", "Scanner", "Master_subject_ID")
  }
  
  return(crosscombat)
}


##############################
#                            #
# gamCombat / neuroHarmonize #
#                            #
##############################

get_gamcombat <- function(forcombat, from = "Ventricles", to = "Brainstem", is_dti = "no"){
  
  #we need a matrix with Scan_ID as row names and ROI index as column names
  dat <- forcombat %>% ungroup %>% select(from:to)
  if (is_dti == "no"){
    rownames(dat) <- forcombat$Scan_ID
  } else {
    rownames(dat) <- forcombat$Merge_ID
  }
  dat <- as.matrix(dat)
  
  
  #making the model matrix
  #this DOES include batch aka SITE
  covars <- forcombat
  covars <- 
    covars %>% 
    ungroup %>% 
    select(SITE = Scanner, Age_at_this_scan, Sex, ICV) 
  
  #use neuroHarmonize aka gamCombat
  #assume Age to be non-linear
  result <- py$harmonizationLearn(dat, covars, smooth_terms = "Age_at_this_scan")
  colnames(result[[2]]) <- colnames(dat)
  
  if (is_dti == "no"){
    gamcombat <- cbind(forcombat$Scan_ID,
                         as.character(forcombat$Scanner), 
                         forcombat$Master_subject_ID, 
                         result[[2]]) %>% #the second list element contains the harmonized data
      as.data.frame()
    colnames(gamcombat)[1:3] <- c("Scan_ID", "Scanner", "Master_subject_ID")
  } else {
    gamcombat <- cbind(forcombat$Merge_ID,
                         as.character(forcombat$Scanner), 
                         forcombat$Master_subject_ID, 
                         result[[2]]) %>% #the second list element contains the harmonized data
      as.data.frame()
    colnames(gamcombat)[1:3] <- c("Merge_ID", "Scanner", "Master_subject_ID")
  }
  
  return(gamcombat)
}

##############################################
#                                            #
# Longitudinal Combat - intercept only       #
#                                            #
##############################################

get_longcombat <- function(forcombat, is_corthick = "no", is_dti = "no"){
  
  
  rois1 <- c("Ventricles", 
             "SupraCortex", 
             "SupraWM", 
             "SupraDGM", 
             "CerebellumGM",
             "CerebellumWM", 
             "Brainstem")
  
  rois2 <- c("Frontal", 
             "Parietal", 
             "Insular", 
             "Occipital", 
             "Temporal", 
             "Hippocampal",
             "WholeCortex")
  
  if (is_corthick == "no"){
    rois <- rois1
  } else {
    rois <- rois2
  }
  
  df <- forcombat %>% select(Master_subject_ID,Time, Scanner, Age_at_this_scan, Sex, ICV, rois)
  batch <<- as.factor(df$Scanner)
  
  longres <- longCombat(
    idvar = "Master_subject_ID",
    timevar = "Time",
    batchvar = "Scanner",
    features = rois,
    formula = "Age_at_this_scan + Sex + Time + ICV",
    ranef = "(1 |Master_subject_ID)",
    data = df,
    niter = 30,
    method = 'REML',
    verbose = FALSE
  )
  
  if (is_dti == "no"){
    longcombat <- cbind(forcombat$Scan_ID, longres[[1]]) %>% as.data.frame()
    colnames(longcombat)[1] <- "Scan_ID"
  } else {
    longcombat <- cbind(forcombat$Merge_ID, longres[[1]]) %>% as.data.frame()
    colnames(longcombat)[1] <- "Merge_ID"
  }
  
  colnames(longcombat) <- gsub(".combat", "", colnames(longcombat))
  
  return(longcombat)
}

##############################################
#                                            #
# Longitudinal Combat - intercept & slope    #
#                                            #
##############################################

get_longcombat_plusslope <- function(forcombat, is_corthick = "no", is_dti = "no"){
  
  
  rois1 <- c("Ventricles", 
             "SupraCortex", 
             "SupraWM", 
             "SupraDGM", 
             "CerebellumGM",
             "CerebellumWM", 
             "Brainstem")
  
  rois2 <- c("Frontal", 
             "Parietal", 
             "Insular", 
             "Occipital", 
             "Temporal", 
             "Hippocampal",
             "WholeCortex")
  
  if (is_corthick == "no"){
    rois <- rois1
  } else {
    rois <- rois2
  }
  
  df <- forcombat %>% select(Master_subject_ID,Time, Scanner, Age_at_this_scan, Sex, ICV, rois)
  batch <<- as.factor(df$Scanner)
  
  longres <- longCombat(
    idvar = "Master_subject_ID",
    timevar = "Time",
    batchvar = "Scanner",
    features = rois,
    formula = "Age_at_this_scan + Sex + Time + ICV",
    ranef = "(1 + Time |Master_subject_ID)",
    data = df,
    niter = 30,
    method = 'REML',
    verbose = FALSE
  )
  
  if (is_dti == "no"){
    longcombat <- cbind(forcombat$Scan_ID, longres[[1]]) %>% as.data.frame()
    colnames(longcombat)[1] <- "Scan_ID"
  } else {
    longcombat <- cbind(forcombat$Merge_ID, longres[[1]]) %>% as.data.frame()
    colnames(longcombat)[1] <- "Merge_ID"
  }
  
  colnames(longcombat) <- gsub(".combat", "", colnames(longcombat))
  
  return(longcombat)
}


#########################################################
#                                                       #
# Extract harmonized data for across-scanner cohort     #
#                                                       #
#########################################################


extract_across <- function(harmonised_all, tlong, is_corthick = "no") {
  
  #pick scans for across scanner comparison and add scan number
  harmonised_across <- harmonised_all %>% filter(Scan_ID %in% tlong$Scan_ID)
  num               <- tlong %>% select(Scan_ID, Scan_num)
  harmonised_across <- merge(num, harmonised_across)
  harmonised_across <- harmonised_across[order(harmonised_across$Scan_num),]
  
  if (is_corthick == "no"){
    
    #calculate CoV and make short dataframe
    harmonised_across_short <- calculate_CoV(harmonised_across)
    
  } else {
    
    harmonised_across_short <- calculate_CoV_corthick(harmonised_across)
  }
  
  return(harmonised_across_short)
}
