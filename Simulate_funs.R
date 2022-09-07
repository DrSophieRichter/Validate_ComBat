#Functions to simulate comparisons between group A and B

##################################################
#                                                #
# Get p-value for Group from mixed model summary #
#                                                #
##################################################

get_p <- function(fit){
  coef <- summary(fit)$coefficients
  p <- coef[which(rownames(coef) == "Group"),which(colnames(coef)=="Pr(>|t|)")]
  return(p)
}

#######################################################
#                                                     #
# Get p-value for Group:Time from mixed model summary #
#                                                     #
#######################################################

get_p_slope <- function(fit){
  coef <- summary(fit)$coefficients
  p <- coef[which(rownames(coef) == "Group:Time"),which(colnames(coef)=="Pr(>|t|)")]
  return(p)
}

###################################################
#                                                 #
# Compare Groups A & B for each ROI (not corthick)#
# for difference in intercept                     #
#                                                 #
###################################################

#Note I do not allow a random slope per subject here
#because the time is almost zero
#so the slope is driven by random noise
#computing it is not meaningful
#and interferes with model convergence
#I do not include the interaction term as I am interested in the main effect

test_dif <- function(current, res){
  #test for group difference and save p-value
  fit <- lmerTest::lmer(Ventricles ~ 
                          Age_at_this_scan + Sex + Group*Time +
                          (1| Master_subject_ID), 
                        data = current)
  res[i,1] <- get_p(fit)
  
  fit <- lmerTest::lmer(SupraCortex ~ 
                          Age_at_this_scan + Sex + Group*Time +
                          (1| Master_subject_ID), 
                        data = current)
  res[i,2] <- get_p(fit)
  
  fit <- lmerTest::lmer(SupraWM ~ 
                          Age_at_this_scan + Sex + Group*Time +
                          (1| Master_subject_ID), 
                        data = current)
  res[i,3] <- get_p(fit)
  
  fit <- lmerTest::lmer(SupraDGM ~ 
                          Age_at_this_scan + Sex + Group*Time +
                          (1| Master_subject_ID), 
                        data = current)
  res[i,4] <- get_p(fit)
  
  fit <- lmerTest::lmer(CerebellumGM ~ 
                          Age_at_this_scan + Sex + Group*Time +
                          (1| Master_subject_ID), 
                        data = current)
  res[i,5] <- get_p(fit)
  
  fit <- lmerTest::lmer(CerebellumWM ~ 
                          Age_at_this_scan + Sex + Group*Time +
                          (1| Master_subject_ID), 
                        data = current)
  res[i,6] <- get_p(fit)
  
  fit <- lmerTest::lmer(Brainstem ~ 
                          Age_at_this_scan + Sex + Group*Time +
                          (1| Master_subject_ID), 
                        data = current)
  res[i,7] <- get_p(fit)
  
  return(res)
}

###################################################
#                                                 #
# Compare Groups A & B for each ROI (not corthick)#
# for difference in slope                         #
#                                                 #
###################################################


test_dif_slope <- function(current, res){
  #test for group difference and save p-value
  fit <- lmerTest::lmer(Ventricles ~ 
                          Age_at_this_scan + Sex + Group*Time + 
                          (1+Time| Master_subject_ID), 
                        data = current)
  res[i,1] <- get_p_slope(fit)
  
  fit <- lmerTest::lmer(SupraCortex ~ 
                          Age_at_this_scan + Sex + Group*Time + 
                          (1+Time| Master_subject_ID), 
                        data = current)
  res[i,2] <- get_p_slope(fit)
  
  fit <- lmerTest::lmer(SupraWM ~ 
                          Age_at_this_scan + Sex + Group*Time + 
                          (1+Time| Master_subject_ID), 
                        data = current)
  res[i,3] <- get_p_slope(fit)
  
  fit <- lmerTest::lmer(SupraDGM ~ 
                          Age_at_this_scan + Sex + Group*Time + 
                          (1+Time| Master_subject_ID), 
                        data = current)
  res[i,4] <- get_p_slope(fit)
  
  fit <- lmerTest::lmer(CerebellumGM ~ 
                          Age_at_this_scan + Sex + Group*Time + 
                          (1+Time| Master_subject_ID), 
                        data = current)
  res[i,5] <- get_p_slope(fit)
  
  fit <- lmerTest::lmer(CerebellumWM ~ 
                          Age_at_this_scan + Sex + Group*Time + 
                          (1+Time| Master_subject_ID), 
                        data = current)
  res[i,6] <- get_p_slope(fit)
  
  fit <- lmerTest::lmer(Brainstem ~ 
                          Age_at_this_scan + Sex + Group*Time + 
                          (1+Time| Master_subject_ID), 
                        data = current)
  res[i,7] <- get_p_slope(fit)
  
  return(res)
}

###################################################
#                                                 #
# Compare Groups A & B for each ROI for CORTHICK  #
# for difference in intercept                     #
#                                                 #
###################################################


test_dif_corthick <- function(current, res){
  #test for group difference and save p-value
  fit <- lmerTest::lmer(Frontal ~ 
                          Age_at_this_scan + Sex + Group*Time +
                          (1| Master_subject_ID), 
                        data = current)
  res[i,1] <- get_p(fit)
  
  fit <- lmerTest::lmer(Parietal ~ 
                          Age_at_this_scan + Sex + Group*Time +
                          (1| Master_subject_ID), 
                        data = current)
  res[i,2] <- get_p(fit)
  
  fit <- lmerTest::lmer(Insular ~ 
                          Age_at_this_scan + Sex + Group*Time +
                          (1| Master_subject_ID), 
                        data = current)
  res[i,3] <- get_p(fit)
  
  fit <- lmerTest::lmer(Occipital ~ 
                          Age_at_this_scan + Sex + Group*Time +
                          (1| Master_subject_ID), 
                        data = current)
  res[i,4] <- get_p(fit)
  
  fit <- lmerTest::lmer(Temporal ~ 
                          Age_at_this_scan + Sex + Group*Time +
                          (1| Master_subject_ID), 
                        data = current)
  res[i,5] <- get_p(fit)
  
  fit <- lmerTest::lmer(Hippocampal ~ 
                          Age_at_this_scan + Sex + Group*Time +
                          (1| Master_subject_ID), 
                        data = current)
  res[i,6] <- get_p(fit)
  
  fit <- lmerTest::lmer(WholeCortex ~ 
                          Age_at_this_scan + Sex + Group*Time +
                          (1| Master_subject_ID), 
                        data = current)
  res[i,7] <- get_p(fit)
  
  return(res)
}

###################################################
#                                                 #
# Compare Groups A & B for each ROI for CORTHICK  #
# for difference in intercept                     #
#                                                 #
###################################################


test_dif_slope_corthick <- function(current, res){
  #test for group difference and save p-value
  fit <- lmerTest::lmer(Frontal ~ 
                          Age_at_this_scan + Sex + Group*Time + 
                          (1+Time| Master_subject_ID), 
                        data = current)
  res[i,1] <- get_p_slope(fit)
  
  fit <- lmerTest::lmer(Parietal ~ Age_at_this_scan + Sex + Group*Time + 
                          (1+Time| Master_subject_ID), 
                        data = current)
  res[i,2] <- get_p_slope(fit)
  
  fit <- lmerTest::lmer(Insular ~ 
                          Age_at_this_scan + Sex + Group*Time + 
                          (1+Time| Master_subject_ID), 
                        data = current)
  res[i,3] <- get_p_slope(fit)
  
  fit <- lmerTest::lmer(Occipital ~ 
                          Age_at_this_scan + Sex + Group*Time + 
                          (1+Time| Master_subject_ID), 
                        data = current)
  res[i,4] <- get_p_slope(fit)
  
  fit <- lmerTest::lmer(Temporal ~ 
                          Age_at_this_scan + Sex + Group*Time + 
                          (1+Time| Master_subject_ID), 
                        data = current)
  res[i,5] <- get_p_slope(fit)
  
  fit <- lmerTest::lmer(Hippocampal ~ 
                          Age_at_this_scan + Sex + Group*Time + 
                          (1+Time| Master_subject_ID), 
                        data = current)
  res[i,6] <- get_p_slope(fit)
  
  fit <- lmerTest::lmer(WholeCortex ~ 
                          Age_at_this_scan + Sex + Group*Time + 
                          (1+Time| Master_subject_ID), 
                        data = current)
  res[i,7] <- get_p_slope(fit)
  
  return(res)
}

##########################################################
#                                                        #
# Simulate FPR using unharmonized data                   #
#                                                        #
##########################################################

sim_unharm <- function(mymax, forcombat, is_corthick = "no", intercept_only = "yes"){
  
  #prepare dataframe to collect results (aka p-values)in
  runs <- seq(1, mymax)
  res <- as.data.frame(matrix(,ncol=7,nrow=length(runs)))
  
  if (is_corthick == "no"){
    colnames(res) <- c("Ventricles", 
                       "SupraCortex", 
                       "SupraWM", 
                       "SupraDGM",
                       "CerebellumGM",
                       "CerebellumWM",  
                       "Brainstem")
    
  } else {
    colnames(res) <- c("Frontal", 
                       "Parietal", 
                       "Insular", 
                       "Occipital",  
                       "Temporal", 
                       "Hippocampal",
                       "WholeCortex")
  }
  
  sim <- forcombat
  sub <- distinct(sim, Master_subject_ID)
  sub$Group <- NA
  n <- nrow(sub)
  
  for (i in runs){
    #generate random group assignement
    i <<- i
    set.seed(i)
    sub$Group <- sample(c(0,1), replace=TRUE, size=n)
    current <- merge(sim, sub, by = "Master_subject_ID", all.x = TRUE, all.y = FALSE)
    
    if (is_corthick == "no" & intercept_only == "yes"){
      res <- test_dif(current, res)
    }else if (is_corthick == "yes" & intercept_only == "yes"){
      res <- test_dif_corthick(current, res)  
    }else if (is_corthick == "no" & intercept_only == "no"){
      res <- test_dif_slope(current, res) 
    }else {
      res <- test_dif_slope_corthick(current, res) 
    }
  }
  
  fp_raw <- apply(res,2,function(x) sum(x < 0.05)) %>% as.data.frame()
  colnames(fp_raw) <- "FPR"
  fp_raw$Method <- "unharm"
  fp_raw$ROI <- row.names(fp_raw)
  
  return(fp_raw)
  
}


#############################################
#                                           #
# Simulate FPR after cross-sectional Combat #            
#                                           #
#############################################

sim_crosscombat <- function(mymax, forcombat, from = "Ventricles", to = "Brainstem", is_corthick = "no", is_fudti = "no", intercept_only = "yes"){
  
  #we need a matrix with Scan_ID (or, in the case of fu_dti data "Merge_ID") 
  #as column names and ROI index as row names
  dat <- forcombat %>% ungroup %>% select(from:to)
  if (is_fudti == "no") {
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
  sub <- distinct(mod, Master_subject_ID)
  sub$Group <- NA
  n <- nrow(sub)
  
  #prep table to collect results (aka simulated p-values) in
  runs <- seq(1, mymax)
  res  <- as.data.frame(matrix(,ncol=7,nrow=length(runs)))
  
  if (is_corthick == "no"){
    colnames(res) <- c("Ventricles", 
                       "SupraCortex", 
                       "SupraWM", 
                       "SupraDGM",
                       "CerebellumGM",  
                       "CerebellumWM", 
                       "Brainstem")
    
  } else {
    colnames(res) <- c("Frontal", 
                       "Parietal", 
                       "Insular", 
                       "Occipital",  
                       "Temporal", 
                       "Hippocampal",
                       "WholeCortex")
  }
  
  #use neuroCombat
  for (i in runs){
    #generate random group assignment
    i <<- i
    set.seed(i)
    sub$Group <- sample(c(0,1), replace=TRUE, size=n)
    current <- merge(mod, sub, by = "Master_subject_ID", all.x = TRUE, all.y = FALSE)
    current <- current %>% ungroup %>% select(Intercept, Age_at_this_scan, Sex, ICV, Group) %>% as.matrix()
    
    
    #use neuroCombat
    result <- neuroCombat(dat = dat, batch = batch, mod = current, verbose = FALSE)
    current <- cbind(forcombat$Master_subject_ID, forcombat$Time, current, t(result[[1]])) %>% as.data.frame()
    colnames(current)[1:2] <- c("Master_subject_ID", "Time")
    current[,-1] <- sapply(current[,-1], as.numeric)
    
    #test for group difference and save p-value
    if (is_corthick == "no" & intercept_only == "yes"){
      res <- test_dif(current, res)
    }else if (is_corthick == "yes" & intercept_only == "yes"){
      res <- test_dif_corthick(current, res)  
    }else if (is_corthick == "no" & intercept_only == "no"){
      res <- test_dif_slope(current, res) 
    }else {
      res <- test_dif_slope_corthick(current, res) 
    }
    
  }
  
  fp_xharm <- apply(res, 2 ,function(x) sum(x < 0.05)) %>% as.data.frame()
  colnames(fp_xharm) <- "FPR"
  fp_xharm$Method <- "crosscombat"
  fp_xharm$ROI <- row.names(fp_xharm)
  
  return(fp_xharm)
}



#############################################
#                                           #
# Simulate FPR after Combat-GAM             #            
#                                           #
#############################################

sim_gamcombat <- function(mymax, forcombat, from = "Ventricles", to = "Brainstem", is_corthick = "no", is_fudti = "no", intercept_only = "yes"){
  
  #we need a matrix with Scan_ID (or, in the cas eof fu_dti data "Merge_ID") 
  #as row names and ROI index as column names
  dat <- forcombat %>% ungroup %>% select(from:to)
  if (is_fudti == "no") {
    rownames(dat) <- forcombat$Scan_ID
  } else {
    rownames(dat) <- forcombat$Merge_ID
  }
  dat <- as.matrix(dat)
  
  #making covars
  #this does include batch
  covars    <- forcombat
  sub       <- distinct(covars, Master_subject_ID)
  sub$Group <- NA
  n         <- nrow(sub)
  
  #prep table to collect results (aka simulated p-values) in
  runs <- seq(1, mymax)
  res  <- as.data.frame(matrix(,ncol=7,nrow=length(runs)))
  
  if (is_corthick == "no"){
    colnames(res) <- c("Ventricles", 
                       "SupraCortex", 
                       "SupraWM", 
                       "SupraDGM",
                       "CerebellumGM",  
                       "CerebellumWM", 
                       "Brainstem")
    
  } else {
    colnames(res) <- c("Frontal", 
                       "Parietal", 
                       "Insular", 
                       "Occipital",  
                       "Temporal", 
                       "Hippocampal",
                       "WholeCortex")
  }
  
  #use neuroCombat
  for (i in runs){
    #generate random group assignment
    i <<- i
    set.seed(i)
    sub$Group <- sample(c(0,1), replace=TRUE, size=n)
    current   <- merge(covars, sub, by = "Master_subject_ID", all.x = TRUE, all.y = FALSE)
    current   <- 
      current %>% 
      ungroup %>% 
      select(SITE = Scanner, Age_at_this_scan, Sex, ICV, Group) 
    
    
    #use neuroCombat
    result  <- py$harmonizationLearn(dat, current, smooth_terms = "Age_at_this_scan")
    colnames(result[[2]]) <- colnames(dat)
    current <- cbind(forcombat$Master_subject_ID, forcombat$Time, current, result[[2]]) %>% as.data.frame()
    colnames(current)[1:2] <- c("Master_subject_ID", "Time")
    current[,-1] <- sapply(current[,-1], as.numeric)
    
    #test for group difference and save p-value
    if (is_corthick == "no" & intercept_only == "yes"){
      res <- test_dif(current, res)
    }else if (is_corthick == "yes" & intercept_only == "yes"){
      res <- test_dif_corthick(current, res)  
    }else if (is_corthick == "no" & intercept_only == "no"){
      res <- test_dif_slope(current, res) 
    }else {
      res <- test_dif_slope_corthick(current, res) 
    }
    
  }
  
  fp_xharm <- apply(res, 2 ,function(x) sum(x < 0.05)) %>% as.data.frame()
  colnames(fp_xharm) <- "FPR"
  fp_xharm$Method <- "gamcombat"
  fp_xharm$ROI <- row.names(fp_xharm)
  
  return(fp_xharm)
}





#############################################
#                                           #
# Simulate FPR after longitudinal Combat    #
#                                           #
#############################################

sim_longcombat <- function(mymax, forcombat, is_corthick = "no", intercept_only = "yes"){
  
  
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
  
  df        <- forcombat %>% select(Master_subject_ID, Time, Scanner, Age_at_this_scan, Sex, ICV, rois)
  batch     <<- as.factor(df$Scanner)
  sub       <- distinct(df, Master_subject_ID)
  sub$Group <- NA
  n         <- nrow(sub)
  
  #Preparae dataframe to collect results (aka simulated p-values) in
  runs <- seq(1, mymax)
  res <- as.data.frame(matrix(,ncol=7,nrow=length(runs)))
  colnames(res) <- rois
  
  for (i in runs){
    
    #generate random group assignment
    i <<- i
    set.seed(i)
    sub$Group <- sample(c(0,1), replace=TRUE, size=n)
    current <- merge(df, sub, by = "Master_subject_ID", all.x = TRUE, all.y = FALSE)
    batch <- current$Scanner
    
    #use longitudinal combat harmonization
    if (intercept_only == "yes"){
      longres <- longCombat(
        idvar = "Master_subject_ID",
        timevar = "Time",
        batchvar = "Scanner",
        features = rois,
        formula = "Age_at_this_scan + Sex + ICV + Group*Time",
        ranef = "(1 |Master_subject_ID)",
        data = current,
        niter = 30,
        method = 'REML',
        verbose = FALSE
      )
    
    }else { #add random slope for time
      longres <- longCombat(
        idvar = "Master_subject_ID",
        timevar = "Time",
        batchvar = "Scanner",
        features = rois,
        formula = "Age_at_this_scan + Sex + ICV + Group*Time",
        ranef = "(1 + Time |Master_subject_ID)",
        data = current,
        niter = 30,
        method = 'REML',
        verbose = FALSE
      )
    }
    current <- cbind(current$Group,current$Age_at_this_scan, current$Sex, longres[[1]]) %>% as.data.frame()
    colnames(current)[1:3] <- c("Group", "Age_at_this_scan", "Sex")
    colnames(current) <- gsub(".combat", "", colnames(current))
    
    #compare groups and save p-value
    if (is_corthick == "no" & intercept_only == "yes"){
      res <- test_dif(current, res)
    }else if (is_corthick == "yes" & intercept_only == "yes"){
      res <- test_dif_corthick(current, res)  
    }else if (is_corthick == "no" & intercept_only == "no"){
      res <- test_dif_slope(current, res) 
    }else {
      res <- test_dif_slope_corthick(current, res) 
    }
  }
  
  fp_longharm <- apply(res, 2, function(x) sum(x < 0.05)) %>% as.data.frame()
  colnames(fp_longharm) <- "FPR"
  fp_longharm$Method <- "longcombat"
  fp_longharm$ROI <- row.names(fp_longharm)
  
  return(fp_longharm)
}


#############################################
#                                           #
# Compare simulations across harm. methods  #
#                                           #
#############################################

compare_fpr <- function(data, metric){
  
  #get summary statistics
  temp  <- data %>%
    group_by(Method) %>%
    get_summary_stats(FPR, type = "full") %>% 
    select(Method, median, min, max)
  
  #convert into percent
  n <-  100/mymax
  fun <- function(x){
    x*n
  }
  temp[,-1] <- lapply(temp[,-1], fun)
  
  #format digits of summary statistics
  fun <- function(x){
    format(round(x, 1), nsmall = 1)
  }
  temp[,-1] <- lapply(temp[,-1], fun)
  
  temp$median_range <- paste0(temp$median, " (", temp$min, "-", temp$max, ")")
  
  #reshape dataframe in preparation for ANOVA
  temp <- temp %>% select(Method, median_range) %>% t() 
  colnames(temp) <- temp[1,]
  temp <- as.data.frame(temp)
  temp$Metric <- metric
  temp <- temp[2,]
  
  res.fried <- data %>% friedman_test(FPR ~ Method |ROI)
  temp$P <- res.fried$p
    
  return(temp)
}

