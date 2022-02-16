calculate_delta <- function(longdata){
  tlong <- longdata
  tlong <- tlong %>% 
    mutate(across(c("Ventricles", 
                    "SupraCortex", 
                    "SupraWM", 
                    "SupraDGM", 
                    "CerebellumGM", 
                    "CerebellumWM",
                    "Brainstem"), 
                  as.numeric))
  
  #taking the absolute values of the difference (omitting minus and plus sign)
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(dif_Ventricles = Ventricles - lag(Ventricles))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(dif_SupraCortex = SupraCortex - lag(SupraCortex))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>%
    mutate(dif_SupraWM = SupraWM - lag(SupraWM))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>%
    mutate(dif_SupraDGM = SupraDGM - lag(SupraDGM))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(dif_CerebellumGM = CerebellumGM - lag(CerebellumGM))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(dif_CerebellumWM = CerebellumWM - lag(CerebellumWM))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(dif_Brainstem = Brainstem - lag(Brainstem))
  

  
  #taking the absolute values of the difference (omitting minus and plus sign)
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(delta_Ventricles = abs(Ventricles - lag(Ventricles)))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(delta_SupraCortex = abs(SupraCortex - lag(SupraCortex)))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>%
    mutate(delta_SupraWM = abs(SupraWM - lag(SupraWM)))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>%
    mutate(delta_SupraDGM = abs(SupraDGM - lag(SupraDGM)))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(delta_CerebellumGM = abs(CerebellumGM - lag(CerebellumGM)))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(delta_CerebellumWM = abs(CerebellumWM - lag(CerebellumWM)))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(delta_Brainstem = abs(Brainstem - lag(Brainstem)))
  
  return(tlong)
}


make_shortdata <- function(longdata){
  
  tshort <- longdata %>% filter(Scan_num == 2)
  
  tshort$perc_Ventricles   <- tshort$delta_Ventricles/tshort$Ventricles*100
  tshort$perc_SupraCortex  <- tshort$delta_SupraCortex/tshort$SupraCortex*100
  tshort$perc_SupraWM      <- tshort$delta_SupraWM/tshort$SupraWM*100
  tshort$perc_SupraDGM     <- tshort$delta_SupraDGM/tshort$SupraDGM*100
  tshort$perc_CerebellumGM <- tshort$delta_CerebellumGM/tshort$CerebellumGM*100
  tshort$perc_CerebellumWM <- tshort$delta_CerebellumWM/tshort$CerebellumWM*100
  tshort$perc_Brainstem    <- tshort$delta_Brainstem/tshort$Brainstem*100
  
  return(tshort)
}



calculate_CoV <- function(longdata){
  tlong <- longdata
  tlong <- tlong %>% 
    mutate(across(c("Ventricles", 
                    "SupraCortex", 
                    "SupraWM", 
                    "SupraDGM", 
                    "CerebellumGM", 
                    "CerebellumWM",
                    "Brainstem"), 
                  as.numeric))
  
  #For each subject calculate the CoV in %
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(CoV_Ventricles = sd(Ventricles)/mean(Ventricles)*100)
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(CoV_SupraCortex = sd(SupraCortex)/mean(SupraCortex)*100)
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>%
    mutate(CoV_SupraWM = sd(SupraWM)/mean(SupraWM)*100)
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>%
    mutate(CoV_SupraDGM = sd(SupraDGM)/mean(SupraDGM)*100)
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(CoV_CerebellumGM = sd(CerebellumGM)/mean(CerebellumGM)*100)
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(CoV_CerebellumWM = sd(CerebellumWM)/mean(CerebellumWM)*100)
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(CoV_Brainstem = sd(Brainstem)/mean(Brainstem)*100)
  
  #select only one line per subject
  tshort <- tlong %>% filter(Scan_num == 2)
  
  return(tshort)
}


make_shortdata <- function(longdata){
  
  tshort <- longdata %>% filter(Scan_num == 2)
  
  tshort$perc_Ventricles   <- tshort$delta_Ventricles/tshort$Ventricles*100
  tshort$perc_SupraCortex  <- tshort$delta_SupraCortex/tshort$SupraCortex*100
  tshort$perc_SupraWM      <- tshort$delta_SupraWM/tshort$SupraWM*100
  tshort$perc_SupraDGM     <- tshort$delta_SupraDGM/tshort$SupraDGM*100
  tshort$perc_CerebellumGM <- tshort$delta_CerebellumGM/tshort$CerebellumGM*100
  tshort$perc_CerebellumWM <- tshort$delta_CerebellumWM/tshort$CerebellumWM*100
  tshort$perc_Brainstem    <- tshort$delta_Brainstem/tshort$Brainstem*100
  
  return(tshort)
}

make_summary_table <- function(shortdata, refdata){
  
  df <- as.data.frame(matrix(,ncol=1,nrow=7))
  colnames(df) <- c("ROI")
  
  df$ROI <- c("Ventricles", 
              "Supratent. WM", 
              "Supratent. Cortex", 
              "Supratent. deep GM", 
              "Cerebellar GM", 
              "Cerebellar WM",
              "Brainstem")
  
  df$Scan_pairs <- nrow(shortdata)
  
  df$Mean <- c(mean(shortdata$CoV_Ventricles),
                 mean(shortdata$CoV_SupraWM),
                 mean(shortdata$CoV_SupraCortex),
                 mean(shortdata$CoV_SupraDGM),
                 mean(shortdata$CoV_CerebellumGM),
                 mean(shortdata$CoV_CerebellumWM),
                 mean(shortdata$CoV_Brainstem))
  
  df$Mean <- format(round(df$Mean,1), nsmall = 1)
  
  df$SD <- c(sd(shortdata$CoV_Ventricles),
               sd(shortdata$CoV_SupraWM),
               sd(shortdata$CoV_SupraCortex),
               sd(shortdata$CoV_SupraDGM),
               sd(shortdata$CoV_CerebellumGM),
               sd(shortdata$CoV_CerebellumWM),
               sd(shortdata$CoV_Brainstem))
  
  df$SD <- format(round(df$SD,1), nsmall = 1)
  
  plusminus <-"\u00b1" 
  Encoding(plusminus)<-"UTF-8"
  
  df$CoV <- paste0(df$Mean, " (", plusminus, df$SD, ")")
  df <- df %>% select(-c(Mean, SD))
  
  df$Delta <- c((mean(shortdata$CoV_Ventricles)   - mean(refdata$CoV_Ventricles)),
                (mean(shortdata$CoV_SupraWM)      - mean(refdata$CoV_SupraWM)),
                (mean(shortdata$CoV_SupraCortex)  - mean(refdata$CoV_SupraCortex)),
                (mean(shortdata$CoV_SupraDGM)     - mean(refdata$CoV_SupraDGM)),
                (mean(shortdata$CoV_CerebellumGM) - mean(refdata$CoV_CerebellumGM)),
                (mean(shortdata$CoV_CerebellumWM) - mean(refdata$CoV_CerebellumWM)),
                (mean(shortdata$CoV_Brainstem)    - mean(refdata$CoV_Brainstem)))
  
  df$Delta <- sprintf("%+.1f", as.numeric(df$Delta))
  
  df$P <- c((t.test(shortdata$CoV_Ventricles,   refdata$CoV_Ventricles)[[3]]),
            (t.test(shortdata$CoV_SupraWM,      refdata$CoV_SupraWM)[[3]]),
            (t.test(shortdata$CoV_SupraCortex,  refdata$CoV_SupraCortex)[[3]]),
            (t.test(shortdata$CoV_SupraDGM,     refdata$CoV_SupraDGM)[[3]]),
            (t.test(shortdata$CoV_CerebellumGM, refdata$CoV_CerebellumGM)[[3]]),
            (t.test(shortdata$CoV_CerebellumWM, refdata$CoV_CerebellumWM)[[3]]),
            (t.test(shortdata$CoV_Brainstem,    refdata$CoV_Brainstem)[[3]]))
  
  df$Padj <- c(p.adjust(df$P, method = "holm")[1],
               p.adjust(df$P, method = "holm")[2],
               p.adjust(df$P, method = "holm")[3],
               p.adjust(df$P, method = "holm")[4],
               p.adjust(df$P, method = "holm")[5],
               p.adjust(df$P, method = "holm")[6],
               p.adjust(df$P, method = "holm")[7])
 
   #https://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/nonpz
  #I will only calculate effect sizes for significant differences
  
  df$Cohensd <- ifelse(df$Padj >= 0.05, "NA", format(round(abs(qnorm(df$Padj)/sqrt(df$Scan_pairs)),2), nsmall = 2))
  
  df$P <- as.numeric(format(round(df$P, 3), nsmall = 3))
  
  df$Padj <- ifelse(df$Padj <= 0.01, 
                    as.numeric(format(round(df$Padj, 3), nsmall = 3)),
                    as.numeric(format(round(df$Padj, 2), nsmall = 2)))
  df$Padj <- ifelse(df$Padj < 0.001, "<0.001", df$Padj)
  df$Padj <- ifelse(df$Padj > 0.99, ">0.99", df$Padj)
  
  df$Size <- ifelse(df$Cohensd == "NA", "ns",df$Cohensd)
  
  return(df)
}

####################################################
#                                                  #
# Make summary table for cohort characteristics    #
#                                                  #
####################################################

make_cohort_table_1 <- function(struc_wide, struc_long, scans, clinical){
  
  df <- df <- as.data.frame(matrix(,ncol=1,nrow=1))
  colnames(df) <- "Subjects"
  
  df$Subjects <- n_distinct(struc_long$Master_subject_ID)
  df$Scans <- nrow(struc_long)
  df$Scanners <- n_distinct(struc_long$Scanner)
  
  df$Time_median <- format(round(median(struc_wide$Days), 1), nsmall = 1) 
  df$Time_min <- format(round(min(struc_wide$Days), 1), nsmall = 1)
  df$Time_max <- format(round(max(struc_wide$Days), 1), nsmall = 1)
  df$Time <- paste0(df$Time_median, " (", df$Time_min, " - ", df$Time_max, ") days")
  
  
  df <- df %>% select(-c(Time_median, Time_min, Time_max))
  
  age.df <- scans %>% select(Scan_ID, Age_at_this_scan)
  struc_long <- merge(struc_long, age.df, by = "Scan_ID", all.x = TRUE, all.y = FALSE)
  
  df$Age <- paste0(median(struc_long$Age_at_this_scan) %>% round(0), 
                   " (",
                   min(struc_long$Age_at_this_scan),
                   " - ",
                   max(struc_long$Age_at_this_scan),
                   ")")
  
  sex.df <- clinical %>% select(Master_subject_ID, Sex)
  struc_long <- merge(struc_long, sex.df, by = "Master_subject_ID", all.x = TRUE, all.y = FALSE)
  
  df$Sex <- paste0(sum(struc_long$Sex=="male")/2, 
                   " (", 
                   ((sum(struc_long$Sex=="male")/2)/df$Subjects*100) %>% round(0), 
                   "%)")
  
  df <- t(df) %>% as.data.frame()
  
  return(df)
}

####################################################
#                                                  #
# Helper function to insert row in flextable       #
#                                                  #
####################################################

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

####################################################
#                                                  #
# format flextable for word                        #
#                                                  #
####################################################

make_flextable <- function(ft){
  
  ft <- ft %>% select(-c(Scan_pairs_within, Scan_pairs_across))
  
  ft <- insertRow(ft, rep("Structural data - Volume", 16), 1)
  ft <- insertRow(ft, rep("Structural data - Cortical thickness", 16), 9)
  ft <- insertRow(ft, rep("DTI data - Mean diffusivity", 16), 17 )
  ft <- insertRow(ft, rep("DTI data - Fractional anisotropy", 16), 25)
  
  
  
  vec <<- c("Structural data - Volume", "Structural data - Cortical thickness", "DTI data - Mean diffusivity", "DTI data - Fractional anisotropy")
  
  colormatrix <- ifelse(ft[, c(6, 10, 14)] == "ns", "lightgreen",
                        ifelse(ft[, c(6, 10, 14)] >= 0.8, "red", 
                        ifelse(ft[, c(6, 10, 14)] >= 0.5, "orange",
                               ifelse(ft[, c(6, 10, 14)] >= 0.2, "yellow", "lightgreen"))))
  
  ft <- 
    ft %>% 
    regulartable() %>%
    set_header_labels(ROI = "Region of interest",   
                      Perc_change_within = "CoV %\nMean (SD)",
                      Perc_change_unharm = "CoV %\nMean (SD)", 
                      Delta_unharm = "\u394\nCoV", 
                      Padj_unharm = "Padj", 
                      Size_unharm = "Scanner\neffect",
                      Perc_change_cross = "CoV %\nMean (SD)", 
                      Delta_cross = "\u394\nCoV", 
                      Padj_cross = "Padj", 
                      Size_cross = "Scanner\neffect",
                      Perc_change_long = "CoV %\nMean (SD)", 
                      Delta_long = "\u394\nCoV", 
                      Padj_long = "Padj", 
                      Size_long = "Scanner\neffect") %>%
    autofit() %>%
    bg(j = c(6, 10, 14), bg=colormatrix) %>%
    bg(i = c(1, 9, 17, 25), bg = "grey") %>%
    merge_h(i =  ~ ROI %in% vec)%>%
    add_header_row(values = c(" "," ", "Unharmonized", "neuroCombat", "longCombat"),
                   colwidths = c(1,1, 4, 4, 4))%>%
    add_header_row(values = c(" ", "Within-scanner", "Across-scanner"),
                   colwidths = c(1, 1, 12)) %>%
    theme_vanilla() %>%
    vline(j = c(1, 2, 6, 10))%>%
    vline_left() %>%
    vline_right(part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    set_table_properties(layout = "autofit")%>% 
    padding(padding = 0, part = "all", padding.right = 3, padding.left = 3) %>%
    align(align = "right", part = "all")%>%
    align(align = "left", j = c(1), part = "all") %>%
    align(align = "center", part = "header", i = c(1, 2))
  
  return(ft)
}
