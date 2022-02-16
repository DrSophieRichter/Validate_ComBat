calculate_delta_corthick <- function(longdata){
  tlong <- longdata
  tlong <- tlong %>% 
    mutate(across(c("Frontal", 
                    "Parietal", 
                    "Insular", 
                    "Occipital", 
                    "Temporal", 
                    "Hippocampal",
                    "WholeCortex"), 
                  as.numeric))
  
  #taking the absolute values of the difference (omitting minus and plus sign)
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(dif_Frontal = Frontal - lag(Frontal))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(dif_Parietal = Parietal - lag(Parietal))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>%
    mutate(dif_Insular = Insular - lag(Insular))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>%
    mutate(dif_Occipital = Occipital - lag(Occipital))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(dif_Temporal = Temporal - lag(Temporal))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(dif_Hippocampal = Hippocampal - lag(Hippocampal))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(dif_WholeCortex = WholeCortex - lag(WholeCortex))
  

  
  #taking the absolute values of the difference (omitting minus and plus sign)
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(delta_Frontal = abs(Frontal - lag(Frontal)))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(delta_Parietal = abs(Parietal - lag(Parietal)))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>%
    mutate(delta_Insular = abs(Insular - lag(Insular)))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>%
    mutate(delta_Occipital = abs(Occipital - lag(Occipital)))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(delta_Temporal = abs(Temporal - lag(Temporal)))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(delta_Hippocampal = abs(Hippocampal - lag(Hippocampal)))
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(delta_WholeCortex = abs(WholeCortex - lag(WholeCortex)))
  
  return(tlong)
}

calculate_CoV_corthick <- function(longdata){
  tlong <- longdata
  tlong <- tlong %>% 
    mutate(across(c("Frontal", 
                    "Parietal", 
                    "Insular", 
                    "Occipital", 
                    "Temporal", 
                    "Hippocampal",
                    "WholeCortex"), 
                  as.numeric))
  
  #taking the absolute values of the difference (omitting minus and plus sign)
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(CoV_Frontal = sd(Frontal)/mean(Frontal)*100)
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(CoV_Parietal = sd(Parietal)/mean(Parietal)*100)
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>%
    mutate(CoV_Insular = sd(Insular)/mean(Insular)*100)
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>%
    mutate(CoV_Occipital = sd(Occipital)/mean(Occipital)*100)
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(CoV_Temporal = sd(Temporal)/mean(Temporal)*100)
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(CoV_Hippocampal = sd(Hippocampal)/mean(Hippocampal)*100)
  
  tlong <- tlong %>% group_by(Master_subject_ID) %>% 
    mutate(CoV_WholeCortex = sd(WholeCortex)/mean(WholeCortex)*100)
  
  tshort <- tlong %>% filter(Scan_num == 2)
  
  return(tshort)
}

make_shortdata_corthick <- function(longdata){
  
  tshort <- longdata %>% filter(Scan_num == 2)
  
  tshort$perc_Frontal   <- tshort$delta_Frontal/tshort$Frontal*100
  tshort$perc_Parietal  <- tshort$delta_Parietal/tshort$Parietal*100
  tshort$perc_Insular      <- tshort$delta_Insular/tshort$Insular*100
  tshort$perc_Occipital     <- tshort$delta_Occipital/tshort$Occipital*100
  tshort$perc_Temporal <- tshort$delta_Temporal/tshort$Temporal*100
  tshort$perc_Hippocampal <- tshort$delta_Hippocampal/tshort$Hippocampal*100
  tshort$perc_WholeCortex    <- tshort$delta_WholeCortex/tshort$WholeCortex*100
  
  return(tshort)
}



make_summary_table_corthick <- function(shortdata, refdata){
  
  df <- as.data.frame(matrix(,ncol=1,nrow=7))
  colnames(df) <- c("ROI")
  
  df$ROI <- c("Frontal", 
              "Insular", 
              "Parietal", 
              "Occipital", 
              "Temporal", 
              "Hippocampal",
              "WholeCortex")
  
  df$Scan_pairs <- nrow(shortdata)
  
  df$Mean <- c(mean(shortdata$CoV_Frontal),
                 mean(shortdata$CoV_Insular),
                 mean(shortdata$CoV_Parietal),
                 mean(shortdata$CoV_Occipital),
                 mean(shortdata$CoV_Temporal),
                 mean(shortdata$CoV_Hippocampal),
                 mean(shortdata$CoV_WholeCortex))
  
  df$Mean <- format(round(df$Mean,1), nsmall = 1)
  
  df$SD <- c(sd(shortdata$CoV_Frontal),
               sd(shortdata$CoV_Insular),
               sd(shortdata$CoV_Parietal),
               sd(shortdata$CoV_Occipital),
               sd(shortdata$CoV_Temporal),
               sd(shortdata$CoV_Hippocampal),
               sd(shortdata$CoV_WholeCortex))
  
  df$SD <- format(round(df$SD,1), nsmall = 1)
  
  plusminus <-"\u00b1" 
  Encoding(plusminus)<-"UTF-8"
  
  df$CoV <- paste0(df$Mean, " (", plusminus, df$SD, ")")
  df <- df %>% select(-c(Mean, SD))
  
  df$Delta <- c((mean(shortdata$CoV_Frontal)   - mean(refdata$CoV_Frontal)),
                (mean(shortdata$CoV_Insular)      - mean(refdata$CoV_Insular)),
                (mean(shortdata$CoV_Parietal)  - mean(refdata$CoV_Parietal)),
                (mean(shortdata$CoV_Occipital)     - mean(refdata$CoV_Occipital)),
                (mean(shortdata$CoV_Temporal) - mean(refdata$CoV_Temporal)),
                (mean(shortdata$CoV_Hippocampal) - mean(refdata$CoV_Hippocampal)),
                (mean(shortdata$CoV_WholeCortex)    - mean(refdata$CoV_WholeCortex)))
  
  df$Delta <- sprintf("%+.1f", as.numeric(df$Delta))
  
  df$P <- c((t.test(shortdata$CoV_Frontal,   refdata$CoV_Frontal)[[3]]),
            (t.test(shortdata$CoV_Insular,      refdata$CoV_Insular)[[3]]),
            (t.test(shortdata$CoV_Parietal,  refdata$CoV_Parietal)[[3]]),
            (t.test(shortdata$CoV_Occipital,     refdata$CoV_Occipital)[[3]]),
            (t.test(shortdata$CoV_Temporal, refdata$CoV_Temporal)[[3]]),
            (t.test(shortdata$CoV_Hippocampal, refdata$CoV_Hippocampal)[[3]]),
            (t.test(shortdata$CoV_WholeCortex,    refdata$CoV_WholeCortex)[[3]]))
  
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
  
  df$Size <- ifelse(df$Cohensd == "NA", "ns", df$Cohensd)
  
  return(df)
}