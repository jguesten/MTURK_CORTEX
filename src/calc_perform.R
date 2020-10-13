#### CALC_PERFORM ####

### FUNCTION TO CALCULATE PERFORMANCE MEASURES ###
# input: trialdata, subjectdata
# output: subjectdata dataframe with added domain_specific performance measures

Calc_perform = function(trialdata, subjectdata) {
  subject_data = split(trialdata, trialdata$ID)
  subject_results = data.frame()
  temp_results = data.frame()
  
  
  cnt = 0
  for (name in names(subject_data)) {
    
    HR_o = double()
    H_cnt_o = double()
    M_cnt_o = double()
    FAR_o = double()
    FA_cnt_o = double()
    CR_cnt_o = double()
    HR_s = double()
    H_cnt_s = double()
    M_cnt_s = double()
    FAR_s = double()
    FA_cnt_s = double()
    CR_cnt_s = double()
    RT_H_o = double()
    RT_M_o = double()
    RT_FA_o = double()
    RT_CR_o = double()
    RT_H_s = double()
    RT_M_s = double()
    RT_FA_s = double()
    RT_CR_s = double()
    
    indi = which(
      subject_data[[name]]$Response_Type_Test == "hit" |
        subject_data[[name]]$Response_Type_Test == "miss" |
        subject_data[[name]]$Response_Type_Test == "false_alarm" |
        subject_data[[name]]$Response_Type_Test == "correct_rejection"
    )
    if (!length(indi) == 0L) {
      cnt = cnt + 1
      #get deviation from mean stimulus FA HR and Perc Corr, for Objects and Scenes separately
      for (i in 1:length(indi)) {
        RT=subject_data[[name]]$RT[indi[i]]
        
        if (subject_data[[name]]$`Stim_Type(object=1, scene=2)`[indi[i]] ==
            1) {
          if (subject_data[[name]]$Response_Type_Test[indi[i]] == "hit" |
              subject_data[[name]]$Response_Type_Test[indi[i]] == "miss") {
            if (subject_data[[name]]$Response_Type_Test[indi[i]] == "hit") {
              H_cnt_o[indi[i]] = 1
              RT_H_o[indi[i]] = RT
            } else if (subject_data[[name]]$Response_Type_Test[indi[i]] ==
                       "miss") {
              M_cnt_o[indi[i]] = 1
              RT_M_o[indi[i]] = RT
            }
          } else if (subject_data[[name]]$Response_Type_Test[indi[i]] ==
                     "false_alarm" |
                     subject_data[[name]]$Response_Type_Test[indi[i]] == "correct_rejection") {
            if (subject_data[[name]]$Response_Type_Test[indi[i]] == "false_alarm") {
              FA_cnt_o[indi[i]] = 1
              RT_FA_o[indi[i]] = RT
            } else if (subject_data[[name]]$Response_Type_Test[indi[i]] ==
                       "correct_rejection") {
              CR_cnt_o[indi[i]] = 1
              RT_CR_o[indi[i]] = RT
            }
          }
        } else if (subject_data[[name]]$`Stim_Type(object=1, scene=2)`[indi[i]] ==
                   2) {
          if (subject_data[[name]]$Response_Type_Test[indi[i]] == "hit" |
              subject_data[[name]]$Response_Type_Test[indi[i]] == "miss") {
            if (subject_data[[name]]$Response_Type_Test[indi[i]] == "hit") {
              H_cnt_s[indi[i]] = 1
              RT_H_s[indi[i]] = RT
              
            } else if (subject_data[[name]]$Response_Type_Test[indi[i]] ==
                       "miss") {
              M_cnt_s[indi[i]] = 1
              RT_M_s[indi[i]] = RT
            }
          } else if (subject_data[[name]]$Response_Type_Test[indi[i]] ==
                     "false_alarm" |
                     subject_data[[name]]$Response_Type_Test[indi[i]] == "correct_rejection") {
            if (subject_data[[name]]$Response_Type_Test[indi[i]] == "false_alarm") {
              FA_cnt_s[indi[i]] = 1
              RT_FA_s[indi[i]] = RT
            } else if (subject_data[[name]]$Response_Type_Test[indi[i]] ==
                       "correct_rejection") {
              CR_cnt_s[indi[i]] = 1
              RT_CR_s[indi[i]] = RT
            }
          }
        }
      }
      
      H_o = sum(H_cnt_o, na.rm = TRUE)
      M_o = sum(M_cnt_o, na.rm = TRUE)
      HR_o = H_o / (H_o + sum(M_cnt_o, na.rm = TRUE))
      FA_o = sum(FA_cnt_o, na.rm = TRUE)
      FAR_o = FA_o / (FA_o + sum(CR_cnt_o, na.rm =TRUE))
      CR_o = sum(CR_cnt_o, na.rm = TRUE)
      H_s = sum(H_cnt_s, na.rm = TRUE)
      M_s = sum(M_cnt_s, na.rm = TRUE)
      HR_s = H_s / (H_s + sum(M_cnt_s, na.rm = TRUE))
      FA_s = sum(FA_cnt_s, na.rm = TRUE)
      FAR_s = FA_s / (FA_s + sum(CR_cnt_s, na.rm =
                                   TRUE))
      CR_s = sum(CR_cnt_s, na.rm = TRUE)
      PR_o = HR_o - FAR_o
      PR_s = HR_s - FAR_s
      RT_H_o = mean(RT_H_o[!is.na(RT_H_o)])
      RT_M_o = mean(RT_M_o[!is.na(RT_M_o)])
      RT_FA_o = mean(RT_FA_o[!is.na(RT_FA_o)])
      RT_CR_o = mean(RT_CR_o[!is.na(RT_CR_o)])
      RT_H_s = mean(RT_H_s[!is.na(RT_H_s)])
      RT_M_s = mean(RT_M_s[!is.na(RT_M_s)])
      RT_FA_s = mean(RT_FA_s[!is.na(RT_FA_s)])
      RT_CR_s = mean(RT_CR_s[!is.na(RT_CR_s)])
      
      # concatenate object and scene vectors
      
      PR = c(PR_o, PR_s)
      RT_H = c(RT_H_o, RT_H_s)
      RT_M = c(RT_M_o, RT_M_s)
      RT_FA = c(RT_FA_o, RT_FA_s)
      RT_CR = c(RT_CR_o, RT_CR_s)
      HR = c(HR_o, HR_s)
      FAR = c(FAR_o, FAR_s)
      
      Domain = c("Object", "Scene")
      subdata_row = as.data.frame(subset(subjectdata,ID==name, drop=TRUE))
      temp_results = tibble(PR,RT_H,RT_M,RT_FA,RT_CR,HR,FAR,Domain)
      temp_results = cbind(rbind(subdata_row,subdata_row),temp_results)
      
      
      if(plyr::empty(subject_results)){
        subject_results=temp_results} else {
          
          subject_results=rbind(subject_results,temp_results) 
          
        }
      
    }
  }
  
  return(subject_results)
}
