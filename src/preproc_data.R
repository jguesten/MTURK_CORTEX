#### PREPROC_DATA ####

### PREPROCESS MTURK DATA BASED ON TRIALDATA AND SUBJECTDATA ###
### APPLIES SEVERAL UNIVARIATE FILTERS BASED ON CALCULATED PERFORMANCE MEASURES ###
### GIVES FILTERED OUTPUT ###
# input: trialdata, subjectdata
# output: list: 
# 1:preprocessed trialdata 
# 2:preprocessed subjectdata

preproc_data = function(trialdata, subjectdata) {
  
  ### IMPORT LIBRARIES ###
  require(univOutl)
  
  ### DATA WRANGLING ###
  
  # change name Anon ID to ID
  colnames(trialdata)[which(colnames(trialdata) == 'Anon ID')] = 'ID'
  colnames(subjectdata)[which(colnames(subjectdata) == 'Anon ID')] = 'ID'
  
  subjectdata$ID = as.numeric(numextract(subjectdata$`ID`))
  trialdata$ID = as.numeric(numextract(trialdata$`ID`))
  
  trialdata = trialdata[-(
    which(
      trialdata$Stimulus == "o_533_1.jpg" |
        trialdata$Stimulus == "o_533_2.jpg" |
        trialdata$Stimulus == "o_059_1.jpg" |
        trialdata$Stimulus == "o_362_1.jpg" |
        trialdata$Stimulus == "o_684_1.jpg"
    )
  ), ]
  
  # change name "Condition" to "Item Set"
  subjectdata=rename(subjectdata, "Item_set" = "Condition")
  
  # code experimental control variables
  subjectdata$Random_Seq = NA
  subjectdata$Random_Seq[which(subjectdata$Counterbalance == 0 |
                                 subjectdata$Counterbalance == 4)] = 1
  subjectdata$Random_Seq[which(subjectdata$Counterbalance == 1 |
                                 subjectdata$Counterbalance == 5)] = 2
  subjectdata$Random_Seq[which(subjectdata$Counterbalance == 2 |
                                 subjectdata$Counterbalance == 6)] = 3
  subjectdata$Random_Seq[which(subjectdata$Counterbalance == 3 |
                                 subjectdata$Counterbalance == 7)] = 4
  subjectdata$Random_Seq = as.factor(subjectdata$Random_Seq)
  
  subjectdata$Version_order = NA
  subjectdata$Version_order[which(
    subjectdata$Counterbalance == 0 |
      subjectdata$Counterbalance == 1 |
      subjectdata$Counterbalance == 2 |
      subjectdata$Counterbalance == 3
  )] = "Original"
  subjectdata$Version_order[which(
    subjectdata$Counterbalance == 4 |
      subjectdata$Counterbalance == 5 |
      subjectdata$Counterbalance == 6 |
      subjectdata$Counterbalance == 7
  )] = "Reverse"
  
  ## Get exclusion measures ##
  
  # short RTs
  trialdata$Response_Type_Test[(which(trialdata$RT < 200 &
                                        trialdata$RT > 0))] = "NA"
  
  # wrong button press during presentation
  trialdata$Prespress = 0
  trialdata$Prespress[which(
    trialdata$`Trial prop.` == "presentation" &
      trialdata$Response > -1 & !trialdata$Response == 39
  )] = 1
  pressperform = table(trialdata$ID[trialdata$`Trial prop.` == "presentation"], trialdata$Prespress[trialdata$`Trial prop.` ==
                                                                                                      "presentation"])
  percentwrong = pressperform[, 2] / (pressperform[, 2] + pressperform[, 1])
  
  #korrigiert: 0 wird auch als Fehler aufgefasst. Nur Presentation Trials werden betrachtet
  
  # performance table
  testdata = subset(trialdata, `Trial prop.` == 'test')
  ID = cbind(unique(testdata$ID))
  sperform = cbind(ID, table(testdata$ID, testdata$Response_Type_Test))
  colnames(sperform)[1] = "ID"
  
  # percent of old button press
  temp <- trialdata[which(trialdata$`Trial prop.` == "test"), ]
  temp <- as.data.frame(table(temp$Response, temp$ID))
  temp <- temp %>% spread(Var1, Freq)
  old_perc = temp$`37` / (temp$`39` + temp$`37`)
  
  # false responses during Reaction Time Task
  RT_data = trialdata[trialdata$`Trial prop.` == "RT_Task", ]
  falsedata = table(RT_data$Response_Type_Presentation, RT_data$ID)
  falsedata = as.data.frame(falsedata)
  falsedata = falsedata[which(falsedata$Var1 == "false"), ]
  falsedata = falsedata %>% spread(Var1, Freq)
  falsetrain = falsedata$false
  
  temp$falsetrain = 0
  temp$falsetrain[match(falsedata$Var2, temp$Var2)] = falsetrain
  falsetrain = temp$falsetrain
  
  # get performance and RT measures for RT task
  RT_perform = table(RT_data$Response_Type_Presentation, RT_data$ID) %>% as.data.frame() %>% spread(Var1, Freq) %>% mutate(pc =
                                                                                                                             correct / (false + correct)) %>% mutate(pc_miss_inc = correct / 12)
  mean(RT_perform$pc_miss_inc, na.rm = T)
  sd(RT_perform$pc_miss_inc, na.rm = T)
  
  RT_RT = RT_data %>% filter(RT > 0) %>% group_by(ID) %>% dplyr::summarize(mean_RT = mean(RT, na.rm = TRUE))
  mean(RT_RT$mean_RT)
  sd(RT_RT$mean_RT)
  
  
  # no responses are counted in addition to false responses
  if (identical(as.vector(temp$Var2), names(percentwrong))) {
    pressperform = cbind(percentwrong, cbind(falsetrain, old_perc))
  }
  if (identical(rownames(pressperform), rownames(sperform))) {
    sperform = cbind(sperform, pressperform)
  }
  
  # percent total responses
  given_responses = c(sperform[, colnames(sperform) == "correct_rejection"] +
                        sperform[, colnames(sperform) == "false_alarm"] + sperform[, colnames(sperform) ==
                                                                                     "hit"] + sperform[, colnames(sperform) == "miss"])
  failed_responses = c(sperform[, colnames(sperform) == "-"])
  total_responses = (given_responses / (given_responses + failed_responses)) *
    100
  sperform = cbind(sperform, total_responses)
  sperform = as.data.frame(sperform)
  
  # corrected hit rate
  sperform$HR = (sperform$hit / (sperform$hit + sperform$miss))
  sperform$FAR = (sperform$false_alarm / (sperform$false_alarm + sperform$correct_rejection))
  sperform$CHR = sperform$HR - sperform$FAR
  
  subjectdata = merge(subjectdata, sperform, no.dupl = T)
  
  # delete all subjects with missing values
  subjectdata = subjectdata[complete.cases(subjectdata), ]
  
  ## exclude data ##
  
  # old perc
  outperc_ratio = LocScaleB(as.numeric(subjectdata$old_perc))
  
  # falsetrain
  outperc_falsetrain = which(subjectdata$falsetrain > 1)
  
  # percentwrong
  outperc_percentwrong = LocScaleB(as.numeric(subjectdata$percentwrong))
  
  # performance (corrected hit rate)
  outperc_CHR = LocScaleB(as.numeric(subjectdata$CHR))
  outperc_CHR = which(subjectdata$CHR < outperc_CHR$bounds[1])
  outperc_CHR
  
  # total responses
  outperc_ttlresp = LocScaleB(as.numeric(subjectdata$total_responses))
  
  # remove outliers
  subjectdata = subjectdata[-(union(
    union(union(
      union(outperc_falsetrain, outperc_ttlresp$outliers),
      outperc_CHR
    ), outperc_ratio$outliers),
    outperc_percentwrong$outliers
  )), ]
  nrow(subjectdata)
  
  # exclude due to same button press
  same_exclude = which(subjectdata$correct_rejection == 0 &
                         subjectdata$miss == 0)
  same_exclude = union(same_exclude,
                       which(subjectdata$hit == 0 & subjectdata$false_alarm == 0))
  if (!identical(same_exclude, integer(0)))
  {
    subjectdata = subjectdata[-same_exclude, ]
  }
  
  # exclude incomplete data
  subjectdata = subjectdata[complete.cases(subjectdata), ]
  
  # remove performance columns
  subjectdata = dplyr::select(subjectdata,-c(Counterbalance,HR,FAR,CHR))
  
  ## apply subject selection to trialdata
  
  goodtrials = c()
  for (i in 1:length(subjectdata$ID)) {
    goodtrials = c(goodtrials, which(trialdata$ID == subjectdata$ID[i]))
  }
  trialdata = trialdata[goodtrials, ]
  
  return(list(trialdata, subjectdata))
}
