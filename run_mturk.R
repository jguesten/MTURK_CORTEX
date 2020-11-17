#### SCRIPT PERFORMS DATA PREPROCESSING, CALCULATES PERFORMANCE MEASURES ####
#### AND RUNS LINEAR MIXED EFFECTS MODELS FOR PERFORMANCE AND RT DATA ####

rm(list = ls(all.names = TRUE))

#### LOAD DEPENDENCIES ####

# load environment paths
source("src/get_env.R")

# libraries
require(tidyverse)
require(lme4)
require(emmeans)

# custom functions
source("src/calc_perform.R")
source("src/numextract.R")
source("src/preproc_data.R")
source("src/sstest.R")


#### IMPORT DATA ####

trialdata <-
  read_delim(
    paste0(getwd(), "/Logfiles/Final_Logs/trialdata.csv"),
    ";",
    escape_double = FALSE,
    trim_ws = TRUE
  )

subjectdata <-
  read_delim(
    paste0(getwd(), "/Logfiles/Final_Logs/subjectdata.csv"),
    ";",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>% as.data.frame()


### PREPROCESS DATA ###

preproc_list = preproc_data(trialdata, subjectdata)
preproc_trialdata = preproc_list[[1]]
preproc_subjectdata = preproc_list[[2]]

### CALCULATE PERFORMANCE MEASURES ###

# Calculate (domain-wise) performance measures
subject_results = Calc_perform(preproc_trialdata, preproc_subjectdata)

# Code categorical predictors as factor
subject_results$Item_set = as.factor(subject_results$Item_set)
subject_results$Domain = as.factor(subject_results$Domain)
subject_results$Browser = as.factor(subject_results$Browser)
subject_results$Handedness = as.factor(subject_results$Handedness)
subject_results$Platform = as.factor(subject_results$Platform)
subject_results$Gender = as.factor(subject_results$Gender)
subject_results$Random_Seq = as.factor(subject_results$Random_Seq)
subject_results$Version_order = as.factor(subject_results$Version_order)


### RUN LINEAR MIXED EFFECTS MODEL ###

## PERFORMANCE MEASURES ##

run_LME = function(data, addstring = NULL) {
  DV = c("HR", "FAR", "PR")
  
  # descriptives
  desc = lapply(DV, function(x) {
    print(cbind(
      rbind(psych::describe(data[[x]][data$Domain == "Object"]), psych::describe(data[[x]][data$Domain ==
                                                                                             "Scene"])),
      domain = c("object", "scene"),
      var = c(x, x)
    ), digits = 4)
  })
  
  # LME
  lmer_mods = lapply(DV, function(x) {
    f = as.formula(
      paste(
        x,
        "~Browser+Platform+Gender+Handedness+Item_set*Version_order*Random_Seq+Age*Domain+(1|ID)"
      )
    )
    model = eval(bquote(lmer(.(f), data = data)))
  })
  
  # simple slopes
  slopes = lapply(lmer_mods, function(x) {
    sstest(x)
  })
  
  # summary
  modsum = lapply(lmer_mods, function(x) {
    summary(x)
  })
  
  # anova
  anova = lapply(lmer_mods, function(x) {
    anova(x, type = "3")
  })
  
  # save output
  # capture.output(
  #   list(desc, anova, modsum, slopes),
  #   file = paste0(result_path, "total_output", addstring, ".doc")
  # )
  lmer_mods
}

# run LME on entire sample
lmer_mods = run_LME(subject_results)

# run LME without msie subs
noie_subs = subject_results[-which(subject_results$Browser == 'msie'),]
run_LME(noie_subs, '_no_msie')

# no floor effect
no_floor_subs = subject_results[-which(subject_results$PR < 0),]
run_LME(no_floor_subs, '_no_floor')


## REACTION TIMES ##

# convert to long format
subject_results_long = subject_results %>% tidyr::pivot_longer(cols =
                                                                 RT_H:RT_CR)
subject_results_long$condition = factor(subject_results_long$condition)

# descriptives
desc_domain = print(cbind(
  rbind(
    psych::describe(subject_results_long[["value"]][subject_results_long$domain ==
                                                      "Object"]),
    psych::describe(subject_results_long[["value"]][subject_results_long$domain ==
                                                      "Scene"])
  ),
  domain = c("object", "scene"),
  var = c("value", "value")
), digits = 4)
desc_cond = as.data.frame(describeBy(subject_results))[match(c("RT_H", "RT_M", "RT_FA", "RT_CR"), rownames(describeBy(subject_results))),]
desc_gender = print(cbind(
  rbind(
    psych::describe(subject_results_long[["value"]][subject_results_long$Gender ==
                                                      "male"]),
    psych::describe(subject_results_long[["value"]][subject_results_long$Gender ==
                                                      "female"])
  ),
  domain = c("male", "female"),
  var = c("value", "value")
), digits = 4)
desc_platform = print(cbind(
  rbind(
    psych::describe(subject_results_long[["value"]][subject_results_long$Platform ==
                                                      "android"]),
    psych::describe(subject_results_long[["value"]][subject_results_long$Platform ==
                                                      "windows"]),
    psych::describe(subject_results_long[["value"]][subject_results_long$Platform ==
                                                      "chromeos"]),
    psych::describe(subject_results_long[["value"]][subject_results_long$Platform ==
                                                      "macos"]),
    psych::describe(subject_results_long[["value"]][subject_results_long$Platform ==
                                                      "linux"])
  ),
  domain = c("android", "windows", "chromeos", "mac", "linux"),
  var = c("value", "value", "value", "value", "value")
), digits = 4)

# LME
lmer_mod_RT = lmer(
  value ~ Age * Domain + Age * condition + Browser + Platform + Gender + Handedness +
    Condition * version_order * Random_Seq + (1 |
                                                ID),
  data = subject_results_long
)

# get slopes for simple contrasts
cond_Age_trends = emtrends(lmer_mod_RT,
                           "condition",
                           var = "Age",
                           options = list())
domain_Age_trends = emtrends(lmer_mod_RT, "Domain", var = "Age", options = list())
cond_Age_pairwise = pairs(emtrends(
  lmer_mod_RT,
  "condition",
  var = "Age",
  options = list()
))

anova_RT = anova(lmer_mod_RT)
capture.output(
  list(
    desc_domain,
    desc_cond,
    desc_platform,
    desc_sex,
    cond_Age_trends,
    cond_Age_pairwise,
    domain_Age_trends,
    lmer_mod_RT,
    anova_RT
  ),
  file = paste0(result_path, "RT_output.doc")
)
