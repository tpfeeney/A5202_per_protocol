## ---------------------------
##
## Script name: recreation of Daar 2011 results
##
## Purpose of script:
##
## Author: Timothy Feeney
##
## Date Created: 2024-04-17
##
## Copyright (c) Timothy Feeney, 2024
## Email: feeney@unc.edu
##
## ---------------------------
##
## Notes: this is meant to recreate, roughly, the results from Daar 2011 as a proof of concept. 
## this will focus on recreating the results for virologic failure using the long dataset made in
## long_dataset_creation.R
##
## ---------------------------

## set working directory for Mac and PC

# setwd("~/")          # working directory (mac)
# setwd("C:/Users/")      # working directory (PC)

## ---------------------------

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation
#memory.limit(30000000)       # this is needed on some PCs to increase memory allowance, but has no impact on macs.

## ---------------------------

## list the packages we will need: 

packages<-c('ggplot2', 'dplyr', 'survival', 'ggsurvfit', 'here', 'forcats')      # list all the packages we need

## ---------------------------

## load up our functions into memory

for (i in packages){
    require(i, character.only=T)
} 

## :::::::::::::::::::::::::::::::::::::::
## Establish Directories to use----
## :::::::::::::::::::::::::::::::::::::::

results_tabs_dir<-paste0('/Users/thefeeney/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/',
                         'Dissertation/DissertationAnalysis/Aim 2/results/DR098/tables/')
results_figs_dir<-paste0('/Users/thefeeney/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/',
                         'Dissertation/DissertationAnalysis/Aim 2/results/DR098/figures/')
data_dir<-paste0('/Users/thefeeney/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/',
                 'Dissertation/DissertationAnalysis/Aim 2/data/')
analysis_dir<-paste0('/Users/thefeeney/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/',
                     'Dissertation/DissertationAnalysis/Aim 2/analysis')

## ---------------------------
## LOAD DATA----
## ---------------------------
longdf<-read.csv(paste0(data_dir,'last_week_long_data2_DR098-2.csv'))

## :::::::::::::::::::::::::::::::::::::::
## DATA DESCRIPTION----
## :::::::::::::::::::::::::::::::::::::::

# Take a look and make sure there is a treatment value for each person
longdf %>%group_by(patid) %>% slice(1) %>%  select(abc_tdf) %>% Hmisc::describe()

#abc_tdf 
#.       n  missing distinct 
#.    1857        0        2 
#
#Value      ABC/3TC TDF/FTC
#Frequency      928     929
#Proportion     0.5     0.5

# Look at the number being censored
cens_df<-longdf %>% group_by(patid) %>% slice(1) %>% ungroup() %>% mutate(cens_indic=if_else(virfail==0 & !(patid%in%protocol_finished),1,0)) %>% {gmodels::CrossTable(.$cens_indic,.$abc_tdf)}
longdf %>% group_by(patid) %>% slice(1) %>% ungroup() %>% filter(!(patid %in% protocol_finished)) %>% {gmodels::CrossTable(.$ever_cens,.$abc_tdf)}
gmodels::CrossTable(longdf %>% group_by(patid) %>% slice(1) %>% select(ever_cens), longdf %>% group_by(patid) %>% slice(1) %>% select(abc_tdf))

longdf %>% group_by(patid) %>% mutate(cens_indic=if_else(!(patid %in% protocol_finished) &
                                                           ever_dev_doses_1==0 &
                                                           #ever_cens==1 & 
                                                           sum(target_cens_1, na.rm=T)>0 ,1,0),
                                      cens_immed=if_else(target_cens_1==1 & Week==0, 1, 0)
                                      ) %>%  slice(1) %>% ungroup() %>% {gmodels::CrossTable(.$cens_indic,.$abc_tdf)}

longdf %>% group_by(patid) %>% filter(target_cens_1==1 & Week==0) %>% select(patid, Week,reasnot, adhsoewk,adhwks,visitdt,abc_tdf) %>% print(n=100) %>% group_by(abc_tdf) %>% summarise(count=n())

cens_early<-
  longdf %>% group_by(patid) %>% filter(target_cens_1==1 & Week==0) %>%
  select(patid, Week,reasnot, adhsoewk,adhwks,visitdt,abc_tdf) %>%
  print(n=100) %>%
  group_by(abc_tdf) %>%
  pull(patid)

longdf %>% group_by(patid) %>% mutate(cens_indic=if_else(!(patid %in% protocol_finished) &
                                                          #!(patid %in% cens_early)&
                                                           #!is.na(adhsoewk)&
                                                           #ever_dev_doses_1==0 &
                                                           ever_outcome==0 &
                                                           ever_cens==1 & 
                                                           sum(target_cens_1, na.rm=T)>0 ,1,0),
                                      cens_immed=if_else(target_cens_1==1 & Week==0, 1, 0)
) %>%  slice(1) %>% ungroup() %>% {gmodels::CrossTable(.$cens_indic,.$abc_tdf)}


# evaluating the total number censored and total number deviating for a single deviation in the prior week

cens_dev_by_txt<-longdf %>% group_by(patid) %>% 
  mutate(cens_indic=if_else(!(patid %in% protocol_finished) & 
                              ever_dev_doses_1==0 &
                              #ever_cens==1 & 
                              sum(target_cens_1, na.rm=T)>0 ,1,0),
         no_miss=if_else(ever_dev_doses_1==0 ,1,0),
         cens_immed=if_else(target_cens_1==1 & Week==0, 1, 0,),
         one_dose=if_else(ever_dev_doses_1==1 & cens_indic==0,1,0),
         dose_deviation_1=case_when(sum(target_cens_1, na.rm = T)>0 & cens_week_1<outcome_week~1,
                                    TRUE~0),
         two_dose=if_else(ever_dev_doses_2==1 & cens_indic==0,1,0),
         dose_deviation_2=case_when(sum(target_cens_2, na.rm = T)>0 & cens_week_2<outcome_week~1,
                                    TRUE~0),
         three_dose=if_else(ever_dev_doses_3==1 & cens_indic==0,1,0),
         dose_deviation_3=case_when(sum(target_cens_3, na.rm = T)>0 & cens_week_3<outcome_week~1,
                                    TRUE~0),
         five_dose=if_else(ever_dev_doses_5==1 & cens_indic==0,1,0),
         dose_deviation_5=case_when(sum(target_cens_5, na.rm = T)>0 & cens_week_5<outcome_week~1,
                                    TRUE~0),
         ten_dose=if_else(ever_dev_doses_10==1 & cens_indic==0,1,0),
         dose_deviation_10=case_when(sum(target_cens_10, na.rm = T)>0 & cens_week_10<outcome_week~1,
                                    TRUE~0)
         
) %>%  mutate(never_miss=if_else(dose_deviation_1==0,1,0)) %>% slice(1)  %>% ungroup() %>%
    group_by(abc_tdf) %>% select(ever_cens,ever_dev_doses_1:ever_dev_doses_10,patid, abc_tdf,never_miss,dose_deviation_1:dose_deviation_10, no_miss,cens_indic,
                                 one_dose, two_dose, three_dose, five_dose,ten_dose) %>%
    summarise(#'Not missed'    =sum(no_miss),
              "Not Missed"=sum(never_miss),
              'Censored'    =sum(cens_indic),
              #'1 Dose'=sum(one_dose),
              "1 Dose"=sum(dose_deviation_1),
              #'2 Dose'=sum(two_dose),
              "2 Dose"=sum(dose_deviation_2),
              #'3 Dose'=sum(three_dose),
              "3 Dose"=sum(dose_deviation_3),
              #'5 Dose'=sum(five_dose),
              "5 Dose"=sum(dose_deviation_5),
              #'10 Dose'=sum(ten_dose),
              "10 Dose"=sum(dose_deviation_10),
              'Total'=length(unique(patid))) %>%
    ungroup() %>% rename('Treatment Group'="abc_tdf")
print(cens_dev_by_txt)


writeLines(kableExtra::kable(cens_dev_by_txt, #print out tex file
                             'latex',
                             booktabs=T,
                             align = "l", #left justify
                             digits=3),   # 3 sigfigs
           paste0(results_tabs_dir,'cens_dev_by_tx_DR098.tex'))

cens_dev_by_txt_combo<-longdf %>% group_by(patid) %>% 
  mutate(cens_indic=if_else(!(patid %in% protocol_finished) & 
                              ever_dev_doses_1==0 &
                              #ever_cens==1 & 
                              sum(target_cens_1, na.rm=T)>0 ,1,0),
         no_miss=if_else(ever_dev_doses_1==0 ,1,0),
         cens_immed=if_else(target_cens_1==1 & Week==0, 1, 0,),
         one_dose=if_else(ever_dev_doses_1==1 & cens_indic==0,1,0),
         dose_deviation_1=case_when(sum(target_cens_1, na.rm = T)>0 & cens_week_1<outcome_week~1,
                                    TRUE~0),
         two_dose=if_else(ever_dev_doses_2==1 & cens_indic==0,1,0),
         dose_deviation_2=case_when(sum(target_cens_2, na.rm = T)>0 & cens_week_2<outcome_week~1,
                                    TRUE~0),
         three_dose=if_else(ever_dev_doses_3==1 & cens_indic==0,1,0),
         dose_deviation_3=case_when(sum(target_cens_3, na.rm = T)>0 & cens_week_3<outcome_week~1,
                                    TRUE~0),
         five_dose=if_else(ever_dev_doses_5==1 & cens_indic==0,1,0),
         dose_deviation_5=case_when(sum(target_cens_5, na.rm = T)>0 & cens_week_5<outcome_week~1,
                                    TRUE~0),
         ten_dose=if_else(ever_dev_doses_10==1 & cens_indic==0,1,0),
         dose_deviation_10=case_when(sum(target_cens_10, na.rm = T)>0 & cens_week_10<outcome_week~1,
                                     TRUE~0)
         ) %>%  mutate(never_miss=if_else(dose_deviation_1==0,1,0)) %>% 
  slice(1) %>% ungroup() %>%
    group_by(abc_tdf) %>% select(ever_cens,ever_target_cens_1:ever_target_cens_10, patid,abc_tdf,dose_deviation_1:dose_deviation_10, cens_indic, one_dose, two_dose,three_dose,five_dose,ten_dose) %>%
    summarise(#'1 Dose +cens'=sum(cens_indic,one_dose),
              '1 Dose +cens'=sum(cens_indic,dose_deviation_1),
              #'2 Dose +cens'=sum(cens_indic, two_dose),
              '2 Dose +cens'=sum(cens_indic, dose_deviation_2),
              #'3 Dose +cens'=sum(cens_indic, three_dose),
              '3 Dose +cens'=sum(cens_indic, dose_deviation_3),
              #'5 Dose +cens'=sum(cens_indic, five_dose),
              '5 Dose +cens'=sum(cens_indic, dose_deviation_5),
              #'10 Dose +cens'=sum(cens_indic, ten_dose),
              '10 Dose +cens'=sum(cens_indic, dose_deviation_10),
              'Total'=length(unique(patid))) %>%
    ungroup() %>% rename('Treatment Group'="abc_tdf")
print(cens_dev_by_txt_combo )


writeLines(kableExtra::kable(cens_dev_by_txt_combo, #print out tex file
                             'latex',
                             booktabs=T,
                             align = "l", #left justify
                             digits=3),   # 3 sigfigs
           paste0(results_tabs_dir,'cens_dev_by_tx_combo_DR098.tex'))

## ---------------------------
## Again recreating Table 1----
## ---------------------------

actg5202dir<-"/Users/thefeeney/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Dissertation/DissertationAnalysis/Aim 2/data/DR098/DR098_data/"

t1df<-longdf %>% select(patid,
                        Week,
                        sex,
                        age,
                        prace,
                        pethnc,
                        new_ethn,
                        strat1,
                        logrna0,
                        basecd4,
                        aidshx,
                        abc_tdf,
                        rxcode,
                        ever_dev1:ever_dev7,
                        ever_dev_doses_1:ever_dev_doses_10,
                        dose_miss_cat,
                        cens_week_1:cens_week_10,
                        target_cens_1:target_cens_10,
                        outcome_week,
                        strat1) %>%
    group_by(patid) %>% 
  mutate(cens_indic=if_else(!(patid %in% protocol_finished) & 
                              ever_dev_doses_1==0 &
                              #ever_cens==1 & 
                              sum(target_cens_1, na.rm=T)>0 ,1,0),
         no_miss=if_else(ever_dev_doses_1==0 ,1,0),
         cens_immed=if_else(target_cens_1==1 & Week==0, 1, 0,),
         one_dose=if_else(ever_dev_doses_1==1 & cens_indic==0,1,0),
         two_dose=if_else(ever_dev_doses_2==1 & cens_indic==0,1,0),
         three_dose=if_else(ever_dev_doses_3==1 & cens_indic==0,1,0),
         five_dose=if_else(ever_dev_doses_5==1 & cens_indic==0,1,0),
         ten_dose=if_else(ever_dev_doses_10==1 & cens_indic==0,1,0),
         one_dose=if_else(ever_dev_doses_1==1 & cens_indic==0,1,0),
         dose_deviation_1=case_when(sum(target_cens_1, na.rm = T)>0 & cens_week_1<outcome_week~1,
                                    TRUE~0),
         two_dose=if_else(ever_dev_doses_2==1 & cens_indic==0,1,0),
         dose_deviation_2=case_when(sum(target_cens_2, na.rm = T)>0 & cens_week_2<outcome_week~1,
                                    TRUE~0),
         three_dose=if_else(ever_dev_doses_3==1 & cens_indic==0,1,0),
         dose_deviation_3=case_when(sum(target_cens_3, na.rm = T)>0 & cens_week_3<outcome_week~1,
                                    TRUE~0),
         five_dose=if_else(ever_dev_doses_5==1 & cens_indic==0,1,0),
         dose_deviation_5=case_when(sum(target_cens_5, na.rm = T)>0 & cens_week_5<outcome_week~1,
                                    TRUE~0),
         ten_dose=if_else(ever_dev_doses_10==1 & cens_indic==0,1,0),
         dose_deviation_10=case_when(sum(target_cens_10, na.rm = T)>0 & cens_week_10<outcome_week~1,
                                     TRUE~0),
         never_miss=if_else(dose_deviation_1==0,1,0)
  ) %>% slice(1) %>% ungroup() %>% 
    mutate(rxcode=factor(rxcode, levels=c("EFV + ABC/3TC",
                                          "ATV/rtv + ABC/3TC",
                                          "EFV + TDF/FTC",
                                          "ATV/rtv + TDF/FTC")),
           newrace=fct_other(prace,keep=c("White","Black or African American","Asian"), other_level = "Other"),
           newraceth=factor(case_when(prace=="White" & pethnc=="Not Hispanic or Latino" ~"White, NH",
                            prace=="White" &  pethnc=="Hispanic or Latino" ~"White, H",
                            prace=="White" &  pethnc=="Unknown"~"White, U",
                            prace=="Black or African American" &  pethnc=="Not Hispanic or Latino"~"Black, NH",
                            prace=="Black or African American" &  pethnc=="Hispanic or Latino"~"Black, H",
                            prace=="Black or African American" &  pethnc=="Unknown"~"Black, U",
                            prace=="Asian" &  pethnc=="Not Hispanic or Latino"~"Asian, NH",
                            prace=="Asian" &  pethnc=="Hispanic or Latino"~"Asian, H",
                            prace=="Asian" &  pethnc=="Unknown"~"Asian, U",
                            prace=="Native Hawaiian or other Pacific Islander" &  pethnc=="Not Hispanic or Latino"~"NHPI, NH",
                            prace=="Native Hawaiian or other Pacific Islander" &  pethnc=="Hispanic or Latino"~"NHPI, H",
                            prace=="Native Hawaiian or other Pacific Islander" &  pethnc=="Unknown"~"NHPI, U",
                            prace=="More than One Race" &  pethnc=="Not Hispanic or Latino"~">1, NH",
                            prace=="More than One Race" &  pethnc=="Hispanic or Latino"~">1, H",
                            prace=="More than One Race" &  pethnc=="Unknown"~">1, U",
                            prace=="Unknown" &  pethnc=="Not Hispanic or Latino"~"U, NH",
                            prace=="Unknown" &  pethnc=="Hispanic or Latino"~"U, H",
                            prace=="Unknown" &  pethnc=="Unknown"~"U, U"),
                            levels=c("White, NH","White, H", "White, U",
                                               "Black, NH", "Black, H", "Black, U",
                                               "Asian, NH","Asian, H","Asian, U",
                                               "NHPI, NH", "NHPI, H", "NHPI, U",
                                               "2+, NH", "2+, H", "2+, U",
                                               "U, NH", "U, H", "U, U"),
                            ordered=TRUE),
           row_sum=select(.,ever_dev1:ever_dev7) %>% rowSums(na.rm=T),
           doses_missed=if_else(ever_dev_doses_1==1,'1','0')
           )

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Main TABLE 1 Stratified by Treatement Level
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


table1_by_txt<-tableone::CreateTableOne(
    vars=c("age",
           "sex",
           "newrace",
           "new_ethn",
           "logrna0",
           "basecd4",
           "aidshx",
           "rxcode",
           "abc_tdf",
           "strat1"),
    strata="abc_tdf",
    factorVars=c("sex",
                 "newrace",
                 "new_ethn",
                 "aidshx",
                 "rxcode",
                 "strat1"),
    data=t1df, test=FALSE, addOverall = T)

t1_by_txt_print<-print(table1_by_txt,showAllLevels = TRUE, nonnormal=c("logrna0","basecd4", "age"),
                       printToggle = F,
                       quote=F,
                       noSpaces = T); t1_by_txt_print

write.csv(t1_by_txt_print, paste0(results_tabs_dir, "table1_by_txt_DR098.csv"), row.names = TRUE)


tab1_by_txt_df<-read.csv(paste0(results_tabs_dir, "table1_by_txt_DR098.csv"))



t1_by_txt<-kableExtra::kable(tab1_by_txt_df %>% slice(1:35), format='latex',
                             booktabs=T,
                             align = "l",
                             digits=3,
                             linesep="",
                             col.names=c("",
                                         "",
                                         "Overall",
                                         "ABC/3TC",
                                         "TDF/FTC")) %>% 
    kableExtra::row_spec(c(1,3,6,7,15,16,22,25,28,30),
                         extra_latex_after = "\\addlinespace[0.1cm]")

writeLines(t1_by_txt,
           paste0(results_tabs_dir,'t1_by_txt_DR098.tex'))



## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Main Table 1 by Times Doses Missed Categorized
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


table1_strat <- function(strat = NULL) {
    filtered_data <- t1df
    
    # Apply filtering if strat is provided
    if (!is.null(strat)) {
        filtered_data <- filtered_data %>% filter(!!rlang::parse_expr(strat))
    }
    
    tableone::CreateTableOne(
        vars=c("age",
               "sex",
               "newrace",
               "new_ethn",
               "logrna0",
               "basecd4",
               "aidshx",
               "rxcode",
               "abc_tdf",
               "strat1"),
        factorVars=c("sex",
                     "newrace",
                     "new_ethn",
                     "aidshx",
                     "rxcode",
                     "strat1"),
        data = filtered_data, 
        test = FALSE, 
        addOverall = TRUE
    )
}


t1zero<-print(table1_strat(strat='never_miss==1'), nonnormal = c("age",'logrna0','basecd4'), quote=F, noSpace=T)
#t1zero<-print(table1_strat(strat='no_miss==1'),nonnormal = c("age",'logrna0','basecd4'), quote=F, noSpace=T)

t1one<-print(table1_strat(strat='dose_deviation_1==1'),nonnormal = c("age",'logrna0','basecd4'), quote=F, noSpace=T)
#t1one<-print(table1_strat(strat='one_dose==1'),nonnormal = c("age",'logrna0','basecd4'), quote=F, noSpace=T)

t1two<-print(table1_strat(strat='dose_deviation_2==1'),nonnormal = c("age",'logrna0','basecd4'), quote=F, noSpace=T)
#t1two<-print(table1_strat(strat='two_dose==1'),nonnormal = c("age",'logrna0','basecd4'), quote=F, noSpace=T)

t1three<-print(table1_strat(strat='dose_deviation_3==1'),nonnormal = c("age",'logrna0','basecd4'), quote=F, noSpace=T)
#t1three<-print(table1_strat(strat='three_dose==1'),nonnormal = c("age",'logrna0','basecd4'), quote=F, noSpace=T)

t1five<-print(table1_strat(strat='dose_deviation_5==1'),nonnormal = c("age",'logrna0','basecd4'), quote=F, noSpace=T)
#t1five<-print(table1_strat(strat='five_dose==1'),nonnormal = c("age",'logrna0','basecd4'), quote=F, noSpace=T)

t1ten<-print(table1_strat(strat='dose_deviation_10==1'),nonnormal = c("age",'logrna0','basecd4'), quote=F, noSpace=T)
#t1ten<-print(table1_strat(strat='ten_dose==1'),nonnormal = c("age",'logrna0','basecd4'), quote=F, noSpace=

table1_by_dose<-cbind(t1zero, t1one, t1two, t1three, t1five)
colnames(table1_by_dose)<-c("None missed", "One Dose", "Two Doses", "Three Doses", "Five Doses")


t1_by_dose_print<-print(table1_by_dose,showAllLevels = TRUE, nonnormal=c("logrna0","basecd4", "age"),
                        printToggle = F,
                        quote=F,
                        noSpaces = T); t1_by_dose_print
t1_10_dose_print<-print(t1ten,showAllLevels = TRUE, nonnormal=c("logrna0","basecd4", "age"),
                        printToggle = F,
                        quote=F,
                        noSpaces = T); t1_10_dose_print

write.csv(t1_by_dose_print, paste0(results_tabs_dir, "table1_by_dose_DR098.csv"), row.names = TRUE)
write.csv(t1_10_dose_print, paste0(results_tabs_dir, "table1_10_dose_DR098.csv"), row.names = TRUE)


tab1_by_dose_df<-read.csv(paste0(results_tabs_dir, "table1_by_dose_DR098.csv"))
tab1_10_dose_df<-read.csv(paste0(results_tabs_dir, "table1_10_dose_DR098.csv"))



t1_by_dose<-kableExtra::kable(tab1_by_dose_df %>% slice(1:35),'latex',
                              booktabs=T,
                              align = "l",
                              digits=3,
                              linesep="",
                              col.names=c(
                                  "",
                                  "0 missed",
                                  "1 missed",
                                  "2 missed",
                                  "3 missed",
                                  "5 missed"),) %>% 
    kableExtra::row_spec(c(1,2,3,4,5,10,25,29,30),
                         extra_latex_after = "\\addlinespace[0.1cm]") 

t1_10_dose<-kableExtra::kable(tab1_10_dose_df %>% slice(1:35),'latex',
                              booktabs=T,
                              align = "l",
                              digits=3,
                              linesep="",
                              col.names=c(
                                "",
                                "10 missed"),) %>% 
  kableExtra::row_spec(c(1,2,3,4,5,10,25,29,30),
                       extra_latex_after = "\\addlinespace[0.1cm]") 

writeLines(t1_by_dose,
           paste0(results_tabs_dir,'t1_by_dose_DR098.tex'))
writeLines(t1_10_dose,
           paste0(results_tabs_dir,'t1_10_dose_DR098.tex'))




 ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## APPENDIX TABLE 1 Stratified by Treatement Level
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


table1_by_txt<-tableone::CreateTableOne(
    vars=c("age",
           "sex",
           "logrna0",
           "basecd4",
           "newrace",
           "newraceth",
           "aidshx",
           "rxcode",
           "abc_tdf",
           "ever_cens"),
    strata="rxcode",
    factorVars=c("sex",
                 "prace",
                 "pethnc",
                 "aidshx",
                 "rxcode",
                 "ever_cens"),
    data=t1df, test=FALSE, addOverall = T)

t1_by_txt_print<-print(table1_by_txt,showAllLevels = TRUE, nonnormal=c("logrna0","basecd4", "age"),
               printToggle = F,
               quote=F,
               noSpaces = T); t1_by_txt_print

write.csv(t1_by_txt_print, "/Users/thefeeney/OneDrive - University of North Carolina at Chapel Hill/Dissertation/DissertationAnalysis/Aim 2/results/tableone_output_exp_DR098.csv", row.names = TRUE)


tab1_by_txt_df<-read.csv("/Users/timf/OneDrive - University of North Carolina at Chapel Hill/Dissertation/DissertationAnalysis/Aim 2/results/tableone_output_exp_DR098.csv")



t1_by_txt<-kableExtra::kable(tab1_by_txt_df %>% slice(1:32), format='latex',
                                  booktabs=T,
                                  align = "l",
                                  digits=3,
                           linesep="",
                           col.names=c("",
                                       "",
                                       "Overall",
                                       "+ATV",
                                       "+EFV",
                                       "+ATV",
                                       "+EFV")) %>% 
    kableExtra::row_spec(c(1,3,6,7,15,16,22,25,28,30),
                         extra_latex_after = "\\addlinespace[0.1cm]") %>% 
    kableExtra::add_header_above(c(" "=3, "ABC/3TC" = 2, "TDF/FTC" = 2))

writeLines(t1_by_txt,
           paste0(results_tabs_dir,'t1_by_txt_exp_DR098.tex'))



## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## APPENDIX Table 1 by Times Doses Missed Categorized
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



table1_strat <- function(strat = NULL, data=NULL) {
    filtered_data <- data
    
    # Apply filtering if strat is provided
    if (!is.null(strat)) {
        filtered_data <- filtered_data %>% filter(!!rlang::parse_expr(strat))
    }
    
    tableone::CreateTableOne(
        vars = c("age", "sex", "newrace","new_ethn","logrna0", "basecd4"
                 ,"aidshx", "strat1", "abc_tdf"),
        factorVars = c("sex", "prace", "pethnc", "aidshx", "strat1"),
        data = filtered_data, 
        test = FALSE, 
        addOverall = TRUE
    )
}


# Table 1 with numbers of those deviating from doses, and making sure to include only those with a censor week prior to the outcome so
# capturing those that were censored prior to outcom
t1zero<-print(table1_strat(strat='dose_deviation_1==0', data=t1df),
              quote=F, noSpace=T, showAllLevels = T, nonnormal = c("logrna0", "basecd4", "age"))
t1one<-print(table1_strat(strat='ever_dev_doses_1==1', data=t1df %>%
                            filter(ever_dev_doses_1==1 & cens_week_1<outcome_week)),
             quote=F, noSpace=T, showAllLevels=T, nonnormal = c("logrna0", "basecd4", "age"))
t1two<-print(table1_strat(strat='ever_dev_doses_2==1', data=t1df %>%
                            filter(ever_dev_doses_2==1 & cens_week_2<outcome_week)),
             quote=F, noSpace=T, showAllLevels=T, nonnormal = c("logrna0", "basecd4", "age"))
t1three<-print(table1_strat(strat='ever_dev_doses_3==1', data=t1df %>%
                              filter(ever_dev_doses_3==1 & cens_week_3<outcome_week)),
               quote=F, noSpace=T, showAllLevels=T, nonnormal = c("logrna0", "basecd4", "age"))
t1five<-print(table1_strat(strat='ever_dev_doses_5==1', data=t1df %>%
                             filter(ever_dev_doses_5==1 & cens_week_5<outcome_week)),
              quote=F, noSpace=T, showAllLevels=T, nonnormal = c("logrna0", "basecd4", "age"))
t1ten<-print(table1_strat(strat='ever_dev_doses_10==1', data=t1df %>% 
                            filter(ever_dev_doses_10==1 & cens_week_10<outcome_week)),
             quote=F, noSpace=T, showAllLevels=T,nonnormal = c("logrna0", "basecd4", "age"))

table1_by_dose<-cbind(t1zero[,1:2], t1one[,2], t1two[,2], t1three[,2], t1five[,2], t1ten[,2])
colnames(table1_by_dose)<-c("Level","None missed","One Dose Missed", "Two Doses Missed", "Three Doses Missed", "Five Doses Missed", "Ten Doses Missed")




t1_by_dose_print<-print(table1_by_dose,showAllLevels = TRUE, nonnormal=c("logrna0","basecd4", "age"),
                        printToggle = F,
                        quote=F,
                        noSpaces = T); t1_by_dose_print

write.csv(t1_by_dose_print, "/Users/thefeeney/OneDrive - University of North Carolina at Chapel Hill/Dissertation/DissertationAnalysis/Aim 2/results/table1_by_dose_exp_DR098.csv", row.names = TRUE)


tab1_by_dose_df<-read.csv("/Users/thefeeney/OneDrive - University of North Carolina at Chapel Hill/Dissertation/DissertationAnalysis/Aim 2/results/table1_by_dose_exp_DR098.csv")



t1_by_dose<-kableExtra::kable(tab1_by_dose_df %>% slice(1:35),'latex',
                              booktabs=T,
                              align = "l",
                              digits=3,
                              linesep="",
                              col.names=c(
                                          "",
                                          "0 missed",
                                          "1 missed",
                                          "2 missed",
                                          "3 missed",
                                          "5 missed",
                                          "10 missed"),) %>% 
    kableExtra::row_spec(c(1,2,3,4,5,10,25,29,30),
                         extra_latex_after = "\\addlinespace[0.1cm]") 

writeLines(t1_by_dose,
           paste0(results_tabs_dir,'t1_by_dose_exp_DR098.tex'))


