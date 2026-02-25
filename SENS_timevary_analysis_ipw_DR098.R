## :::::::::::::::::::::::::::::::::::::::
##
## Script name: Sensitivity Analysis of per protocol effects while including baseline and time-varying variables and IPW.
##
## Purpose of script:
##
## Author: Timothy Feeney
##
## Date Created: 2025-12-20
##
## Copyright (c) Timothy Feeney, 2025
## Email: tfeeney@unc.edu
##
## :::::::::::::::::::::::::::::::::::::::
##
## Notes: This is an analysis of virologic failure + death in ABC/3TC v TDF/FTC looking at a 
## range of protocol deviations while accounting for protocol deviation and censoring based on 
## baseline variables; here we us the data set created with file SENS_long_dataset_creation_DR098.R 
## Here in this sensitivity analysis we modulated resolved ties between time of deviating from 
## a missed dose and the outcome of interest by setting the deviation indicator to `1` a week prior
##
## :::::::::::::::::::::::::::::::::::::::

## set working directory for Mac and PC

# setwd("~/")          # working directory (mac)
# setwd("C:/Users/")      # working directory (PC)

## :::::::::::::::::::::::::::::::::::::::
## HEADER----
## :::::::::::::::::::::::::::::::::::::::

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation
#memory.limit(30000000)       # this is needed on some PCs to increase memory allowance, but has no impact on macs.

## :::::::::::::::::::::::::::::::::::::::

## list the packages we will need: 

packages<-c('ggplot2', 'dplyr', 'survival', 'ggsurvfit', 'here',
            'survminer', 'broom', 'kableExtra', 'gmodels', "Hmisc",
            "ggtext")      # list all the packages we need

## :::::::::::::::::::::::::::::::::::::::
## Preliminaries:----
### load up our functions into memory----

for (i in packages){
  if (!require(i, character.only = TRUE, quietly = TRUE)) { #check if package is installed
    install.packages(i)                                   # if not install it
    if (require(i, character.only = TRUE, quietly = TRUE)) { #check if install worked
      print(paste0(i, " installed"))                       # if it worked let me know
    } else {
      print(paste0(i, " failed to install"))              # if not, also let me know
    }
  } else {
    require(i, character.only=T, quietly=T)                            #otherwise if it is installed, load it
    print(paste0(i," attached"))                            #let me know it was loaded/attached
  }
  
} 

### adding function for rounding----
round_numeric <- function(x) {
  if (is.numeric(x)) {
    return(round(x, 4))  # Round to 2 decimal places
  } else {
    return(x)
  }
}

### Function for Risk extraction----

wt_model<- function(data, response, exposure, covariates,  spline_vars = NULL,
                    probs = c(0.1, 0.25, 0.5, 0.75, 0.9),
                    time_var=NULL, time_probs=c(0.10,0.33,0.5,0.66,0.9)) {
  
  df<-data
  # Helper function to create terms
  create_term <- function(var) {
    if (var %in% spline_vars) {
      # Use week_probs if the variable is 'Week', otherwise use probs
      final_probs <- if (var == time_var) time_probs else probs
      
      # Create spline term for specified covariates
      paste0(
        "splines::ns(", var, 
        ", knots = quantile(", deparse(substitute(data)), "$", var, 
        ", probs = c(", paste(final_probs, collapse = ", "), "), na.rm = TRUE))"
      )
    } else {
      # Use the covariate as-is for non-spline covariates
      var
    }
  }
  
  # Generate terms for covariates
  covariate_terms <- lapply(covariates, create_term)
  
  # Generate interaction terms
  interaction_terms <- lapply(covariate_terms, function(term) {
    paste0(exposure, " * ", term)
  })
  
  # Combine all terms into a formula
  den_formula <- as.formula(paste("(as.factor(",response,")==0)","~",
                                  paste(unlist(interaction_terms), collapse = " + ")))
  num_formula <- as.formula(paste(
    "(as.factor(", response, ") == 0) ~ ",  # Left-hand side
    exposure,                             # Exposure variable
    "* splines::ns(Week",                # Interaction with spline of Week
    ", knots = quantile(", deparse(substitute(data)), "$Week",  # Use dynamic reference for data frame
    ", probs = c(", paste(time_probs, collapse = ", "), "), na.rm = TRUE))"
  ))
  
  print(den_formula); print(num_formula)
  
  den_model<-glm(den_formula,
                 family = binomial(link='logit'),
                 data=data)
  num_model<-glm(num_formula,
                 family = binomial(link='logit'),
                 data=data)
  print(summary(den_model));print(summary(num_model))
  
  df$num_cens<-predict(num_model, newdata=df, type = 'response')
  df$den_cens<-predict(den_model, newdata=df, type = 'response')
  
  df$wt_cens <- ifelse(df[[response]]==1, 0,
                       df$num_cens/df$den_cens)
  
  final<-df %>%
    group_by(patid) %>% 
    mutate(
      !!paste0(response, "_wt") := cumprod(wt_cens),
      !!paste0(response, "_den_cens_wt") := cumprod(den_cens),
      !!paste0(response, "_num_cens_wt") := cumprod(num_cens),
      !!paste0(response, "_num_cens_wt_nonsum") := num_cens,
      !!paste0(response, "_den_cens_wt_nonsum") := den_cens,
      !!paste0(response, "_wt_nonsum") := num_cens
    ) %>% 
    ungroup()
  
  return(final)
  
}

wt_model_indic<- function(data, response, exposure, covariates,  spline_vars = NULL,
                          probs = c(0.1, 0.25, 0.5, 0.75, 0.9),
                          time_var=NULL, time_probs=c(0.10,0.33,0.5,0.66,0.9)) {
  
  df<-data
  # Helper function to create terms
  create_term <- function(var) {
    if (var %in% time_var) {
      
      paste0("factor(", var,")")
      
    } else if (var %in% spline_vars) {
      # Use week_probs if the variable is 'Week', otherwise use probs
      final_probs <- if (var == time_var) time_probs else probs
      
      # Create spline term for specified covariates
      paste0(
        "splines::ns(", var, 
        ", knots = quantile(", deparse(substitute(data)), "$", var, 
        ", probs = c(", paste(final_probs, collapse = ", "), "), na.rm = TRUE))"
      )
    } else {
      # Use the covariate as-is for non-spline covariates
      var
    }
  }
  
  # Generate terms for covariates
  covariate_terms <- lapply(covariates, create_term)
  
  # Generate interaction terms
  interaction_terms <- lapply(covariate_terms, function(term) {
    paste0(exposure, " * ", term)
  })
  
  # Combine all terms into a formula
  den_formula <- as.formula(paste("(as.factor(",response,")==0)","~",
                                  paste(unlist(interaction_terms), collapse = " + ")))
  num_formula <- as.formula(paste0(
    "(as.factor(", response, ") == 0) ~ ",  # Left-hand side
    exposure,                             # Exposure variable
    "* factor(",time_var,")")
  )
  
  print(den_formula); print(num_formula)
  
  den_model<-glm(den_formula,
                 family = binomial(link='logit'),
                 data=data)
  num_model<-glm(num_formula,
                 family = binomial(link='logit'),
                 data=data)
  print(summary(den_model));print(summary(num_model))
  
  df$num_cens<-predict(num_model, newdata=df, type = 'response')
  df$den_cens<-predict(den_model, newdata=df, type = 'response')
  
  df$wt_cens <- ifelse(df[[response]]==1, 0,
                       df$num_cens/df$den_cens)
  
  final<-df %>%
    group_by(patid) %>% 
    mutate(
      !!paste0(response, "_wt") := cumprod(wt_cens),
      !!paste0(response, "_den_cens_wt") := cumprod(den_cens),
      !!paste0(response, "_num_cens_wt") := cumprod(num_cens)
    ) %>% 
    ungroup()
  
  return(final)
  
}

weights_display<-function(df, response = NULL, save=FALSE, directory=NULL, filename=NULL) {
  
  if (is.null(response)){
    stop("Error: 'response term needed.")
  }
  
  response_columns <- c(
    paste0(response, "_num_cens_wt"),
    paste0(response, "_den_cens_wt"),
    paste0(response, "_wt")
  )
  print(response_columns)
  all(response_columns %in% colnames(df))
  # Ensure all required columns exist
  if (all(response_columns %in% colnames(df))) {
    # Compute statistics using dplyr and store results as a data.frame
    wts <- df %>% filter(if_all(dplyr::all_of(response_columns),~ . != 0)) %>% 
      dplyr::reframe(
        Statistic = c("Mean", "Max", "Min", "Median"),
        `Num Weight` = c(
          mean(.data[[response_columns[1]]], na.rm = TRUE),
          max(.data[[response_columns[1]]], na.rm = TRUE),
          min(.data[[response_columns[1]]], na.rm = TRUE),
          median(.data[[response_columns[1]]], na.rm = TRUE)
        ),
        `Denom Weight` = c(
          mean(.data[[response_columns[2]]], na.rm = TRUE),
          max(.data[[response_columns[2]]], na.rm = TRUE),
          min(.data[[response_columns[2]]], na.rm = TRUE),
          median(.data[[response_columns[2]]], na.rm = TRUE)
        ),
        `Tot Weight` = c(
          mean(.data[[response_columns[3]]], na.rm = TRUE),
          max(.data[[response_columns[3]]], na.rm = TRUE),
          min(.data[[response_columns[3]]], na.rm = TRUE),
          median(.data[[response_columns[3]]], na.rm = TRUE)
        )
      ) %>%
      as.data.frame()  # Convert to data.frame
  } else {
    stop("Some required columns are missing from the data frame.")
  }
  
  if (save) {
    if (is.null(filename)) {
      stop("Error: 'filename' required.")
    }
    if (is.null(directory)) {
      stop("Error: 'directory' required.")
    }
    
    # Create the directory if it doesn't exist
    if (!dir.exists(directory)) {
      stop("Error: 'directory' does not exist.")
    }
    
    # Save the data to the specified file
    filepath <- file.path(directory, filename)
    
    writeLines(kableExtra::kable(wts, #print out tex file
                                 'latex',
                                 booktabs=T,
                                 align = "l", #left justify
                                 digits=3),   # 3 sigfigs
               filepath)
    
    message("Data saved to: ", filepath)
  }
  
  wts_rounded<-wts %>% mutate_all(round_numeric); print(wts_rounded)
  
}

risk_xtract<-function(df, name, setype="A", save=FALSE, directory=NULL, filename=NULL){ #type A for non-itt SE, else type B
  surv_ob<-broom::tidy(df) %>% dplyr::select(time,strata,estimate,std.error, n.risk) %>%
    mutate(surv=estimate,
           risk=1-estimate, #risk is 1-survival
           time=time-1) %>% #adjust the time period by one week for each for correct time 
    select(-estimate)
  if (save) {
    if (is.null(filename)) {
      stop("Error: 'filename' required.")
    }
    if (is.null(directory)) {
      stop("Error: 'directory' required.")
    }
    
    # Create the directory if it doesn't exist
    if (!dir.exists(directory)) {
      stop("Error: 'directory' does not exist.")
    }
    
    # Save the data to the specified file
    filepath <- file.path(directory, filename)
    write.csv(surv_ob, filepath, row.names = FALSE)
    message("Data saved to: ", filepath)
  }
  
  df4896<-surv_ob %>% dplyr::filter(time==48 | time==96)
  w48_tdf_r<-df4896[3,6]; w48_abc_r<-df4896[1,6]; w48_rd<-w48_tdf_r-w48_abc_r
  w48_tdf_se<-df4896[3,3]; w48_abc_se<-df4896[1,3]; 
  
  w96_tdf_r<-df4896[4,6]; w96_abc_r<-df4896[2,6]; w96_rd<-w96_tdf_r-w96_abc_r
  w96_tdf_se<-df4896[4,3]; w96_abc_se<-df4896[2,3]; 
  
  if (setype=="A"){ 
    w48_rd_se<-sqrt((df4896[1,5]*df4896[1,6]/df4896[1,4])+
                      (df4896[3,5]*df4896[3,6]/df4896[3,4]))
    w96_rd_se<-sqrt((df4896[2,5]*df4896[2,6]/df4896[2,4])+
                      (df4896[4,5]*df4896[4,6]/df4896[4,4]))
    # 96 Week IPCW SE; Time Varying using \sqrt{\frac{p_1*(1-p_1)}{n_1}+\frac{p_2*(1-p_2)}{n_2}}
  } else {
    w48_rd_se<-sqrt(w48_tdf_se^2 + w48_abc_se^2)
    w96_rd_se<-sqrt(w96_tdf_se^2 + w96_abc_se^2)
  }
  
  df<-data.frame(Protocol=as.character(name),
                 TDF_wk48_risk=w48_tdf_r[[1]],
                 ABC_wk48_risk=w48_abc_r[[1]],
                 Wk48_RD=w48_rd[[1]],
                 Wk48_RD_SE=w48_rd_se[[1]],
                 TDF_wk96_risk=w96_tdf_r[[1]],
                 ABC_wk96_risk=w96_abc_r[[1]],
                 Wk96_RD=w96_rd[[1]],
                 Wk96_RD_SE=w96_rd_se[[1]])
  return(df)
}

weights_boxplot <- function(data = NULL,
                            weights = NULL,
                            response=NULL,
                            Weeks_to_include = c(0, 4, 8, 16, 24, 36, 48, 60, 72, 84, 96)) {
  
  # Check for required inputs
  if (is.null(data) || is.null(weights)) {
    stop("Error: 'data' and 'weights' arguments are required")
  }
  
  wts_plot<-ggplot(df_tvary %>% 
                     filter(cens_wt>0, Week<97),
                   aes(x = factor(Week), y=cens_wt)) +
    geom_boxplot(aes(middle=mean(cens_wt,na.rm=T)))+geom_hline(yintercept = 1) +
    geom_hline(yintercept = 0.99, linetype="dashed")+
    geom_hline(yintercept = 0.98, linetype="dotdash")+
    labs(title=paste0(response," Weights over time"), 
         x="Weeks", y="Censoring Weights")
  
  return(wts_plot)
}


plot_survival_curves <- function(data_list,
                                 labels,
                                 type="step",
                                 strata_labels,
                                 color_palette = NULL, 
                                 title = "Survival Curve",
                                 xlabel="Weeks",
                                 ylabel="Risk",
                                 xlim = c(0, 100),
                                 ylim = c(0, 0.15),
                                 x_breaks=c(0, 16, 24, 36, 48, 60, 72, 84, 96, 144, 180),
                                 y_breaks=c(0, 0.1, 0.2, 0.3, 0.4)) {
  if (length(data_list) != length(labels)) {
    stop("The number of datasets must match the number of labels.")
  }
  
  # Prepare color palette
  if (is.null(color_palette)) {
    color_palette <- scales::hue_pal()(length(data_list))
  }
  names(color_palette) <- labels
  
  # Combine all datasets into one, adding labels
  combined_data <- bind_rows(
    lapply(seq_along(data_list), function(i) {
      broom::tidy(data_list[[i]]) %>%
        mutate(Protocol = labels[i])
    })
  )
  
  if (type=="line"){
    # Plot the combined data
    plot <- ggplot(combined_data, aes(x = time, y = 1 - estimate, color = Protocol, linetype = strata)) +
      geom_line(linewidth = 1) +
      scale_linetype_manual(name = "Treatment\nAssignment", 
                            values = c("solid", "dotdash"), 
                            labels = strata_labels) +
      scale_color_manual(name = "Protocol", 
                         values = color_palette) +
      guides(
        linetype = guide_legend(order = 1),
        color    = guide_legend(order = 2)
      )+
      labs(color = "Protocol",
           y = ylabel,
           x = xlabel,
           title = title) +
      theme_classic() +
      coord_cartesian(ylim = ylim,
                      xlim = xlim) +
      scale_x_continuous(breaks = x_breaks,expand=c(0,0)) +
      scale_y_continuous(breaks = y_breaks,expand=c(0,0))+
      theme(text = element_text(size = 25))
  } else{
    # Plot the combined data
    plot <- ggplot(combined_data, aes(x = time, y = 1 - estimate, color = Protocol, linetype = strata)) +
      geom_step(linewidth = 1) +
      scale_linetype_manual(name = "Strata", 
                            values = c("solid", "dotdash"), 
                            labels = strata_labels) +
      scale_color_manual(name = "Protocol", 
                         values = color_palette) +
      labs(color = "Protocol",
           y = ylabel,
           x = xlabel,
           title = title) +
      theme_classic() +
      coord_cartesian(ylim = ylim,
                      xlim = xlim) +
      scale_x_continuous(breaks = x_breaks,expand=c(0,0)) +
      scale_y_continuous(breaks = y_breaks,expand=c(0,0))+
      theme(text = element_text(size = 25))
  }
  
  return(plot)
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
## :::::::::::::::::::::::::::::::::::::::
## LOAD DATA----
## :::::::::::::::::::::::::::::::::::::::

df_sens<-read.csv(paste0(data_dir,"sens_last_week_long_data2_DR098-2.csv"))


## ::::::::::::::::::::::::::::::::::::::::::::::
## Sensitivity ITT analysis
## ::::::::::::::::::::::::::::::::::::::::::::::

## :::::::::::::::::::::::::::::::::::::::
## Survival, ITT and IPCW NO LTFU----
## :::::::::::::::::::::::::::::::::::::::

sens_tvary_itt<-ggsurvfit::survfit2(Surv(time=Week, time2=end, event=outcome)~abc_tdf,
                               data=df_sens,
                               robust=T,
                               id=patid); summary(tvary_itt)

## :::::::::::::::::::::::::::::::::::::::
## RISKS for ITT----
## :::::::::::::::::::::::::::::::::::::::

# Extract Risks and save this data for later use

## ITT 

sens_itt_risks<-risk_xtract(sens_tvary_itt, name="ITT", setype="B",
                       save=T, directory=data_dir,
                       filename="sens_tvary_ittdf_DR098.csv"); print(sens_itt_risks)


## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## ONE missed does weights----
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Make the weights 

df_sens<-wt_model_indic(data=df_sens,
                         response="sens_target_cens_1",
                         exposure="abc_tdf",
                         covariates=c('n_race_col','new_ethn', 'sex', 'strat1', 'cd4', 'logrna', 'lastcd4',
                                      'lastlogrna','age','basecd4','Week' ),
                         spline_vars = c('cd4', 'logrna', 'lastcd4', 'lastlogrna','age','basecd4','Week' ),
                         time_var="Week")


## :::::::::::::::::::::::::::::::::::::::
### Weights Table ----
## :::::::::::::::::::::::::::::::::::::::

weights_display(df_sens,response="sens_target_cens_1",save=T,
                directory = results_tabs_dir, filename ='sens_timevary_weights1_sum_DR098.tex' )

## :::::::::::::::::::::::::::::::::::::::
### Survival analysis ABC v TDF  results----
## ::1:::::::::::::::::::::::::::::::::::::


sens_tvary_weighted1<-ggsurvfit::survfit2(Surv(time=Week, time2=end, event=sens_target_outcome_1)~abc_tdf,
                                     data=df_sens,
                                     weights=sens_target_cens_1_wt,
                                     robust=T,
                                     id=patid)


sens_tvary_week1_surv<-plot_survival_curves(list(sens_tvary_itt,
                                            sens_tvary_weighted1),
                                       ylim=c(0,0.18),
                                       xlim=c(0,97),
                                       type="line",
                                       title = "",
                                       y_breaks=c(.05,.1,.15),
                                       xlabel="Weeks after randomization",
                                       ylabel="Risk of Death/Viral Failure",
                                       labels=c("ITT", "1 Dose Weighted"),
                                       strata_labels=c("ABC/3TC","TDF/FTC"),
                                       color_palette = c("ITT"= "black",
                                                         "1 Dose Weighted"="gray60"));print(sens_tvary_week1_surv)

ggsave(paste0(results_figs_dir,"sens_tvary_week1_surv.png"),
       plot=sens_tvary_week1_surv, width=12, height=8, units='in')

## :::::::::::::::::::::::::::::::::::::::::
### Risk Differences ABC v TDF  results----
## :::::::::::::::::::::::::::::::::::::::::

sens_dev1_r<-risk_xtract(sens_tvary_weighted1, name = "1 Dose",
                    setype="B",
                    save=T, directory = data_dir,
                    filename = "sens_tvary_1dose_DR098.csv"); print(sens_dev1_r)


