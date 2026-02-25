## :::::::::::::::::::::::::::::::::::::::
##
## Script name: Analysis of per protocol effects while including baseline and time-varying variables and IPW.
##
## Purpose of script:
##
## Author: Timothy Feeney
##
## Date Created: 2024-04-29
##
## Copyright (c) Timothy Feeney, 2024
## Email: feeney@unc.edu
##
## :::::::::::::::::::::::::::::::::::::::
##
## Notes: This is an analysis of virologic failure + death in ABC/3TC v TDF/FTC looking at a 
## range of protocol deviations while accounting for protocol deviation and censoring based on 
## baseline variables; here we us the data set created with file long_dataset_creation.R
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
             y = "Risk",
             x = "Weeks",
             title = title) +
        theme_classic() +
        coord_cartesian(ylim = ylim,
                        xlim = xlim) +
        scale_x_continuous(breaks = x_breaks,expand=c(0,0)) +
        scale_y_continuous(breaks = y_breaks,expand=c(0,0))+
        theme(text = element_text(size = 25))
    } else if (type=="linecomb"){
      # Plot the combined data
      plot <- ggplot(combined_data, aes(x = time, y = 1 - estimate, color = Protocol,
                                        #linetype = strata
                                        )) +
        geom_line(linewidth = 1) +
        #scale_linetype_manual(name = "Treatment\nAssignment", 
        #                      values = c("solid", "dotdash"), 
        #                      labels = strata_labels) +
        scale_color_manual(name = "Protocol", 
                           values = color_palette) +
        #guides(
        #  linetype = guide_legend(order = 1),
        #  color    = guide_legend(order = 2)
        #)+
        labs(color = "Protocol",
             y = "Risk",
             x = "Weeks",
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
             y = "Risk",
             x = "Weeks",
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

df<-read.csv(paste0(data_dir,"last_week_long_data2_DR098-2.csv"))
df_tvary<-df %>% mutate(end=Week+1,
                        outcome_old=outcome,
                        outcome=if_else(outcome_week+1<end, NA, outcome),
                        n_race_col=case_when(new_race=="White"~"White",
                                             new_race=="Black/AA"~"AA",
                                             new_race=="Other" |
                                                 new_race==">1" |
                                                 new_race=="Asian"~"Other",
                                             new_race=="Unknown"~"Unknown",
                                             TRUE~NA))

## DATA DESCRIPTION----
## :::::::::::::::::::::::::::::::::::::::

# see Baseline_analysis_ipw.R for descriptive analysis

### Survival analysis of censoring ----
### this will allow us to see how the risk of censoring changes over time

cens_surv<-ggsurvfit::survfit2(Surv(time=Week,time2=end, event=cens)~abc_tdf,
                               data=df_tvary,
                               robust=T,
                               id=patid); summary(cens_surv)

cens_surv_agg<-ggsurvfit::survfit2(Surv(time=Week,time2=end, event=cens)~1,
                               data=df_tvary,
                               robust=T,
                               id=patid); summary(cens_surv_agg)


miss1_surv<-ggsurvfit::survfit2(Surv(time=Week, time2=end, event=target_cens_1)~abc_tdf,
                                        data=df_tvary,
                                        robust=T,
                                        id=patid); summary(miss1_surv)

miss1_surv_agg<-ggsurvfit::survfit2(Surv(time=Week, time2=end, event=target_cens_1)~1,
                                data=df_tvary,
                                robust=T,
                                id=patid); summary(miss1_surv_agg)

miss2_surv<-ggsurvfit::survfit2(Surv(time=Week, time2=end, event=target_cens_2)~abc_tdf,
                                data=df_tvary,
                                robust=T,
                                id=patid); summary(miss2_surv)

miss2_surv_agg<-ggsurvfit::survfit2(Surv(time=Week, time2=end, event=target_cens_2)~1,
                                data=df_tvary,
                                robust=T,
                                id=patid); summary(miss2_surv_agg)

miss10_surv<-ggsurvfit::survfit2(Surv(time=Week,time2=end, event=target_cens_10)~abc_tdf,
                                data=df_tvary,
                                robust=T,
                                id=patid); summary(miss10_surv)

miss10_surv_agg<-ggsurvfit::survfit2(Surv(time=Week,time2=end, event=target_cens_10)~1,
                                 data=df_tvary,
                                 robust=T,
                                 id=patid); summary(miss10_surv_agg)

tvary_censoring_surv<-plot_survival_curves(list(cens_surv,miss1_surv, miss2_surv, miss10_surv),
                     x_breaks=c(0, 16, 24, 36, 48, 60, 72, 84, 96, 144, 180),
                     y_breaks=c(0.05, 0.15, 0.25, 0.35),
                     type="line",
                     labels=c("LTFU", "1 Dose Protocol", "2 Dose Protocol","10 Dose Protocol"),
                     strata_labels=c("ABC/3TC","TDF/FTC"),
                     #color_palette = c("LTFU" = "black",
                     #                  "1 Dose Protocol" = "#DF8F44FF",
                     #                  "2 Dose Protocol" = "blue",
                     #                  "10 Dose Protocol" = "#AD002AFF"
                     #                  ),
                                       title = "",
                     xlim = c(0,97),
                     ylim=c(0,.36))+
                     xlab("Weeks Post Randomization")+
                     ylab("Risk of Censoring")+
  theme(legend.position = c(0.16,0.78),
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1))+
  scale_color_manual(
    values = c(
      "LTFU" = "black",
      "1 Dose Protocol" = "#DF8F44FF",
      "2 Dose Protocol" = "blue",
      "10 Dose Protocol" = "#AD002AFF"
    ),
    limits = c("LTFU", "1 Dose Protocol", "2 Dose Protocol", "10 Dose Protocol") # <-- force order
  ); print(tvary_censoring_surv)

tvary_censoring_surv_agg<-plot_survival_curves(list(cens_surv_agg,miss1_surv_agg, miss2_surv_agg, miss10_surv_agg),
                                           x_breaks=c(0, 16, 24, 36, 48, 60, 72, 84, 96, 144, 180),
                                           y_breaks=c(0.05, 0.15, 0.25, 0.35),
                                           type="linecomb",
                                           labels=c("LTFU", "1 Dose Protocol", "2 Dose Protocol","10 Dose Protocol"),
                                           strata_labels=c("ABC/3TC","TDF/FTC"),
                                           #color_palette = c("LTFU" = "black",
                                           #                  "1 Dose Protocol" = "#DF8F44FF",
                                           #                  "2 Dose Protocol" = "blue",
                                           #                  "10 Dose Protocol" = "#AD002AFF"
                                           #                  ),
                                           title = "",
                                           xlim = c(0,97),
                                           ylim=c(0,.36))+
  xlab("Weeks Post Randomization")+
  ylab("Risk of Censoring")+
  theme(legend.position = c(0.16,0.78),
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1))+
  scale_color_manual(
    values = c(
      "LTFU" = "black",
      "1 Dose Protocol" = "#DF8F44FF",
      "2 Dose Protocol" = "blue",
      "10 Dose Protocol" = "#AD002AFF"
    ),
    limits = c("LTFU", "1 Dose Protocol", "2 Dose Protocol", "10 Dose Protocol") # <-- force order
  ); print(tvary_censoring_surv_agg)
                    

ggsave(paste0(results_figs_dir,"tvary_censor_surv_DR098.png"),
       plot=tvary_censoring_surv, width=12, height=8, units='in')

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::
## MAKING IPCW NO LTFU----
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::


# Make the weights 

#### OLD WEIGHTS CODE ####

tvary_den_cens_model<-glm((as.factor(cens)==0) ~ abc_tdf*n_race_col+
                              abc_tdf*sex+
                              abc_tdf*strat1 +
                              abc_tdf*splines::ns(cd4,
                                                  knots=quantile(df_tvary$cd4,
                                                                 probs=c(0.1,0.25,0.5,0.75,0.9),
                                                                 na.rm = T)) +
                              abc_tdf*splines::ns(logrna,
                                                  knots=quantile(df_tvary$logrna,
                                                                 probs=c(0.1,0.25,0.5,0.75,0.9),
                                                                 na.rm=T)) +
                              abc_tdf*splines::ns(lastcd4,
                                                  knots=quantile(df_tvary$lastcd4,
                                                                 probs=c(0.1,0.25,0.5,0.75,0.9),
                                                                 na.rm=T))+
                              abc_tdf*splines::ns(lastlogrna,
                                                  knots=quantile(df_tvary$lastlogrna,
                                                                 probs=c(0.1,0.25,0.5,0.75,0.9),
                                                                 na.rm=T)) +
                              abc_tdf*splines::ns(age,
                                                  knots=quantile(df_tvary$age,
                                                                 probs=c(0.1,0.25,0.5,0.75,0.9))) +
                              abc_tdf*splines::ns(basecd4,
                                                  knots=quantile(df_tvary$basecd4,
                                                                 probs=c(0.1,0.25,0.5,0.75,0.9))) +
                              abc_tdf*factor(Week),
                          family = binomial(link='logit'),
                          data=df_tvary)

tvary_num_cens_model<-glm((as.factor(cens)==0) ~ 
                              abc_tdf*strat1+abc_tdf*factor(Week),
                          family = binomial(link='logit'),
                          data=df_tvary)

df_tvary$num_cens<-predict(tvary_num_cens_model, newdata=df_tvary, type = 'response')
df_tvary$den_cens<-predict(tvary_den_cens_model, newdata=df_tvary, type = 'response')

df_tvary$wt_cens <- ifelse(df_tvary$cens==1, 0,
                           df_tvary$num_cens/df_tvary$den_cens)

df_tvary<-df_tvary %>%
    group_by(patid) %>% 
    mutate(cens_wt=cumprod(wt_cens),
           den_cens_wt=cumprod(den_cens),
           num_cens_wt=cumprod(num_cens)) %>% 
    ungroup()

#### END ####

#this function recreates the analysis directly above; note the spline on weeks has slightly different knot cutoffs

df_tvary<-wt_model_indic(data=df_tvary,
                   response="cens",
                   exposure="abc_tdf",
                   covariates=c('n_race_col','new_ethn', 'sex', 'strat1', 'cd4', 'logrna',
                                'lastcd4', 'lastlogrna','age','basecd4','Week' ),
                   spline_vars = c('cd4', 'logrna', 'lastcd4', 'lastlogrna','age',
                                   'basecd4','Week' ),
                   time_var="Week")

## :::::::::::::::::::::::::::::::::::::::
### IPCW Weights Table ----
## :::::::::::::::::::::::::::::::::::::::

weights_display(df_tvary,response='cens', save = T,
                directory = results_tabs_dir, filename = 'timevary_cens_weights_DR098.tex')

weights_boxplot(data=df_tvary,
                weights='cens_wt',
                Weeks_to_include = c(0,4,8,16,24,36,48,60,72,84,96))+
    geom_hline(yintercept = 0.985,
               linetype="dotdash")

## :::::::::::::::::::::::::::::::::::::::
## Survival, ITT and IPCW NO LTFU----
## :::::::::::::::::::::::::::::::::::::::

tvary_itt<-ggsurvfit::survfit2(Surv(time=Week, time2=end, event=outcome)~abc_tdf,
                               data=df_tvary,
                               robust=T,
                               id=patid); summary(tvary_itt)


tvary_ipcw<-ggsurvfit::survfit2(Surv(time=Week, time2=end, event=outcome)~abc_tdf,
                                data=df_tvary,
                                weights=cens_wt,
                                robust=T,
                                id=patid)

## :::::::::::::::::::::::::::::::::::::::
## RISKS for ITT And IPCW----
## :::::::::::::::::::::::::::::::::::::::

# Extract Risks and save this data for later use

## ITT 

itt_risks<-risk_xtract(tvary_itt, name="ITT", setype="B",
                       save=T, directory=data_dir,
                       filename="tvary_ittdf_DR098.csv"); print(itt_risks)

## IPCW

ipcw_risks<-risk_xtract(tvary_ipcw, name ="IPCW", setype="A",
                        save=T, directory = data_dir,
                        filename = "tvary_ipcwdf_DR098.csv"); print(ipcw_risks)

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## ONE missed does weights----
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Make the weights 

df_tvary<-wt_model_indic(data=df_tvary,
                   response="target_cens_1",
                 exposure="abc_tdf",
                 covariates=c('n_race_col','new_ethn', 'sex', 'strat1', 'cd4', 'logrna', 'lastcd4',
                              'lastlogrna','age','basecd4','Week' ),
                 spline_vars = c('cd4', 'logrna', 'lastcd4', 'lastlogrna','age','basecd4','Week' ),
                 time_var="Week")


## :::::::::::::::::::::::::::::::::::::::
### Weights Table ----
## :::::::::::::::::::::::::::::::::::::::

weights_display(df_tvary,response="target_cens_1",save=T,
                directory = results_tabs_dir, filename ='timevary_weights1_sum_DR098.tex' )

miss1_wt_bxplt<-weights_boxplot(data=df_tvary,
                                response='target_cens_1_wt',
                                weights='target_cens_1_wt',
                                Weeks_to_include = c(0,4,8,16,24,36,48,60,72,84,96));
print(miss1_wt_bxplt)


ggsave(paste0(results_figs_dir,"miss1_wt_bxplt.png"),
       plot=miss1_wt_bxplt, width=12, height=5, units='in')
## :::::::::::::::::::::::::::::::::::::::
### Survival analysis ABC v TDF  results----
## ::1:::::::::::::::::::::::::::::::::::::


tvary_weighted1<-ggsurvfit::survfit2(Surv(time=Week, time2=end, event=target_outcome_1)~abc_tdf,
                                     data=df_tvary,
                                     weights=target_cens_1_wt,
                                     robust=T,
                                     id=patid)


tvary_week1_surv<-plot_survival_curves(list(tvary_itt,
                                            tvary_ipcw,
                                            tvary_weighted1),
                                       ylim=c(0,0.18),
                                       xlim=c(0,97),
                                       type="line",
                                       labels=c("ITT", "IPCW", "1 Dose Weighted"),
                                       strata_labels=c("ABC","TDF"),
                                       color_palette = c("ITT"= "#AD002AFF",
                                                         "LTFU"="#925E9FFF",
                                                         "1 Dose Weighted"="#DF8F44FF"));print(tvary_week1_surv)

ggsave(paste0(results_figs_dir,"tvary_week1_surv.png"),
       plot=tvary_week1_surv, width=12, height=5, units='in')

## :::::::::::::::::::::::::::::::::::::::::
### Risk Differences ABC v TDF  results----
## :::::::::::::::::::::::::::::::::::::::::

dev1_r<-risk_xtract(tvary_weighted1, name = "1 Dose",
                    setype="B",
                    save=T, directory = data_dir,
                    filename = "tvary_1dose_DR098.csv"); print(dev1_r)


## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## TWO DOSE deviation weights----
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

df_tvary<-wt_model(data=df_tvary,response="target_cens_2",
                  exposure="abc_tdf",
                  covariates=c('n_race_col','new_ethn', 'sex', 'strat1', 'cd4',
                               'logrna', 'lastcd4', 'lastlogrna','age','basecd4','Week' ),
                  spline_vars = c('cd4', 'logrna', 'lastcd4', 'lastlogrna',
                                  'age','basecd4','Week' ),
                  time_var="Week")


## :::::::::::::::::::::::::::::::::::::::
### Weights Table ----
## :::::::::::::::::::::::::::::::::::::::

weights_display(df_tvary,response="target_cens_2",save=T, 
                directory = results_tabs_dir, filename ='timevary_weights2_sum_DR098.tex' )

weights_boxplot(data=df_tvary, weights='target_cens_2_wt', Weeks_to_include = c(0,4,8,16,24,36,48,60,72,84,96))


## :::::::::::::::::::::::::::::::::::::::::
### Survival analysis ABC v TDF  results----
## :::::::::::::::::::::::::::::::::::::::::

tvary_weighted2<-ggsurvfit::survfit2(Surv(time=Week, time2=end, event=target_outcome_2)~abc_tdf,
                                     data=df_tvary,
                                     weights=target_cens_2_wt,
                                     robust=T,
                                     id=patid)

tvary_week2_surv<-plot_survival_curves(list(tvary_itt, tvary_ipcw, tvary_weighted2),
                                       y_breaks=c(.05,.1,.15),
                                       ylim=c(0,.18),
                                       type="line",
                                       labels=c("ITT", "IPCW", "2 Dose Weighted"),
                                       strata_labels=c("ABC","TDF"),
                                       color_palette = c("ITT"= "#AD002AFF",
                                                         "LTFU"="#925E9FFF",
                                                         "2 Dose Weighted"="#DF8F44FF"),
                                       title = "2 Doses Missed");print(tvary_week2_surv)

ggsave(paste0(results_figs_dir,"tvary_week2_surv.png"),
       plot=tvary_week2_surv, width=12, height=5, units='in')

## :::::::::::::::::::::::::::::::::::::::::
### Risk Differences ABC v TDF  results----
## :::::::::::::::::::::::::::::::::::::::::

dev2_r<-risk_xtract(tvary_weighted2, name = "2 Doses",
                    save=T, directory = data_dir,
                    filename = "tvary_2dose_DR098.csv"); print(dev2_r)

## :::::::::::::::::::::::::::::::::::::::
## THREE DOSE deviation weights----
## :::::::::::::::::::::::::::::::::::::::

df_tvary<-wt_model(data=df_tvary,response="target_cens_3",
                  exposure="abc_tdf",
                  covariates=c('n_race_col','new_ethn', 'sex', 'strat1', 'cd4', 'logrna',
                               'lastcd4', 'lastlogrna','age','basecd4','baserna','Week' ),
                  spline_vars = c('cd4', 'logrna', 'lastcd4', 'lastlogrna','age',
                                  'basecd4','baserna','Week' ),
                  time_var="Week")


## :::::::::::::::::::::::::::::::::::::::
### Weights Table ----
## :::::::::::::::::::::::::::::::::::::::

weights_display(df_tvary,response="target_cens_3",
                save=T, directory = results_tabs_dir,
                filename ='timevary_weights3_sum_DR098.tex' )

weights_boxplot(data=df_tvary,
                weights='target_cens_3_wt',
                Weeks_to_include = c(0, 16, 24, 36, 48, 60, 72, 84, 96, 144, 180))
## :::::::::::::::::::::::::::::::::::::::::
### Survival analysis ABC v TDF  results----
## :::::::::::::::::::::::::::::::::::::::::

tvary_weighted3<-ggsurvfit::survfit2(Surv(time=Week, time2=end, event=target_outcome_3)~abc_tdf,
                                     data=df_tvary,
                                     weights=target_cens_3_wt,
                                     robust=T,
                                     id=patid)
tvary_week3_surv<-plot_survival_curves(list(tvary_itt, tvary_ipcw, tvary_weighted3),
                                       y_breaks=c(.05,.1,.15),
                                       ylim=c(0,.18),
                                       labels=c("ITT", "IPCW", "2 Dose Weighted"),
                                       strata_labels=c("ABC","TDF"),
                                       color_palette = c("ITT"= "#AD002AFF",
                                                         "LTFU"="#925E9FFF",
                                                         "3 Dose Weighted"="#DF8F44FF"),
                                       title = "3 Doses Missed"); print(tvary_week3_surv)

ggsave(paste0(results_figs_dir,"tvary_week3_surv.png"),
       plot=tvary_week3_surv, width=12, height=5, units='in')

## :::::::::::::::::::::::::::::::::::::::::
### Risk Differences ABC v TDF  results----
## :::::::::::::::::::::::::::::::::::::::::
dev3_r<-risk_xtract(tvary_weighted3, name = "3 Doses",
                    save=T, directory = data_dir,
                    filename = "tvary_3dose_DR098.csv"); print(dev3_r)

## :::::::::::::::::::::::::::::::::::::::
## FIVE DOSE deviation weights----
## :::::::::::::::::::::::::::::::::::::::

df_tvary<-wt_model(data=df_tvary,response="target_cens_5",
                   exposure="abc_tdf",
                   covariates=c('n_race_col','new_ethn', 'sex', 'strat1', 'cd4', 'logrna',
                                'lastcd4', 'lastlogrna','age','basecd4','baserna','Week' ),
                   spline_vars = c('cd4', 'logrna', 'lastcd4', 'lastlogrna','age',
                                   'basecd4','baserna','Week' ),
                   time_var="Week")


## :::::::::::::::::::::::::::::::::::::::
### Weights Table ----
## :::::::::::::::::::::::::::::::::::::::

weights_display(df_tvary,response="target_cens_5",
                save=T, directory = results_tabs_dir,
                filename ='timevary_weights5_sum_DR098.tex' )

weights_boxplot(data=df_tvary,
                weights='target_cens_5_wt',
                Weeks_to_include = c(0, 16, 24, 36, 48, 60, 72, 84, 96, 144, 180))
## :::::::::::::::::::::::::::::::::::::::::
### Survival analysis ABC v TDF  results----
## :::::::::::::::::::::::::::::::::::::::::

tvary_weighted5<-ggsurvfit::survfit2(Surv(time=Week, time2=end, event=target_outcome_5)~abc_tdf,
                                     data=df_tvary,
                                     weights=target_cens_5_wt,
                                     robust=T,
                                     id=patid)
tvary_week5_surv<-plot_survival_curves(list(tvary_itt, tvary_ipcw, tvary_weighted5),
                                       y_breaks=c(.05,.1,.15),
                                       ylim=c(0,.18),
                                       labels=c("ITT", "IPCW", "2 Dose Weighted"),
                                       strata_labels=c("ABC","TDF"),
                                       color_palette = c("ITT"= "#AD002AFF",
                                                         "LTFU"="#925E9FFF",
                                                         "5 Dose Weighted"="#DF8F44FF"),
                                       title = "5 Doses Missed"); print(tvary_week5_surv)

ggsave(paste0(results_figs_dir,"tvary_week5_surv.png"),
       plot=tvary_week5_surv, width=12, height=5, units='in')

## :::::::::::::::::::::::::::::::::::::::::
### Risk Differences ABC v TDF  results----
## :::::::::::::::::::::::::::::::::::::::::
dev5_r<-risk_xtract(tvary_weighted5, name = "5 Doses",
                    save=T, directory = data_dir,
                    filename = "tvary_5dose_DR098.csv"); print(dev5_r)

## :::::::::::::::::::::::::::::::::::::::
## TEN DOSE deviation weights----
## :::::::::::::::::::::::::::::::::::::::

df_tvary<-wt_model(data=df_tvary,response="target_cens_10",
                   exposure="abc_tdf",
                   covariates=c('n_race_col','new_ethn', 'sex', 'strat1', 'cd4', 'logrna',
                                'lastcd4', 'lastlogrna','age','basecd4','baserna','Week' ),
                   spline_vars = c('cd4', 'logrna', 'lastcd4', 'lastlogrna','age',
                                   'basecd4','baserna','Week' ),
                   time_var="Week")


## :::::::::::::::::::::::::::::::::::::::
### Weights Table ----
## :::::::::::::::::::::::::::::::::::::::

weights_display(df_tvary,response="target_cens_10",
                save=T, directory = results_tabs_dir,
                filename ='timevary_weights10_sum_DR098.tex' )

weights_boxplot(data=df_tvary,
                weights='target_cens_10_wt',
                Weeks_to_include = c(0, 16, 24, 36, 48, 60, 72, 84, 96, 144, 180))
## :::::::::::::::::::::::::::::::::::::::::
### Survival analysis ABC v TDF  results----
## :::::::::::::::::::::::::::::::::::::::::

tvary_weighted10<-ggsurvfit::survfit2(Surv(time=Week, time2=end, event=target_outcome_10)~abc_tdf,
                                     data=df_tvary,
                                     weights=target_cens_10_wt,
                                     robust=T,
                                     id=patid)
tvary_week10_surv<-plot_survival_curves(list(tvary_itt, tvary_ipcw, tvary_weighted10),
                                        y_breaks=c(.05,.1,.15),
                                        ylim=c(0,.18),
                                       labels=c("ITT", "IPCW", "10 Dose Weighted"),
                                       strata_labels=c("ABC","TDF"),
                                       color_palette = c("ITT"= "#AD002AFF",
                                                         "LTFU"="#925E9FFF",
                                                         "10 Dose Weighted"="#DF8F44FF"),
                                       title = "10 Doses Missed"); print(tvary_week10_surv)

ggsave(paste0(results_figs_dir,"tvary_week10_surv.png"),
       plot=tvary_week10_surv, width=12, height=5, units='in')

## :::::::::::::::::::::::::::::::::::::::::
### Risk Differences ABC v TDF  results----
## :::::::::::::::::::::::::::::::::::::::::
dev10_r<-risk_xtract(tvary_weighted10, name = "10 Doses",
                     save=T, directory = data_dir,
                     filename = "tvary_10dose_DR098.csv"); print(dev10_r)


## :::::::::::::::::::::::::::::::::::::::
## Aggregate graph----
## :::::::::::::::::::::::::::::::::::::::
labels<-factor(c("ITT",
                 #"LTFU Only", 
                 "10 Dose Protocol",
                 "3 Dose Protocol",
                 "2 Dose Protocol",
                 "1 Dose Protocol"),
               levels = c("ITT",
                          #"LTFU Only",
                          "10 Dose Protocol",
                          "3 Dose Protocol",
                          "2 Dose Protocol",
                          "1 Dose Protocol"))


tvary_agg_base_surv<-plot_survival_curves(list(tvary_itt,
                                               #tvary_ipcw,
                                               tvary_weighted10,
                                               tvary_weighted3,
                                               tvary_weighted2,
                                               tvary_weighted1),
                                          y_breaks=c(.05,.1,.15),
                     labels=labels,
                     strata_labels=c("ABC/3TC","TDF/FTC"),
                     color_palette = c("ITT"= "black",
                                       #"LTFU"="#D31A38",
                                       "10 Dose Protocol"="forestgreen",
                                       "3 Dose Protocol"="blue",
                                       "2 Dose Protocol"="#925E9FFF",
                                       "1 Dose Protocol"="#DF8F44FF"
                                       ),
                     type="line",
                     title="",
                     ylim = c(0,0.18))+
                    xlab("Weeks from Randomization")+
                    ylab("Risk of Death/Viral Failure")+
  theme(legend.box.background = element_rect(),
        legend.box.margin = margin(.5,.5,.5,.5),
        #legend.text = element_text(size = 10),
        #legend.title = element_text(size=12,face = "bold"),
        legend.position = c(0.14, 0.78));print(tvary_agg_base_surv)

 ggsave(paste0(results_figs_dir,"tvary_agg_surv_DR098.png"),
       plot=tvary_agg_base_surv, width=12, height=8, units='in')
 
 
 ## :::::::::::::::::::::::::::::::::::::::
 ## Aggregate graph WITH IPCW curve for supplement----
 ## :::::::::::::::::::::::::::::::::::::::
 labels<-factor(c("ITT",
                  "LTFU Only", 
                  "10 Dose Protocol",
                  "3 Dose Protocol",
                  "2 Dose Protocol",
                  "1 Dose Protocol"),
                levels = c("ITT",
                           "LTFU Only",
                           "10 Dose Protocol",
                           "3 Dose Protocol",
                           "2 Dose Protocol",
                           "1 Dose Protocol"))
 
 
 tvary_agg_base_surv<-plot_survival_curves(list(tvary_itt,
                                                tvary_ipcw,
                                                tvary_weighted10,
                                                tvary_weighted3,
                                                tvary_weighted2,
                                                tvary_weighted1),
                                           y_breaks=c(.05,.1,.15),
                                           labels=labels,
                                           strata_labels=c("ABC","TDF"),
                                           color_palette = c("ITT"= "black",
                                                             "LTFU"="#D31A38",
                                                             "10 Dose Protocol"="forestgreen",
                                                             "3 Dose Protocol"="blue",
                                                             "2 Dose Protocol"="#925E9FFF",
                                                             "1 Dose Protocol"="#DF8F44FF"
                                           ),
                                           type="line",
                                           title="",
                                           ylim = c(0,0.18))+
   xlab("Weeks from Randomization")+
   ylab("Risk of Death/Viral Failure")+
   theme(legend.box.background = element_rect(),
         legend.box.margin = margin(.5,.5,.5,.5),
         #legend.text = element_text(size = 10),
         #legend.title = element_text(size=12,face = "bold"),
         legend.position = c(0.14, 0.78));print(tvary_agg_base_surv)
 
 ggsave(paste0(results_figs_dir,"tvary_agg_surv_DR098_supp.png"),
        plot=tvary_agg_base_surv, width=12, height=8, units='in')

 ## :::::::::::::::::::::::::::::::::::::::
## Risk Table----
## :::::::::::::::::::::::::::::::::::::::

tvary_RDs_cbind<-rbind(itt_risks,#ipcw_risks,
                 dev10_r,
                 dev5_r,
                 dev3_r,
                 dev2_r,
                 dev1_r);print(tvary_RDs)

tvary_RDs<-tvary_RDs_cbind %>% mutate(wk48ci=paste0(round(Wk48_RD-1.96*Wk48_RD_SE,4)," ,",round(Wk48_RD+1.96*Wk48_RD_SE,4) ),
                                      wk96ci=paste0(round(Wk96_RD-1.96*Wk96_RD_SE,4)," ,", round(Wk96_RD+1.96*Wk96_RD_SE,4))) %>% 
  select(Protocol,TDF_wk48_risk, ABC_wk48_risk,  Wk48_RD, Wk48_RD_SE,wk48ci, TDF_wk96_risk, ABC_wk96_risk, Wk96_RD, Wk96_RD_SE, wk96ci); print(tvary_RDs)

tvary_RDs %>% mutate_all(round_numeric)

tvary_RD_table<-kableExtra::kable(tvary_RDs,
                                  col.names=c("Protocol",
                                              "TDF Risk",
                                              "ABC Risk",
                                              "Risk Diff",
                                              "RD SE",
                                              "ABC Risk",
                                              "TDF Risk",
                                              "Risk Diff",
                                              "RD SE"),
                                  'latex',
                                  booktabs=T,
                                  align = "l",
                                  digits=3,
                                  linesep="") %>%
    row_spec(2, extra_latex_after = "\\addlinespace") %>% 
    kableExtra::add_header_above(c(" "=1, "48 Weeks" = 4, "96 Weeks" = 4))

writeLines(tvary_RD_table,
           paste0(results_tabs_dir,'tvary_risk_diff_table_DR098.tex'))

tv_Wk48_forest_plot<-ggplot(data=tvary_RDs, aes(y=factor(Protocol,
                                                         levels=c("1 Dose", "2 Doses","3 Doses","5 Doses","10 Doses","IPCW","ITT")),
                           x=Wk48_RD, xmin=Wk48_RD+1.96*Wk48_RD_SE, xmax=Wk48_RD-1.96*Wk48_RD_SE))+
    geom_point(size=5)+ #change size of point
    geom_errorbarh(height=0.2,linewidth=1)+#chnage size of lines, size for thickness
    geom_vline(xintercept=0, color='black', linetype='dashed', alpha=0.5)+ #null line
    theme_classic()+ #removes the garbage theme components
    labs(x="Risk Difference",
         y="Number of Doses Missed",
         title = "TDF/FTC v ABC/3TC at 48 Weeks")+
    geom_text(aes(x = 0.02,#text position
                  label = paste0(sprintf("%.3f",signif(Wk96_RD, 2)), "; ", sprintf("%.3f",signif(Wk96_RD_SE, 2)))), 
              hjust = 0, size = 8) +
    annotate("richtext", label="<b>Estimate; SE</b>", x=0.038, y=7.4, size=8)+
    coord_cartesian(xlim = c(with(tvary_RDs, min(Wk48_RD+1.96*Wk48_RD_SE))-.08, 0.08))+
    theme(text=element_text(size=25)); print(tv_Wk48_forest_plot)

tv_Wk96_forest_plot<-ggplot(data=tvary_RDs, aes(y=factor(Protocol,
                                                         levels=c("1 Dose", "2 Doses","3 Doses","5 Doses","10 Doses","IPCW","ITT")),
                                             x=Wk96_RD, xmin=Wk96_RD+1.96*Wk96_RD_SE, xmax=Wk96_RD-1.96*Wk96_RD_SE))+
    geom_point(size=3)+ #change size of point
    geom_errorbarh(height=0.2,linewidth=1)+#chnage size of lines, size for thickness
    geom_vline(xintercept=0, color='black', linetype='dashed', alpha=0.5)+ #null line
    theme_classic()+ #removes the garbage theme components
    labs(x="Risk Difference",
         y="Number of Doses Missed",
         title = "TDF/FTC v ABC/3TC at 96 Weeks")+
    geom_text(aes(x = 0.02, #text position 
                  label = paste0(sprintf("%.3f",signif(Wk96_RD, 2)), "; ", sprintf("%.3f",signif(Wk96_RD_SE, 2)))), 
              hjust = 0, size = 8) +
    annotate("richtext", label="<b>Estimate; SE</b>", x=0.038, y=7.4, size=8)+
    coord_cartesian(xlim = c(with(tvary_RDs, min(Wk96_RD+1.96*Wk96_RD_SE))-.08, 0.08))+
    theme(text=element_text(size=25)); print(tv_Wk96_forest_plot)

ggsave(paste0(results_figs_dir,"tvary_wk48_forest_DR098.png"),
       plot=tv_Wk48_forest_plot, width=12, height=9, units='in')

ggsave(paste0(results_figs_dir,"tvary_wk96_forest_DR098.png"),
       plot=tv_Wk96_forest_plot, width=12, height=9, units='in')

## :::::::::::::::::::::::::::::::::::::::
## Aggregate graph----
## Graph combining baseline and time varying adjustment for censoring and protocol
## Deviation. Here We look at 1 does compared to IPCW and ITT
## :::::::::::::::::::::::::::::::::::::::


## :::::::::::::::::::::::::::::::::::::::
## survival curves for SER 2024 presentation ----
## :::::::::::::::::::::::::::::::::::::::

ser_surv1<-ggplot() +
    geom_step(data = tvary_ipcwdf %>% filter(strata=="abc_tdf=ABC/3TC"),
              aes(x = time, y = risk,
                  #linetype = strata,
                  color = "IPCW"), linewidth=1) +
    geom_step(data = tvary_wdf1 %>% filter(strata=="abc_tdf=ABC/3TC"),
              aes(x = time, y = risk,
                  #linetype = strata,
                  color = "1 doses"),linewidth=1)+
    geom_step(data = tvary_wdf2 %>% filter(strata=="abc_tdf=ABC/3TC"),
              aes(x = time, y = risk,
                  #linetype = strata,
                  color = "2 doses"),linewidth=1)+
    geom_step(data = tvary_wdf3 %>% filter(strata=="abc_tdf=ABC/3TC"),
              aes(x = time, y = risk,
                  #linetype = strata,
                  color = "3 doses")) +
    geom_step(data = tvary_ittdf %>% filter(strata=="abc_tdf=ABC/3TC"),
              aes(x = time, y = risk,
                  #linetype = strata,
                  color = "ITT"),linewidth=1)+
    #scale_linetype_manual(name = "Strata", 
    #                      values = c("solid", "dotdash"), 
    #                      labels = c("ABC")) +
    scale_color_manual(name = "Protocol", 
                       values = c("ITT"="#925E9FFF","IPCW" = "black", "1 doses" = "#DF8F44FF",
                                  "2 doses"= "darkblue" , "3 doses" = "#AD002AFF"),
                       breaks=c( "ITT","IPCW",  "3 doses" ,"2 doses", "1 doses")) +
    labs(color = "Protocol",
         y="Risk",
         x="Weeks",
         title="Risk of Composite Outcome; ABC/3TC arm")+
    theme_classic()+
    coord_cartesian(ylim=c(0,0.15), xlim=c(0,50))+
    scale_x_continuous(breaks = c(0,16,24,36,48,60,72,84,96,144,180))+
    theme(text=element_text(size=25))

print(ser_surv1)

ggsave(paste0(results_figs_dir,"ser_surv1.png"),
       plot=ser_surv1, width=12, height=8, units='in')

ser_surv2<-ggplot() +
    geom_step(data = tvary_ipcwdf,
              aes(x = time, y = risk, linetype = strata, color = "IPCW"), linewidth=1) +
    geom_step(data = tvary_wdf1,
              aes(x = time, y = risk, linetype = strata, color = "1 doses"),linewidth=1)+
    geom_step(data = tvary_wdf2,
              aes(x = time, y = risk, linetype = strata, color = "2 doses"),linewidth=1)+
    geom_step(data = tvary_wdf3,
              aes(x = time, y = risk, linetype = strata, color = "3 doses")) +
    geom_step(data = tvary_ittdf,
              aes(x = time, y = risk, linetype = strata, color = "ITT"),linewidth=1)+
    scale_linetype_manual(name = "Strata", 
                          values = c("solid", "dotdash"), 
                          labels = c("ABC", "TDF")) +
    scale_color_manual(name = "Protocol", 
                       values = c("ITT"="#925E9FFF","IPCW" = "black", "1 doses" = "#DF8F44FF",
                                  "2 doses"= "darkblue" , "3 doses" = "#AD002AFF"),
                       breaks=c( "ITT","IPCW",  "3 doses" ,"2 doses", "1 doses")) +
    labs(color = "Protocol",
         y="Risk",
         x="Weeks",
         title="Risk of Composite Outcome")+
    theme_classic()+
    coord_cartesian(ylim=c(0,0.15), xlim=c(0,50))+
    scale_x_continuous(breaks = c(0,16,24,36,48,60,72,84,96,144,180))+
    theme(text=element_text(size=25))

print(ser_surv2)

ggsave(paste0(results_figs_dir,"ser_surv2.png"),
       plot=ser_surv2, width=12, height=8, units='in')

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## ONE DOSE DEVIATION
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## ONE TIME deviation weights----
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Make the weights 

df_dos1<-wt_model(data=df_tvary,response="target_cens_1",
                  exposure="abc_tdf",
                  covariates=c('new_race','pethnc', 'sex', 'strat1', 'cd4', 'logrna', 'lastcd4', 'lastlogrna','age','basecd4','baserna','Week' ),
                  spline_vars = c('cd4', 'logrna', 'lastcd4', 'lastlogrna','age','basecd4','baserna','Week' ),
                  time_var="Week")


## :::::::::::::::::::::::::::::::::::::::
### Weights Table ----
## :::::::::::::::::::::::::::::::::::::::

weights_display(df_dos1,response="target_cens_1",save=T, directory = results_tabs_dir,
                filename ='timevary_weights1dose_sum_DR098.tex' )


## :::::::::::::::::::::::::::::::::::::::
### Survival analysis ABC v TDF  results----
## :::::::::::::::::::::::::::::::::::::::


tvary_weighted1<-ggsurvfit::survfit2(Surv(time=Week, time2=end, event=outcome)~abc_tdf,
                                     data=df_dos1,
                                     weights=target_cens_1_wt,
                                     robust=T,
                                     id=patid)
tvary_week1_surv<-survminer::ggsurvplot_combine(list(tvary_weighted1,tvary_ipcw,tvary_itt),
                                                data=df_tvary%>% filter(outcome_week<207),
                                                fun="event",
                                                linetype=c('solid','dashed',
                                                           'solid','dashed',
                                                           'solid','dotted'),
                                                palette='jco',
                                                censor=FALSE,
                                                ylim=c(0,1),
                                                legend.title="Treatment Arms",
                                                legend.labs=c('Weighted ABC/3TC',
                                                              'Weighted TDF/FTC',
                                                              'IPCW ABC/3TC',
                                                              'IPCW TDF/FTC',
                                                              'ITT ABC/3TC',
                                                              'ITT TDF/FTC'),
                                                legend=c(0.8,0.8),
                                                title="1 Missed Dose, Time-varying Covariates")
print(week1_surv)

ggsave(paste0(results_figs_dir,"tvary_week1_surv.png"),
       plot=tvary_week1_surv$plot, width=12, height=5, units='in')

## :::::::::::::::::::::::::::::::::::::::::
### Risk Differences ABC v TDF  results----
## :::::::::::::::::::::::::::::::::::::::::

dev1_r<-risk_xtract(tvary_weighted1, name = "Dev1"); print(dev1_r)

## :::::::::::::::::::::::::::::::::::::::::
### Evaluating the range of proportion of missed doses----
## :::::::::::::::::::::::::::::::::::::::::

props<-df_tvary %>% mutate(
  dose_deviation_1=case_when(sum(target_cens_1, na.rm = T)>0 & cens_week_1<outcome_week~1,
                                               TRUE~0),
dose_deviation_2=case_when(sum(target_cens_2, na.rm = T)>0 & cens_week_2<outcome_week~1,
                           TRUE~0),
dose_deviation_3=case_when(sum(target_cens_3, na.rm = T)>0 & cens_week_3<outcome_week~1,
                           TRUE~0),
dose_deviation_5=case_when(sum(target_cens_5, na.rm = T)>0 & cens_week_5<outcome_week~1,
                           TRUE~0),
dose_deviation_10=case_when(sum(target_cens_10, na.rm = T)>0 & cens_week_10<outcome_week~1,
                            TRUE~0),
never_miss=if_else(dose_deviation_1==0,1,0)
 ) %>%  slice

prop_1<-df_tvary %>% filter(sum(target_cens_1, na.rm = T)>0 & cens_week_1<outcome_week) %>%
  filter(Week==cens_week_1) %>%
  mutate(tot_doses=if_else(rxcode=="ATV/rtv + ABC/3TC"|rxcode=="ATV/rtv + TDF/FTC",cens_week_1*7*4,cens_week_1*7*3),
                           tdose=tot_doses_miss, pdose=tdose/tot_doses, perc=pdose*100) %>%
  select(patid, tot_doses, tdose, pdose, perc)

prop_2<-df_tvary %>% filter(sum(target_cens_2, na.rm = T)>0 & cens_week_3<outcome_week) %>% filter(Week==cens_week_2) %>%
  mutate(tot_doses=if_else(rxcode=="ATV/rtv + ABC/3TC"|rxcode=="ATV/rtv + TDF/FTC",cens_week_2*7*4,cens_week_2*7*3),
         tdose=tot_doses_miss, pdose=tdose/tot_doses, perc=pdose*100) %>%
  select(patid, tot_doses, tdose, pdose, perc)

prop_3<-df_tvary %>% filter(sum(target_cens_3, na.rm = T)>0 & cens_week_3<outcome_week) %>% filter(Week==cens_week_3) %>%
  mutate(tot_doses=if_else(rxcode=="ATV/rtv + ABC/3TC"|rxcode=="ATV/rtv + TDF/FTC",cens_week_3*7*4,cens_week_3*7*3),
         tdose=tot_doses_miss, pdose=tdose/tot_doses, perc=pdose*100) %>%
  select(patid, tot_doses, tdose, pdose, perc)

prop_5<-df_tvary %>% filter(sum(target_cens_5, na.rm = T)>0 & cens_week_5<outcome_week) %>% filter(Week==cens_week_5) %>%
  mutate(tot_doses=if_else(rxcode=="ATV/rtv + ABC/3TC"|rxcode=="ATV/rtv + TDF/FTC",cens_week_5*7*4,cens_week_5*7*3),
         tdose=tot_doses_miss, pdose=tdose/tot_doses, perc=pdose*100) %>%
  select(patid, tot_doses, tdose, pdose, perc)

prop_10<-df_tvary %>% filter(sum(target_cens_10, na.rm = T)>0 & cens_week_10<outcome_week) %>% filter(Week==cens_week_10) %>%
  mutate(tot_doses=if_else(rxcode=="ATV/rtv + ABC/3TC"|rxcode=="ATV/rtv + TDF/FTC",cens_week_10*7*4,cens_week_10*7*3),
         tdose=tot_doses_miss, pdose=tdose/tot_doses, perc=pdose*100) %>%
  select(patid, tot_doses, tdose, pdose, perc)

p1<-prop_1 %>% pull(perc);p2<-prop_2 %>% pull(perc); p3<-prop_3 %>% pull(perc);p5<-prop_5 %>% pull(perc);p10<-prop_10 %>% pull(perc)
props<-list("1 Dose Missed"=p1,
            "2 Dose Missed"=p2,
            "3 Dose Missed"=p3,
            "5 Dose Missed"=p5,
            "10 Dose Missed"=p10)
props_stacked<-stack(props) %>% mutate(Protocol=ind)

prop_graph<-ggplot(props_stacked, aes(x = Protocol, y = values, fill=Protocol)) +
  geom_boxplot() +
  scale_fill_grey(start = 0.7, end = 0.7)+
  labs(x = "Protocol", y = "Percentage of Expected Doses Missed (%)") +
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x  = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12)); print(prop_graph)

ggsave(paste0(results_figs_dir,"prop_missed_doses.png"),
       plot=prop_graph, width=10, height=7, units='in')
