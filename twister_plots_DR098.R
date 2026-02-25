
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
## NOTES ----
## :::::::::::::::::::::::::::::::::::::::
## Notes: Creating of twister plots to visualize the risk differences between ABC and TDF
## arms in each protocol and to visualize the differences within arms (either ABC of TDF)
## within each protocol
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
## PACKAGES ----
## :::::::::::::::::::::::::::::::::::::::
## list the packages we will need: 

packages<-c('ggplot2', 'dplyr', 'survival', 'ggsurvfit', 'here',
            'survminer', 'data.table', 'tidyr')      # list all the packages we need

## :::::::::::::::::::::::::::::::::::::::
## FUNCTIONS ----
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

###############################################################################
# Twister Plots
#
# Alex Breskin (2021/5/14), Paul Zivich (2021/6/28)
##############################################################################

library(ggplot2)
library(rlang)

#### Step Ribbon ####
StatStepribbon <- ggproto("StatStepribbon",
                          Stat,
                          compute_group=function(., data, scales, direction = "hv",
                                                 yvars = c( "ymin", "ymax" ), ...)
                          {
                            direction <- match.arg( direction, c( "hv", "vh" ) )
                            data <- as.data.frame( data )[ order( data$x ), ]
                            n <- nrow( data )
                            
                            if ( direction == "vh" ) {
                              xs <- rep( 1:n, each = 2 )[ -2 * n ]
                              ys <- c( 1, rep( 2:n, each = 2 ) )
                            } else {
                              ys <- rep( 1:n, each = 2 )[ -2 * n ]
                              xs <- c( 1, rep( 2:n, each = 2))
                            }
                            
                            data.frame(
                              x = data$x[ xs ]
                              , data[ ys, yvars, drop=FALSE ]
                              , data[ xs, setdiff( names( data ), c( "x", yvars ) ), drop=FALSE ]
                            )
                          },
                          required_aes=c( "x", "ymin", "ymax" ),
                          default_geom=GeomRibbon,
                          default_aes=aes( x=..x.., ymin = ..y.., ymax=Inf )
)

stat_stepribbon = function( mapping=NULL, data=NULL, geom="ribbon",
                            position="identity") {
  layer(stat=StatStepribbon, mapping=mapping, data=data, geom=geom, position=position )
}


#### Twister Plot####

#' Twister Plot
#'
#' @param dat A \code{data.frame} with the risk difference, upper and lower confidence limits, and times
#' @param xvar The variable name for the risk difference. Defaults to RD.
#' @param lcl  The variable name for the lower confidence limit of the risk difference. Defaults to RD_LCL.
#' @param ucl  The variable name for the upper confidence limit of the risk difference. Defaults to RD_UCL.
#' @param yvar The variable name for time. Defaults to "t".
#' @param xlab The x-axis label. Defaults to "Risk Difference".
#' @param ylab The y-axis label. Defaults to "Days".
#' @param treat_labs A vector containing the names of the treatment groups. Defaults to c("Treat", "Control")
#'
#' @return a \code{ggplot} object
#' @export
#'
#' @examples
#' twister_plot(dat, treat_labs = c("Vaccine", "Placebo"))
#' 
twister_plot <- function(dat,
                         xvar = RD,
                         lcl = RD_LCL,
                         ucl = RD_UCL,
                         yvar = t,
                         ylab = "Weeks",
                         xlab = "Risk difference",
                         reference_line = 0.0,
                         log_scale = FALSE,
                         treat_labs = c("TDF/FTC", "ABC/3TC")){
  
  base_breaks <- function(n = 10){
    function(x) {
      axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
    }
  }
  
  `%>%` <- magrittr::`%>%`
  pull <- dplyr::pull
  
  t_lim <- max(dat %>% pull({{yvar}}))
  if (log_scale) {
    x_lim <- max(abs(log(dat %>% pull({{lcl}}))), 
                 abs(log(dat %>% pull({{ucl}}))))
    y_scale = scale_y_continuous(limits = c(exp(-x_lim), exp(x_lim)), 
                                 trans="log", 
                                 breaks=base_breaks())
    text_loc = c(-x_lim/2, x_lim/2)
  } else {
    x_lim <- max(abs(dat %>% pull({{lcl}})), abs(dat %>% pull({{ucl}})))
    y_scale = scale_y_continuous(limits = c(-0.1, 0.02),
                                 breaks = c(-0.05, -0.02, 0, 0.02))
    text_loc = c(-x_lim/2, x_lim/2)
  }
  
  p <- ggplot(data = dat, aes(x = {{yvar}}, y = {{xvar}})) + 
    geom_line() +
    geom_ribbon(
      aes(ymin = {{lcl}}, ymax = {{ucl}}),
      alpha = 0.2,
    ) +
    geom_hline(yintercept = reference_line, linetype = "dotted") +
    y_scale + 
    scale_x_continuous(limits = c(0, t_lim), expand = c(0, 0)) +
    coord_flip(clip = "off") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          plot.margin = unit(c(2, 1, 1, 1), "lines")) +
    geom_text(data = head(dat, 1), 
              label = sprintf("Favors %s", treat_labs[1]), 
              x = t_lim+ 3, y = text_loc[1]) +
   #geom_text(data = head(dat, 1), 
   #          label = sprintf("Favors %s", treat_labs[2]), 
   #          x = t_lim+ 3, y = text_loc[2]) +
    xlab(ylab) +
    ylab(xlab)
  p
}



    ## :::::::::::::::::::::::::::::::::::::::


#### My Twister Augmentation
twister<-function(df, name, save=NULL,directory=NULL, plotname=NULL){
  temp<-df %>% filter(strata=="abc_tdf=ABC/3TC") %>%  mutate(ip_abc_r = 1-surv ) %>%
    fill(. , "ip_abc_r", .direction = "down") %>%  select(-surv) %>% arrange(time) %>% 
    left_join(df %>% filter(strata=="abc_tdf=TDF/FTC"), by="time") %>%  
    mutate(ip_tdf_r=1-surv) %>% fill(., c("ip_tdf_r","strata.y" ), .direction="down") %>%
    filter(if_all(everything(), ~ !is.na(.))) %>% 
    mutate(rd=ip_tdf_r-ip_abc_r,
           rdse=sqrt(std.error.x^2+std.error.y^2),
           rdmin=rd-rdse,
           rdmax=rd+rdse,
           dfname=name) %>%
    select(-c(surv,ip_abc_r,ip_tdf_r, strata.x,strata.y)) %>%
    filter(time<=96)
  
  assign(paste0(name,"_df"), temp, envir = .GlobalEnv)
  
  p<-twister_plot(temp, xvar=rd, lcl=rdmin, ucl=rdmax, yvar=time)+
    theme(element_text(size=20))
  assign(paste0(name,"_twister"), p, envir = .GlobalEnv)
  
  if (save) {
    if (is.null(plotname)) {
      stop("Error: 'plotname' required.")
    }
    if (is.null(directory)) {
      stop("Error: 'directory' required.")
    }
    
    # Create the directory if it doesn't exist
    if (!dir.exists(directory)) {
      stop("Error: 'directory' does not exist.")
    }
    
    # Save the data to the specified file
    filepath <- file.path(directory, plotname)
    ggsave(filepath, plot=p, width=12, height=8, units='in')

    message("Data saved to: ", filepath)
  }
  
}
#produces a  dataset called 'name_df' and a plot called 'name_twister'
## :::::::::::::::::::::::::::::::::::::::

## Establish Directories to use ----
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

tvary_ittdf<-read.csv(paste0(data_dir,"tvary_ittdf_DR098.csv")) 
tvary_ipcwdf<-read.csv(paste0(data_dir,"tvary_ipcwdf_DR098.csv")) 
tvary_1dosedf<-read.csv(paste0(data_dir,"tvary_1dose_DR098.csv"))
tvary_2dosedf<-read.csv(paste0(data_dir,"tvary_2dose_DR098.csv"))
tvary_3dosedf<-read.csv(paste0(data_dir,"tvary_3dose_DR098.csv"))
tvary_5dosedf<-read.csv(paste0(data_dir,"tvary_5dose_DR098.csv"))
tvary_10dosedf<-read.csv(paste0(data_dir,"tvary_10dose_DR098.csv"))
       
## :::::::::::::::::::::::::::::::::::::::

twister(df=tvary_ittdf, name="itt",save=T, directory=results_figs_dir, plotname="itt_twister.png")
#twister(df=tvary_ipcwdf, name="ipcw",save=T, directory=results_figs_dir, plotname="ipcw_twister.png")
twister(df=tvary_1dosedf, name="1dose",save=T, directory=results_figs_dir, plotname="1dose_twister.png")
twister(df=tvary_2dosedf, name="2dose",save=T, directory=results_figs_dir, plotname="2dose_twister.png")
twister(df=tvary_3dosedf, name="3dose",save=T, directory=results_figs_dir, plotname="3dose_twister.png")
twister(df=tvary_10dosedf, name="10dose",save=T, directory=results_figs_dir, plotname="10dose_twister.png")

#all plots together, one panel for each
multiple<-ggarrange(itt_twister,
                    #ipcw_twister,
                    `1dose_twister`,
          `2dose_twister`,`3dose_twister`, `10dose_twister`, ncol = 3, nrow=2,
          labels = c("A. ITT",
                     #"B. LTFU",
                     "B. 1 Dose Protocol",
                     "C. 2 Dose Protocol", "D. 3 Dose Protocol", "E. 10 Dose Protocol"),
          hjust = c(-1,
                    #-.50,
                    -.30,-.30,-.3,-.3))

filepath <- file.path(results_figs_dir, "multiple_twisters.png")
ggsave(filepath, plot=multiple, width=12, height=8, units='in')


 ## Lines combined on one plot
all_prot_data<-rbind(itt_df,
                     #ipcw_df,
                     `1dose_df`,`2dose_df`,`3dose_df`,`10dose_df`)
all_prot_data$dfname<-factor(all_prot_data$dfname,
                             levels = c("itt",
                                        #"ipcw",
                                        "10dose" ,"3dose","2dose", "1dose"),
                             labels=c("ITT",
                                      #"LTFU",
                                      "10 Dose Protocol","3 Dose Protocol","2 Dose Protocol","1 Dose Protocol"))

combined<-ggplot()+
  geom_line(data=all_prot_data, aes(x=time, y=rd, linetype=dfname), linewidth=0.9)+
  #geom_hline(yintercept = 0, linetype = "solid") +
  scale_y_continuous(limits = c(-0.065, 0.0001),
                     breaks = c(0,-0.09, -0.05, -0.04,-0.03,-0.02),
                     expand=c(0,0),
                     labels = function(x) {
                       ifelse(x == 0,        # for zero
                              "0",           # print it as “0”
                              sprintf("%.2f", # otherwise format to 2 decimals
                                      x))
                     }) + 
  scale_x_continuous(limits = c(0, 98), expand = c(0, 0),
                     breaks = c(16,24,36,48,60,72,84,96)) +
  coord_flip(clip = "off") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.margin = unit(c(2, 1, 1, 1), "lines"),
        axis.title = element_text(size = 18),       # Axis titles
        axis.text = element_text(size = 18),        # Axis tick labels
        legend.title = element_text(size = 18),     # Legend title
        legend.text = element_text(size = 18),     # Legend labels)+
        legend.position = c(.75,.8),
        legend.box.background = element_rect(),
        legend.box.margin = margin(.5,.5,.5,.5))+
  labs(linetype="Protocol")+
  ylab("Risk difference comparing TDF/FTC to ABC/3TC")+
  xlab("Time from randomization (Weeks)")+
  geom_text(data=head(all_prot_data,1),
            label = "← Favors TDF/FTC", 
            x = 100, y = -0.055/2,
            size=6);print(combined, width=10, height=8, units='in')

filepath <- file.path(results_figs_dir, "twisters_combined.png")
ggsave(filepath, plot=combined, width=10, height=8, units='in')


combined_alt<-ggplot()+
  geom_line(data=all_prot_data, aes(x=time, y=rd, linetype=dfname), linewidth=0.9)+
  geom_hline(yintercept = 0, linetype = "solid") +
  scale_y_continuous(limits = c(-0.065, 0.055),
                     breaks = c(0,-0.05, -0.04,-0.03,-0.02, 0.02, 0.03, 0.04, 0.05, 0.09),
                     expand=c(0,0),
                     labels = function(x) {
                       ifelse(x == 0,        # for zero
                              "0",           # print it as “0”
                              sprintf("%.2f", # otherwise format to 2 decimals
                                      x))
                     }) + 
  scale_x_continuous(limits = c(0, 98), expand = c(0, 0),
                     breaks = c(16,24,36,48,60,72,84,96)) +
  coord_flip(clip = "off") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.margin = unit(c(2, 1, 1, 1), "lines"),
        axis.title = element_text(size = 18),       # Axis titles
        axis.text = element_text(size = 18),        # Axis tick labels
        legend.title = element_text(size = 18),     # Legend title
        legend.text = element_text(size = 18) )+     # Legend labels)+
  labs(linetype="Protocol")+
  ylab("Risk difference comparing TDF/FTC to ABC/3TC")+
  xlab("Time from randomization (Weeks)")+
  geom_text(data=head(all_prot_data,1),
            label = "← Favors TDF/FTC", 
            x = 101, y = -0.055/2,
            size=6)+
geom_text(data=head(all_prot_data,1),
          label = "Favors ABC/3TC →", 
          x = 101, y = 0.055/2,
          size=6);combined_alt

filepath <- file.path(results_figs_dir, "twisters_combined_alt.png")
ggsave(filepath, plot=combined_alt, width=10, height=8, units='in')
   