######################################################################################
#SIR.fun calculate the Standardized Incidence Rate Ratios SIR of a desired event based on the
#historical background rates (BGR) and Incidence data of exposed individuals (Obs.DF).
#
#BGR should be a data.frame with 3 columns:
#     "AGE" which is a factor variable of the desired age groups ([0-4],[5-9], ..).
#     "SEX" which is a factor variable of sex ("M", "F").
#     "Rate" which is a numerical variable of the historical incidence per 100.000 PYRS
#
#Obs.DF should be a data.frame with at least 4 columns: 
#     "AGE" and "SEX" should be identical to the one in BGR
#     "PYRS" is a numerical variable of the follow-up time in the exposed individuals
#     "Obs" is a numerical variable of the recorded cases during follow-up in the exposed individuals
#     "..." additional covariates that could be of interest such as vaccination status or something similar
#
#The estimate of SIR can be stratified by defining "grp_vars" as a vector of the covariate names
#that we wish to stratify for. This could be "AGE" or "SEX" but also additional covariates, such as vaccination status,
#if it is included in Obs.DF 
############################################################################################

SIR.fun <- function(BGR,Obs.DF, grp_vars = NULL){
  
  #Take the data.frame that include the total observed PYRS of the exposed variable (COVID19 infection)
  #and combine it with the historical rates.
  SIR.DF <- Obs.DF  %>%
    left_join(BGR, by = c("AGE","SEX")) %>%
    group_by_at(.,vars(grp_vars)) %>%
    mutate(Expected = Rate*PYRS/100000) %>%
    dplyr::summarise(SIR = sum(Obs)/sum(Expected, na.rm = TRUE),
                     SIR_upr = sum(Obs)/(qchisq(0.025,2*sum(Expected, na.rm = TRUE))/2),
                     SIR_low = sum(Obs)/(qchisq(0.975,2*(sum(Expected, na.rm = TRUE)+1))/2),.groups = "keep")
  return(SIR.DF)
}

############################################################
# Example of utilization

#To illustrate the function we have simulated some background rates of an arbitrary event A
#and some observations from a period afterwards where a population have been exposed to something.
#In the example Females that have been exposed has 1.1 increased risk of having event A and 
#Males that have been exposed have 1.5 increased risk.


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Install necessary packages
source("Package.R")

#Load simulation data
BGR_sim <- read.csv("Data/Background_rates.csv")
Obs.DF <- read.csv("Data/ObservationalData.csv")


# Here we call the function with different arguments of grp_vars to illustrate
#the possibility of stratified SIR estimates.
SIR.fun(BGR_sim, Obs.DF, grp_vars = NULL)
SIR.fun(BGR_sim, Obs.DF, grp_vars = c("SEX"))
SIR.fun(BGR_sim, Obs.DF, grp_vars = c("AGE"))
SIR.fun(BGR_sim, Obs.DF, grp_vars = c("AGE", "SEX"))

############################################################
