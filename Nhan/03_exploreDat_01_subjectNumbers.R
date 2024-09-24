rm(list = ls())

library(tidyverse)
library(openxlsx)
library(ggpubr)
library(webr)

# load data =======================================================================
load("processedDat/cohortDat.RData")

# Main code / data analysis -------------------------
CMV_cohorts_dat <- cohortDat$donor_info %>% 
  group_by(Cohort, CMV_IgG_Serology, SEX_BIRTH) %>% 
  summarise(n = n()) %>% drop_na() %>%
  mutate(CMV_IgG_Serology = ifelse(CMV_IgG_Serology == 1, "CMV+", "CMV-"),
         SEX_BIRTH = ifelse(SEX_BIRTH == 0, "male", "female"))

CMV_cohorts_dat %>% 
  ggplot(aes(x = Cohort, y = n, fill = CMV_IgG_Serology)) + 
  geom_bar(stat= "identity") + facet_wrap(~SEX_BIRTH) + 
  theme_bw() + theme(legend.position = "top")

cohortDat$donor_info %>%
  group_by(Cohort, SEX_BIRTH) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

cohortDat$donor_info %>%
  group_by(SEX_BIRTH) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

range(cohortDat$donor_info$AGE)
mean(cohortDat$donor_info$AGE)

range((cohortDat$donor_info %>% filter(Cohort == "Discovery"))$AGE)
mean((cohortDat$donor_info %>% filter(Cohort == "Discovery"))$AGE)

range((cohortDat$donor_info %>% filter(Cohort == "Validation"))$AGE)
mean((cohortDat$donor_info %>% filter(Cohort == "Validation"))$AGE)

# plot data with CMV status information -----------------------------------------------
plotDat <- cohortDat$donor_info %>% 
  mutate(CMV_status = ifelse(CMV_IgG_Serology == 1, "CMV+", 
                             ifelse(CMV_IgG_Serology == 0, "CMV-", NA)),
         SEX_BIRTH = ifelse(SEX_BIRTH == 0, "male", "female")) %>% 
           drop_na(CMV_status)

# number participant per CMV status
plotDat %>%
  group_by(Cohort, CMV_status) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# age range
range((plotDat %>% filter(Cohort == "Discovery", CMV_status == "CMV+"))$AGE)
range((plotDat %>% filter(Cohort == "Discovery", CMV_status == "CMV-"))$AGE)
range((plotDat %>% filter(Cohort == "Validation", CMV_status == "CMV+"))$AGE)
range((plotDat %>% filter(Cohort == "Validation", CMV_status == "CMV-"))$AGE)

#  number participant per CMV status with age median and mean
plotDat %>%
  group_by(Cohort, CMV_status) %>%
  summarise(age_median = median(AGE))

plotDat %>%
  group_by(Cohort, CMV_status) %>%
  summarise(age_mean = mean(AGE))

# number participant persex
plotDat %>%
  group_by(Cohort, SEX_BIRTH) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# donut chart ---------------------------------------------------------------------
PieDonut(plotDat, aes(CMV_status), 
         r0 = 0.5, start = -120,
         title = "Distribution of gender per season")

PieDonut(plotDat, aes(CMV_status), r0 = 0.5)

cowplot::plot_grid(
  plotDat %>% filter(Cohort == "Discovery") %>% 
    PieDonut(aes(CMV_status), r0 = 0.5) + labs(title = "Discovery"),
  plotDat %>% filter(Cohort == "Validation") %>% 
    PieDonut(aes(CMV_status), r0 = 0.5) + labs(title = "Validation"),
  nrow = 1
)


