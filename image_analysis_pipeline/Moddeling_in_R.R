## creating the moddeling for our data

library(pacman)
pacman::p_load(tidyverse, lme4, ggplot2,here)
setwd("~/Documents/DCL/DCLmeetsRETUNE/image_analysis_pipeline")

#read data
df <- read_csv('~/Documents/DCL/DCLmeetsRETUNE/image_analysis_pipeline/counted_cells.csv')

########### run a lm #############

LM_model <- lm(cc_normalized_per_volume ~ group_id, df)
summary(LM_model)

# make it a function
run_lm <- function(df, IV, DP){
  LM_model <- lm(IV ~ DP, df)
  summary(LM_model)
}

#testing
run_lm(df, df$cc_normalized_per_volume, df$group_id)

########### run a lmer #############

LMER_model <- lmer(cc_normalized_per_volume ~ group_id + (1|subject_id), df)
summary(LMER_model)

# make it a function
run_lmer <- function(df, IV, DP, subject){
  LMER_model <- lm(IV ~ DP, df)
  summary(LMER_model)
}

#testing
run_lmer(df, df$cc_normalized_per_volume, df$group_id, df$subject_id)

########### run a glmer #############
df_pos <- filter(df, df$cc_normalized_per_volume > 0)
GLMER_model <- glmer(cc_normalized_per_volume ~ group_id + (1|subject_id), df_pos, family = inverse.gaussian(link = "1/mu^2"))
summary(GLMER_model)

# make it a function
run_glmer <- function(df, IV, DP, subject, distribution){
  LMER_model <- lm(IV ~ DP, df, distribution)
  summary(GLMER_model)
  }

#testing - doens't run as it doens't like the distribution argument - maybe not make it a function?
run_glmer(df, df$cc_normalized_per_volume, df$group_id, df$subject_id, inverse.gaussian(link = "1/mu^2"))
