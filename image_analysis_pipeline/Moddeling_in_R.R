## creating the moddeling for our data

library(pacman)
pacman::p_load(tidyverse, lme4, ggplot2,here, plyr)
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
run_glmer <- function(df, IV, DP, subject, family){
  LMER_model <- lm(IV ~ DP, df, family)
  summary(GLMER_model)
  }

#testing - doens't run as it doens't like the distribution argument - maybe not make it a function?
run_glmer(df, df$cc_normalized_per_volume, df$group_id, df$subject_id, family = inverse.gaussian(link = "1/mu^2"))

############# visualize the effects    #############
df_plot <- df
df_plot$group_id <- ifelse(df_plot$group_id == "wt", "Wildtype", df_plot$group_id)
df_plot$group_id <- ifelse(df_plot$group_id == "tg", "Transgenic", df_plot$group_id)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}  
df_plot$group_id <- as.factor(df_plot$group_id)

df_sum <- data_summary(df_plot, "cc_normalized_per_volume", "group_id")

plot_results <- function(data, IV, DV, color_list, xlab, ylab, text_size, font){
  ggplot(data=df_sum, aes(x=DV, y=IV, fill = DV)) +
    geom_bar(stat="identity", width = 0.8) +
    theme_minimal() + xlab(xlab) + ylab(ylab) + 
    geom_errorbar(aes(ymin=IV-sd, ymax=IV+sd), width=.2, position=position_dodge(.9)) + 
    scale_fill_manual(values=color_list) +
    theme(legend.position="none") +
    theme(text=element_text(size= text_size, family= font))
}
plot_results(df_sum, IV = df_sum$cc_normalized_per_volume, DV = df_sum$group_id, color_list = c('#1A85FF','#D41159'), xlab = "Group ID", ylab = "Normalized per volume", text_size = 12, font = "Arial")
