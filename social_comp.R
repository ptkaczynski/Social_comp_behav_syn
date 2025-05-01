# Social plasticity analysis ####

# 1 Load data & packages ####

library(tidyverse)
library(igraph)
library(brms)
library(tidybayes)
library(reshape2)
library(ggdist)
library(modelr)
library(rstan)
library(bisonR)
library(assortnet)
library(bayestestR)
library(cmdstanr)
library(lme4)
library(viridis)
library(patchwork)
library(gridExtra)

save.image("results.Rdata")

# 2 Table S1: Summary of data ####

names(act_data)
sort(unique(act_data$Behaviour))

xx = filter(act_data, !is.na(duration))
xx = filter(xx, !Behaviour == "Out_of_sight")

# overall summary of for each group 

Group <- c(rep("Blue", 10), rep("Green",10))
Season <- c(
  "Mating",
  "Mating",
  "Mating",
  "Mating/Non-mating 1",
  "Non-mating 1",
  "Non-mating 1",
  "Non-mating 2",
  "Non-mating 2",
  "Non-mating 2",
  "Whole study",
  "Mating",
  "Mating",
  "Mating",
  "Mating",
  "Non-mating 1",
  "Non-mating 1",
  "Non-mating 2",
  "Non-mating 2",
  "Non-mating 2",
  "Whole study"
)

Time_window <- c(
  "1",
  "2",
  "3",
  "4*",
  "5",
  "6",
  "7",
  "8",
  "9",
  "-",
  "1",
  "2",
  "3",
  "4",
  "5",
  "6",
  "7",
  "8",
  "9",
  "-"
)
names(act_data)
weather <- act_data %>%
  group_by(visit_year) %>%
  summarise(
    yr_effort = length(visit_year)) 


yy = aggregate(xx$duration, by = list(xx$Subject, xx$Group),FUN=sum)
colnames(yy) = c("subject", "group", "total")
mean(yy$total[yy$group == "B"])/3600
mean(yy$total[yy$group == "G"])/3600
sd(yy$total[yy$group == "B"])/3600
sd(yy$total[yy$group == "G"])/3600

names(act_data)

mean(xx$fruit_avail[xx$Group == "B"])
mean(xx$fruit_avail[xx$Group == "G"])

xx$Tourist = as.numeric(as.character(xx$Tourist))
xx = filter(xx, !is.na(Tourist))
mean(xx$Tourist[xx$Group == "B"])
mean(xx$Tourist[xx$Group == "G"])

xx = filter(xx, !is.na(daily_min_temp))
mean(xx$daily_min_temp[xx$Group == "B"])
mean(xx$daily_min_temp[xx$Group == "G"])

# now per time window

xx = filter(act_data, !is.na(duration))
xx = filter(xx, !Behaviour == "Out_of_sight")

yy = aggregate(xx$duration, by = list(xx$Subject, xx$Group, xx$time_window),FUN=sum)
colnames(yy) = c("subject", "group", "time_window", "total")
zz = aggregate(yy$total, by = list(yy$group, yy$time_window),FUN=mean)
window_d$mean_obs = zz$x/3600
zz = aggregate(yy$total, by = list(yy$group, yy$time_window),FUN=sd)
window_d$sd_obs = zz$x/3600

tableS1 = window_d
names(tableS1)

tableS1$time_win_min_temp = round(tableS1$time_win_min_temp, 2)
tableS1$time_win_fruit = round(tableS1$time_win_fruit, 2)
tableS1$time_win_anthro = round(tableS1$time_win_anthro, 2)
tableS1$season = round(tableS1$season, 0)
tableS1$mean_obs = round(tableS1$mean_obs, 2)
tableS1$sd_obs = round(tableS1$sd_obs, 2)

colnames(tableS1) = c("Group", "Time window", "Mean minimum temp (C)", 
                      "Mean daily food availability", "Mean anthropogenic scan", "Season", 
                      "Mean obs time (hrs) per subject",
                      "SD obs time (hrs) per subject")

write.csv(tableS1, "C:/Users/Megaport/OneDrive/Project_data/social_comp/results/tables1.csv")

xx = filter(act_data, Group == "B")
yy = aggregate(xx$frequency, by = list(xx$Season, xx$time_window, xx$rrdate), FUN = sum)
colnames(yy) =  c("season", "time_window", "rrdate", "xx")
yy$freq = 1
zz = aggregate(yy$freq, by = list(yy$season, yy$time_window), FUN = sum)

rm(xx, yy, zz)

# 3 Packages for scan analysis ####

options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')

# 4 Check main variables ####

# Following Julia's comments, we will now check how the model results hold up using:
# 15-minute interval (main model reported in the manuscript)
# 5-minute interval (supplementary model 1)
# whole focal sample, i.e. 30 minute interval (supplementary model 2)

# first we only want where we have all our predictors
# second we want to check which of our predictors have reasonable variation

names(scan_d15)

# need to add tourists in proximity

names(focal_data)

focal_data$humans_in_prox = focal_data$tourists.in.1mS+focal_data$tourists.in.5mS+focal_data$tourists.in.10mS+
  focal_data$berbers.in.1mS+focal_data$berbers.in.5mS+focal_data$berbers.in.10mS
range(focal_data$humans_in_prox)

scan_d15$human_prox = NA

for (i in 1:nrow(scan_d15)){
  
  xx=which(focal_data$obs == scan_d15$obs[i])
  if(length(xx)>0){
    scan_d15$human_prox[i]=as.character(focal_data$humans_in_prox[xx[1]])}
}

scan_d5$human_prox = NA

for (i in 1:nrow(scan_d5)){
  
  xx=which(focal_data$obs == scan_d5$obs[i])
  if(length(xx)>0){
    scan_d5$human_prox[i]=as.character(focal_data$humans_in_prox[xx[1]])}
}

scan_d30$human_prox = NA

for (i in 1:nrow(scan_d30)){
  
  xx=which(focal_data$obs == scan_d30$obs[i])
  if(length(xx)>0){
    scan_d30$human_prox[i]=as.character(focal_data$humans_in_prox[xx[1]])}
}

##

pd15 = scan_d15
names(pd15)

pd15 = pd15 %>%
  mutate(subject_ID = as.character(Subject),
         time = as.numeric(as.character(time.n)),
         social_choice = as.numeric(as.character(groom_scan)),
         obs_ID = as.character(obs),
         group = as.character(Group),
         sex = as.character(Sex),
         season = as.numeric(as.character(Season)),
         anthro = as.numeric(as.character(Tourist_scan)),
         anthro_prox = as.numeric(as.character(human_prox)),
         temp = as.numeric(as.character(daily_min_temp)),
         food = as.numeric(as.character(fruit_avail)),
         own_rank = as.numeric(as.character(rank)),
         date_obs = as.character(rdate),
         conspecific_prox = as.numeric(as.character(adults.in.prox)),
         rank_diff = as.numeric(as.character(max_rank_diff)),
         previous_agg = as.numeric(as.character(previous_aggs)),
         behavioural_syndrome = as.numeric(as.character(excitable)))

pd5 = scan_d5
names(pd5)

pd5 = pd5 %>%
  mutate(subject_ID = as.character(Subject),
         time = as.numeric(as.character(time.n)),
         social_choice = as.numeric(as.character(groom_scan)),
         obs_ID = as.character(obs),
         group = as.character(Group),
         sex = as.character(Sex),
         season = as.numeric(as.character(Season)),
         anthro = as.numeric(as.character(Tourist_scan)),
         anthro_prox = as.numeric(as.character(human_prox)),
         temp = as.numeric(as.character(daily_min_temp)),
         food = as.numeric(as.character(fruit_avail)),
         own_rank = as.numeric(as.character(rank)),
         date_obs = as.character(rdate),
         conspecific_prox = as.numeric(as.character(adults.in.prox)),
         rank_diff = as.numeric(as.character(max_rank_diff)),
         previous_agg = as.numeric(as.character(previous_aggs)),
         behavioural_syndrome = as.numeric(as.character(excitable)))

pd30 = scan_d30
names(pd30)

pd30 = pd30 %>%
  mutate(subject_ID = as.character(Subject),
         time = as.numeric(as.character(time.n)),
         social_choice = as.numeric(as.character(groom_scan)),
         obs_ID = as.character(obs),
         group = as.character(Group),
         sex = as.character(Sex),
         season = as.numeric(as.character(Season)),
         anthro = as.numeric(as.character(Tourist_scan)),
         anthro_prox = as.numeric(as.character(human_prox)),
         temp = as.numeric(as.character(daily_min_temp)),
         food = as.numeric(as.character(fruit_avail)),
         own_rank = as.numeric(as.character(rank)),
         date_obs = as.character(rdate),
         conspecific_prox = as.numeric(as.character(adults.in.prox)),
         rank_diff = as.numeric(as.character(max_rank_diff)),
         previous_agg = as.numeric(as.character(previous_aggs)),
         behavioural_syndrome = as.numeric(as.character(excitable)))

names(pd15)

pd15 = pd15[,c(24:40)]

names(pd5)

pd5 = pd5[,c(24:40)]

names(pd30)

pd30 = pd30[,c(24:40)]

##

str(pd15)

pd15 = filter(pd15, !is.na(temp))
pd15 = subset(pd15, !is.na(food))
range(pd15$time) # fine

str(pd5)

pd5 = filter(pd5, !is.na(temp))
pd5 = subset(pd5, !is.na(food))
range(pd5$time) # fine

str(pd30)

pd30 = filter(pd30, !is.na(temp))
pd30 = subset(pd30, !is.na(food))
range(pd30$time) # fine

# now histogram of predictors

str(pd15)

hist(pd15$time) #fine
hist(pd15$social_choice)
length(pd15$social_choice[pd15$social_choice == 1])#164; 11% of scans
hist(pd15$anthro)
hist(pd15$anthro_prox)
hist(pd15$temp)
hist(pd15$food)
hist(pd15$conspecific_prox)
hist(pd15$rank_diff)
hist(pd15$previous_agg)
length(pd15$previous_agg[pd15$previous_agg == 1 & pd15$social_choice == 1]) #2
hist(pd15$behavioural_syndrome)

hist(pd5$time) #fine
hist(pd5$social_choice)
length(pd5$social_choice[pd5$social_choice == 1])#66; 4% of scans

hist(pd30$time) #fine
hist(pd30$social_choice)
length(pd30$social_choice[pd30$social_choice == 1])#248; 17% of scans

# no variation for the previous agg variable so let's ditch it

pd15$previous_agg = NULL
pd5$previous_agg = NULL
pd30$previous_agg = NULL

# 5 Create time proportion variable ####

## Time needs to be a vector of numbers that specifies times as the
## proportion of a 24-hour day that has elapsed at the observed time. For
## example, noon yields a value of 0.5 whereas an observation as 6:00 PM
## yields a value of 0.75. Current our value is seconds in a day, the total of 
## which is 24*3600, 86400

pd15$time_prop = pd15$time/86400
range(pd15$time_prop) # looks ok to me.

pd5$time_prop = pd5$time/86400
range(pd5$time_prop) # looks ok to me.

pd30$time_prop = pd30$time/86400
range(pd30$time_prop) # looks ok to me.

# 6 Summary stats of scans ####

str(pd15)
length(levels(as.factor(pd15$obs_ID))) #1468

sort(unique(act_data$Behaviour))
check = filter(act_data, !Behaviour == "Out_of_sight")
names(check)
sum(check$duration/3600)
length(unique(check$observation))
check_tab = aggregate(check$duration, by=list(check$Subject), FUN=sum)
sd(check_tab$x/3600)

# 7 Standardize variables ####

names(pd15)

pd15$syndrome_z <- (pd15$behavioural_syndrome - mean(pd15$behavioural_syndrome))/sd(pd15$behavioural_syndrome)
range(pd15$syndrome_z)

range(pd15$season) # currently three factors, just want to distinguish between mating & non-mating
pd15$seasonx = ifelse(pd15$season == 1, 1, 0)
range(pd15$seasonx)

pd15$anthrop_z <- (pd15$anthro_prox - mean(pd15$anthro_prox))/sd(pd15$anthro_prox)
range(pd15$anthrop_z)

pd15$food_z <- (pd15$food - mean(pd15$food))/sd(pd15$food)
range(pd15$food_z)

pd15$temp_z <- (pd15$temp - mean(pd15$temp))/sd(pd15$temp)
range(pd15$temp_z)

pd15$consprox_z <- (pd15$conspecific_prox - mean(pd15$conspecific_prox))/sd(pd15$conspecific_prox)
range(pd15$consprox_z)

pd15$rankdiff_z <- (pd15$rank_diff - mean(pd15$rank_diff))/sd(pd15$rank_diff)
range(pd15$rankdiff_z)

pd15$rank_z <- (pd15$own_rank - mean(pd15$own_rank))/sd(pd15$own_rank)
range(pd15$rank_z)

pd15$time_z <- (pd15$time_prop - mean(pd15$time_prop))/sd(pd15$time_prop)
range(pd15$time_z)

# have just realised that I don't need four different sheets

pd = pd15
names(pd)

pd$social_choice15 = as.numeric(as.character(pd$social_choice))
pd$social_choice10 = as.numeric(as.character(pd10$social_choice))
pd$social_choice30 = as.numeric(as.character(pd30$social_choice))
pd$social_choice5 = as.numeric(as.character(pd5$social_choice))

str(pd)

# 8 Create pairs plots to check predictor variables ####

names(pd)

covees = pd[,c(5,6,18:26)]
colnames(covees) = c("group", "sex", "seasonx","syndrome_z","anthrop_z",           
                     "food_z","temp_z","consprox_z","rankdiff_z","rank_z",              
                     "time_z")

library(car)

mod = lm(social_choice ~ seasonx + syndrome_z + anthrop_z + food_z + temp_z + consprox_z +  
           rankdiff_z + rank_z + time_z + sex + group, data = pd)

vif(mod)

'   seasonx syndrome_z  anthrop_z     food_z     temp_z consprox_z rankdiff_z     rank_z     time_z        sex      group 
  1.257196   2.838315   1.039386   1.345249   1.050205   1.123592   1.531976   3.238556   1.039255   4.411720   1.570015 '

rm(mod)

library(ggplot2)
library(GGally)
windows()
ggpairs(covees) # all <=0.6

rm(covees)

# 9 Check for model complexity ####

# following Julia comments, only really makes sense to include the immediate social environment as
# predictors; and also we only want the models with all the interactions - the single effects are
# not what we're really interested in.

names(pd)

modc = glmer(social_choice ~ syndrome_z*anthrop_z + syndrome_z*consprox_z +  
               syndrome_z*rankdiff_z + temp_z + rank_z + time_z +I(time_z^2) + sex + group +
               (anthrop_z+temp_z+consprox_z+rankdiff_z|subject_ID)+
               (1|date_obs), data = pd, family="binomial")

round(length(residuals(modc))/(length(fixef(modc))+length(summary(modc)$varcor)+1),1)
#91.8 so no issue with complexity, but it did fail to converge

# 10 Fit our model ####

names(pd)

mprior = get_prior(social_choice15 ~ 
                     syndrome_z*anthrop_z + 
                     syndrome_z*consprox_z +  
                     syndrome_z*rankdiff_z + 
                     temp_z + rank_z + time_z +I(time_z^2) + sex + group +
                     (syndrome_z + anthrop_z + consprox_z + rankdiff_z + temp_z + rank_z + time_z +I(time_z^2)|subject_ID)+
                     (1|date_obs), data = pd, 
                   family = bernoulli(link="logit"))

mprior$prior[2:14] <- "normal(0,1)"

make_stancode(social_choice15 ~ 
                syndrome_z*anthrop_z + 
                syndrome_z*consprox_z +  
                syndrome_z*rankdiff_z + 
                temp_z + rank_z + time_z +I(time_z^2) + sex + group +
                (syndrome_z + anthrop_z + consprox_z + rankdiff_z + temp_z + rank_z + time_z +I(time_z^2)|subject_ID)+
                (1|date_obs), data = pd, 
              family = bernoulli(link="logit"), prior = mprior)

mod15 = brm(backend = 'cmdstanr', formula = social_choice15 ~ 
              syndrome_z*anthrop_z + 
              syndrome_z*consprox_z +  
              syndrome_z*rankdiff_z + 
              temp_z + rank_z + time_z +I(time_z^2) + sex + group +
              (syndrome_z + anthrop_z + consprox_z + rankdiff_z + temp_z + rank_z + time_z +I(time_z^2)|subject_ID)+
              (1|date_obs), data = pd, 
            family = bernoulli(link="logit"), prior = mprior, 
            chains = 3, cores = 3, iter=5000, warmup = 2000, thin = 2, 
            control = list(adapt_delta = 0.95), sample_prior = "yes")

summary(mod15)
plot(mod15)
round(posterior_summary(mod15, probs = c(0.05, 0.95)),3)[c(1:27),]

saveRDS(mod15, file = "C:/Users/Megaport/OneDrive/Project_data/social_comp/results/mod15.RDS")

# 11 Other interval models ####

# 5 minute

mprior2 = get_prior(social_choice5 ~ 
                      syndrome_z*anthrop_z + 
                      syndrome_z*consprox_z +  
                      syndrome_z*rankdiff_z + 
                      temp_z + rank_z + time_z +I(time_z^2) + sex + group +
                      (syndrome_z + anthrop_z + consprox_z + rankdiff_z + temp_z + rank_z + time_z +I(time_z^2)|subject_ID)+
                      (1|date_obs), data = pd, 
                    family = bernoulli(link="logit"))

mprior2$prior[2:14] <- "normal(0,1)"

make_stancode(social_choice5 ~ 
                syndrome_z*anthrop_z + 
                syndrome_z*consprox_z +  
                syndrome_z*rankdiff_z + 
                temp_z + rank_z + time_z +I(time_z^2) + sex + group +
                (syndrome_z + anthrop_z + consprox_z + rankdiff_z + temp_z + rank_z + time_z +I(time_z^2)|subject_ID)+
                (1|date_obs), data = pd, 
              family = bernoulli(link="logit"), prior = mprior4)

mod5 = brm(backend = 'cmdstanr', 
           formula = social_choice5 ~ syndrome_z*anthrop_z + 
             syndrome_z*consprox_z +  
             syndrome_z*rankdiff_z + 
             temp_z + rank_z + time_z +I(time_z^2) + sex + group +
             (syndrome_z + anthrop_z + consprox_z + rankdiff_z + temp_z + rank_z + time_z +I(time_z^2)|subject_ID)+
             (1|date_obs), data = pd, 
           family = bernoulli(link="logit"), prior = mprior3, 
           chains = 3, cores = 3, iter=5000, warmup = 2000, thin = 2, 
           control = list(adapt_delta = 0.95), sample_prior = "yes")

round(posterior_summary(mod5, probs = c(0.05, 0.95)),3)

saveRDS(mod5, file = "C:/Users/Megaport/OneDrive/Project_data/social_comp/results/mod5.RDS")

# 30 minute

mprior3 = get_prior(social_choice30 ~  
                      syndrome_z*anthrop_z + 
                      syndrome_z*consprox_z +  
                      syndrome_z*rankdiff_z + 
                      temp_z + rank_z + time_z +I(time_z^2) + sex + group +
                      (syndrome_z + anthrop_z + consprox_z + rankdiff_z + temp_z + rank_z + time_z +I(time_z^2)|subject_ID)+
                      (1|date_obs), data = pd, 
                    family = bernoulli(link="logit"))

mprior3$prior[2:14] <- "normal(0,1)"

make_stancode(social_choice30 ~  
                syndrome_z*anthrop_z + 
                syndrome_z*consprox_z +  
                syndrome_z*rankdiff_z + 
                temp_z + rank_z + time_z +I(time_z^2) + sex + group +
                (syndrome_z + anthrop_z + consprox_z + rankdiff_z + temp_z + rank_z + time_z +I(time_z^2)|subject_ID)+
                (1|date_obs), data = pd, 
              family = bernoulli(link="logit"), prior = mprior3)

mod30 = brm(backend = 'cmdstanr', 
            formula = social_choice30 ~  
              syndrome_z*anthrop_z + 
              syndrome_z*consprox_z +  
              syndrome_z*rankdiff_z + 
              temp_z + rank_z + time_z +I(time_z^2) + sex + group +
              (syndrome_z + anthrop_z + consprox_z + rankdiff_z + temp_z + rank_z + time_z +I(time_z^2)|subject_ID)+
              (1|date_obs), data = pd, 
            family = bernoulli(link="logit"), prior = mprior3, 
            chains = 3, cores = 3, iter=5000, warmup = 2000, thin = 2, 
            control = list(adapt_delta = 0.95), sample_prior = "yes")

round(posterior_summary(mod30, probs = c(0.05, 0.95)),3)

saveRDS(mod30, file = "C:/Users/Megaport/OneDrive/Project_data/social_comp/results/red_mod30.RDS")

save.image("C:/Users/Megaport/Desktop/working.Rdata")

# 12 Model 1 Estimates ####

round(posterior_summary(mod15, probs = c(0.05, 0.95)),3)[c(1:27),]

m1.postb <- posterior_samples(mod15)
names(m1.postb)

sum(m1.postb$b_syndrome_z > 0) / length(m1.postb$b_syndrome_z)
sum(m1.postb$b_anthrop_z < 0) / length(m1.postb$b_anthrop_z)
sum(m1.postb$b_consprox_z > 0) / length(m1.postb$b_consprox_z)
sum(m1.postb$b_rankdiff_z > 0) / length(m1.postb$b_rankdiff_z)
sum(m1.postb$b_temp_z > 0) / length(m1.postb$b_temp_z)
sum(m1.postb$b_rank_z < 0) / length(m1.postb$b_rank_z)
sum(m1.postb$b_time_z < 0) / length(m1.postb$b_time_z)
sum(m1.postb$b_Itime_zE2 < 0) / length(m1.postb$b_Itime_zE2)
sum(m1.postb$b_sexM < 0) / length(m1.postb$b_sexM)
sum(m1.postb$b_groupG > 0) / length(m1.postb$b_groupG)

sum(m1.postb$'b_syndrome_z:anthrop_z' > 0) / length(m1.postb$'b_syndrome_z:anthrop_z')
sum(m1.postb$'b_syndrome_z:consprox_z' < 0) / length(m1.postb$'b_syndrome_z:consprox_z')
sum(m1.postb$'b_syndrome_z:rankdiff_z' > 0) / length(m1.postb$'b_syndrome_z:rankdiff_z')

# 13 Supplementary Model Estimates ####

# 5-minute model

round(posterior_summary(mod5, probs = c(0.05, 0.95)),3)[c(1:27),]

m1.postb <- posterior_samples(mod5)
names(m1.postb)

sum(m1.postb$b_syndrome_z > 0) / length(m1.postb$b_syndrome_z)
sum(m1.postb$b_anthrop_z < 0) / length(m1.postb$b_anthrop_z)
sum(m1.postb$b_consprox_z > 0) / length(m1.postb$b_consprox_z)
sum(m1.postb$b_rankdiff_z > 0) / length(m1.postb$b_rankdiff_z)
sum(m1.postb$b_temp_z > 0) / length(m1.postb$b_temp_z)
sum(m1.postb$b_rank_z < 0) / length(m1.postb$b_rank_z)
sum(m1.postb$b_time_z < 0) / length(m1.postb$b_time_z)
sum(m1.postb$b_Itime_zE2 < 0) / length(m1.postb$b_Itime_zE2)
sum(m1.postb$b_sexM < 0) / length(m1.postb$b_sexM)
sum(m1.postb$b_groupG > 0) / length(m1.postb$b_groupG)

sum(m1.postb$'b_syndrome_z:anthrop_z' < 0) / length(m1.postb$'b_syndrome_z:anthrop_z')
sum(m1.postb$'b_syndrome_z:consprox_z' < 0) / length(m1.postb$'b_syndrome_z:consprox_z')
sum(m1.postb$'b_syndrome_z:rankdiff_z' > 0) / length(m1.postb$'b_syndrome_z:rankdiff_z')

# 30-minute model

round(posterior_summary(mod30, probs = c(0.05, 0.95)),3)[c(1:27),]

m1.postb <- posterior_samples(mod30)
names(m1.postb)

sum(m1.postb$b_syndrome_z > 0) / length(m1.postb$b_syndrome_z)
sum(m1.postb$b_anthrop_z < 0) / length(m1.postb$b_anthrop_z)
sum(m1.postb$b_consprox_z > 0) / length(m1.postb$b_consprox_z)
sum(m1.postb$b_rankdiff_z > 0) / length(m1.postb$b_rankdiff_z)
sum(m1.postb$b_temp_z > 0) / length(m1.postb$b_temp_z)
sum(m1.postb$b_rank_z < 0) / length(m1.postb$b_rank_z)
sum(m1.postb$b_time_z < 0) / length(m1.postb$b_time_z)
sum(m1.postb$b_Itime_zE2 < 0) / length(m1.postb$b_Itime_zE2)
sum(m1.postb$b_sexM < 0) / length(m1.postb$b_sexM)
sum(m1.postb$b_groupG < 0) / length(m1.postb$b_groupG)

sum(m1.postb$'b_syndrome_z:anthrop_z' > 0) / length(m1.postb$'b_syndrome_z:anthrop_z')
sum(m1.postb$'b_syndrome_z:consprox_z' < 0) / length(m1.postb$'b_syndrome_z:consprox_z')
sum(m1.postb$'b_syndrome_z:rankdiff_z' > 0) / length(m1.postb$'b_syndrome_z:rankdiff_z')

# 14 PPC of social choice model ####

pp_check(mod15)
conditional_effects(mod15) # splits by one sd

length(pd$syndrome_z[pd$syndrome_z <= -1]) #142
length(pd$syndrome_z[pd$syndrome_z > -1 & pd$syndrome_z <1]) #1065
length(pd$syndrome_z[pd$syndrome_z >=1]) #261

# so we need a ppc for relationship between:
# time of day
# sex
# excitable*number of conspecifics

pred_vals = predict(mod15, summary = FALSE)
predictions = as.data.frame(t(pred_vals))

# function to extract mode from each row (via liz)
mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

memory.limit()
memory.limit(size=56000)

smode = mode(predictions)

fit_mod15 <- 
  fitted(mod15) %>%
  as_tibble() %>%
  bind_cols(pd)

fit_mod15$predicted = smode$V1

# create new plotting object with 100 draws from posterior

# now proper plot by randomly selecting 100 draws from the posterior

xx = floor(runif(100, 1, 3000))

# extract those from our predictions object

zz = predictions[,xx]

colnames(zz) = c("pred1", "pred2", "pred3", "pred4", "pred5", "pred6", "pred7", "pred8", "pred9", "pred10",
                 "pred11", "pred12", "pred13", "pred14", "pred15", "pred16", "pred17", "pred18", "pred19", "pred20",
                 "pred21", "pred22", "pred23", "pred24", "pred25", "pred26", "pred27", "pred28", "pred29", "pred30",
                 "pred31", "pred32", "pred33", "pred34", "pred35", "pred36", "pred37", "pred38", "pred39", "pred40",
                 "pred41", "pred42", "pred43", "pred44", "pred45", "pred46", "pred47", "pred48", "pred49", "pred50",
                 "pred51", "pred52", "pred53", "pred54", "pred55", "pred56", "pred57", "pred58", "pred59", "pred60",
                 "pred61", "pred62", "pred63", "pred64", "pred65", "pred66", "pred67", "pred68", "pred69", "pred70",
                 "pred71", "pred72", "pred73", "pred74", "pred75", "pred76", "pred77", "pred78", "pred79", "pred80",
                 "pred81", "pred82", "pred83", "pred84", "pred85", "pred86", "pred87", "pred88", "pred89", "pred90",
                 "pred91", "pred92", "pred93", "pred94", "pred95", "pred96", "pred97", "pred98", "pred99", "pred100")

# create new plotting object

names(fit_mod15)
plot_data = fit_mod15[, c(1,3,4,21,30,10,20,22,18,27,31,35)]
names(plot_data)
colnames(plot_data) = c("Estimate", "lower", "upper", "time", "time_z", "sex",
                        "behavioural_syndrome", "syndrome_z", "cons_prox", "cons_prox_z",
                        "observed","predicted")

plot_data = cbind(plot_data, zz)

time_ppc = ggplot(plot_data)+
  geom_jitter(aes(x=time, y=observed+0.025),height=0.05, alpha=1/10)+
  geom_jitter(aes(x=time, y=predicted), height=0.05, alpha=1/10, colour = "#E16462FF")+
  geom_smooth(aes(x=time, y=observed), method = "glm", colour = 'black', 
              formula = y ~ poly(x,2),
              method.args = list(family = "binomial"), 
              se = TRUE)+
  geom_smooth(aes(x=time, y=predicted), method = "glm", colour = '#E16462FF', 
              formula = y ~ poly(x,2),
              method.args = list(family = "binomial"), 
              se = TRUE)+
  labs(x = "Time of day (as proportion)",
       y = "Probability of initiating grooming")
time_ppc = time_ppc+theme_classic()
time_ppc = time_ppc+theme(legend.position="none")
time_ppc = time_ppc+ggtitle("(a)")
time_ppc

mean_male = mean(plot_data$Estimate[plot_data$sex == "M"])
sd_male = sd(plot_data$Estimate[plot_data$sex == "M"])
lower_male = mean(plot_data$lower[plot_data$sex == "M"])
upper_male = mean(plot_data$upper[plot_data$sex == "M"])
mean_fem = mean(plot_data$Estimate[plot_data$sex == "F"])
sd_fem = sd(plot_data$Estimate[plot_data$sex == "F"])
lower_fem = mean(plot_data$lower[plot_data$sex == "F"])
upper_fem = mean(plot_data$upper[plot_data$sex == "F"])

sex = c("M", "F")
means2 = c(mean_male, mean_fem)
lowers2 = c(lower_male, lower_fem)
uppers2 = c(upper_male, upper_fem)
sds2 = c(sd_male,sd_fem)

plot_data2 = cbind(means2,lowers2,uppers2,sds2)
plot_data2 = as.data.frame(plot_data2)
plot_data2$sex = sex

sex_ppc = ggplot(plot_data)+
  geom_jitter(aes(x=sex, y=observed+0.025),height=0.05, alpha=1/10)+
  geom_jitter(aes(x=sex, y=predicted), height=0.05, alpha=1/10, colour = "#E16462FF")+
  geom_errorbar(data=plot_data2, aes(x=sex,ymin=means2-sds2, ymax=means2+sds2), width=0.5)+
  geom_point(data=plot_data2, aes(x=sex,y=means2))+
  labs(x = "Sex", y = "Probability of initiating grooming")
sex_ppc = sex_ppc+theme_classic()
#season_ppc = season_ppc+theme(legend.position="none")
sex_ppc = sex_ppc+ggtitle("(b)")
sex_ppc


# check with just two levels how it looks
names(plot_data)
round(mean(plot_data$syndrome_z),3)
length(plot_data$syndrome_z[plot_data$syndrome_z <=0])
length(plot_data$syndrome_z[plot_data$syndrome_z >0])

bold1 = filter(plot_data, syndrome_z <=0)
bold1$behavioural_syndrome2 = "Low"
bold2 = filter(plot_data, syndrome_z >0)
bold2$behavioural_syndrome2 = "High"
test=rbind(bold1,bold2)
plot_data = test
rm(bold1, bold2, test)
plot_data = plot_data %>%
  mutate(behavioural_syndrome2 = gsub('Low', 'Shy', behavioural_syndrome2),
         behavioural_syndrome2 = gsub('High', 'Bold', behavioural_syndrome2))

bs_ppc = ggplot(plot_data)+
  geom_jitter(aes(x=cons_prox, y=observed+0.025),height=0.05, alpha=1/10)+
  geom_jitter(aes(x=cons_prox, y=predicted), height=0.05, alpha=1/10, colour = "#E16462FF")+
  geom_smooth(aes(x=cons_prox, y=observed), method = "glm", colour = "black",
              method.args = list(family = "binomial"), 
              se = TRUE)+
  geom_smooth(aes(x=cons_prox, y=predicted), method = "glm", colour = "#E16462FF", 
              method.args = list(family = "binomial"), 
              se = TRUE)+
  facet_wrap(~behavioural_syndrome2)+labs(x = "Conspecifics within 10m",
                                          y = "Probability of initiating grooming")
bs_ppc = bs_ppc+theme_classic()
bs_ppc = bs_ppc+theme(legend.position="none")
bs_ppc = bs_ppc+ggtitle("(c)")
bs_ppc

mod1_ppc = time_ppc+sex_ppc+bs_ppc

rm(bs_ppc)
rm(plot_data2, sex_ppc, time_ppc, mod1_ppc)
rm(lower_fem, lower_male, lowers2)
rm(mean_fem, mean_male, means2)

# 13 Figure 1 - main mod1 result ####

# so now we want to use the draws to plot between bold and shy

# probably best to split the data, make two plots and put them together

shy = filter(plot_data, behavioural_syndrome2 == "Shy")
bold = filter(plot_data, behavioural_syndrome2 == "Bold")

names(bold)
fig1a = ggplot(bold)+
  geom_line(colour = "black", aes(x=cons_prox, y=observed),
            stat = "smooth", method = "glm", method.args = list(family = 'binomial'),
            se = FALSE, size=1, alpha=1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred2),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred3),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred4),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred5),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred6),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred7),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred8),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred9),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred10),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred11),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred12),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred13),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred14),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred15),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred16),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred17),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred18),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred19),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred20),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred21),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred22),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred23),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred24),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred25),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred26),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred27),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred28),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred29),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred30),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred31),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred32),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred33),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred34),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred35),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred36),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred37),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred38),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred39),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred40),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred41),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred42),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred43),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred44),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred45),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred46),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred47),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred48),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred49),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred50),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred51),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred52),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred53),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred54),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred55),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred56),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred57),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred58),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred59),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred60),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred61),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred63),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred64),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred65),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred66),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred67),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred68),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred69),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred70),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred71),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred72),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred73),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred74),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred75),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred76),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred77),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred78),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred79),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred80),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred81),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred82),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred83),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred84),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred85),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred86),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred87),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred88),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred89),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred90),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred91),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred92),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred93),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred94),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred95),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred96),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred97),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred98),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred99),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred100),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)
fig1a = fig1a + coord_cartesian(xlim = c(0, 5), ylim = c(0, 0.55))
fig1a = fig1a + labs(x = "Conspecifics within 10m",
                     y = "Probability of initiating grooming")
fig1a = fig1a+theme_classic()
fig1a = fig1a+ggtitle("Bold")   
fig1a

fig1b = ggplot(shy)+
  geom_line(colour = "black", aes(x=cons_prox, y=observed),
            stat = "smooth", method = "glm", method.args = list(family = 'binomial'),
            se = FALSE, size=1, alpha=1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred2),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred3),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred4),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred5),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred6),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred7),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred8),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred9),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred10),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred11),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred12),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred13),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred14),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred15),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred16),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred17),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred18),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred19),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred20),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred21),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred22),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred23),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred24),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred25),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred26),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred27),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred28),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred29),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred30),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred31),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred32),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred33),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred34),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred35),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred36),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred37),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred38),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred39),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred40),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred41),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred42),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred43),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred44),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred45),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred46),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred47),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred48),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred49),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred50),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred51),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred52),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred53),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred54),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred55),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred56),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred57),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred58),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred59),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred60),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred61),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred63),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred64),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred65),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred66),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred67),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred68),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred69),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred70),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred71),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred72),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred73),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred74),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred75),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred76),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred77),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred78),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred79),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred80),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred81),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred82),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred83),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred84),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred85),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred86),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred87),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred88),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred89),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred90),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred91),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred92),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred93),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred94),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred95),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred96),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred97),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred98),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred99),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#E16462FF", aes(x=cons_prox, y=pred100),
            stat = "smooth", method = "glm", method.args = list(family = "binomial"),
            se = FALSE, size=1, alpha=0.1)
fig1b = fig1b + coord_cartesian(xlim = c(0, 5), ylim = c(0, 0.55))
fig1b = fig1b + labs(x = "Conspecifics within 10m",
                     y = "")
fig1b = fig1b+theme_classic()
fig1b = fig1b+ggtitle("Shy")   
fig1b

fig1 = fig1a+fig1b
fig1

rm(fig1, fig1a, fig1b)
rm(bold,shy,plot_data)

# 14 Plasticity in Networks ####

# 15 Construct networks and extract metrics ####

library(igraph)
library(assortnet)

# create object for the global metrics

x = rep("B",9)
y = rep("G",9)
z = c(x,y)
x = seq(1,9)
y = seq(1,9)
z2 = c(x,y)

glob_mets = cbind(z,z2)
colnames(glob_mets) = c("group","window")
rm(x,xx,y,z,z2)

# create list of behavioural syndrome scores for each assortment calculation, e.g. each time window

ind_win_d$season = round(ind_win_d$season,0)
names(ind_win_d)
scores = ind_win_d[,c(1,2,3,5,9)]

xx = scores %>%
  group_by(Group,season,Subject) %>%
  summarize(scorem = mean(win_exc))

b1_4 = as.vector(xx$scorem[xx$Group == "B" & xx$season == 1])
b4_6 = as.vector(xx$scorem[xx$Group == "B" & xx$season == 2])
b7_9 = as.vector(xx$scorem[xx$Group == "B" & xx$season == 3])

g1_4 = as.vector(xx$scorem[xx$Group == "G" & xx$season == 1])
g4_6 = as.vector(xx$scorem[xx$Group == "G" & xx$season == 2])
g7_9 = as.vector(xx$scorem[xx$Group == "G" & xx$season == 3])

rm(xx, scores)

netb1 = graph.adjacency(blue_groom1,mode="directed",weighted=TRUE,diag=FALSE)
netb2 = graph.adjacency(blue_groom2,mode="directed",weighted=TRUE,diag=FALSE)
netb3 = graph.adjacency(blue_groom3,mode="directed",weighted=TRUE,diag=FALSE)
netb4 = graph.adjacency(blue_groom4,mode="directed",weighted=TRUE,diag=FALSE)
netb5 = graph.adjacency(blue_groom5,mode="directed",weighted=TRUE,diag=FALSE)
netb6 = graph.adjacency(blue_groom6,mode="directed",weighted=TRUE,diag=FALSE)
netb7 = graph.adjacency(blue_groom7,mode="directed",weighted=TRUE,diag=FALSE)
netb8 = graph.adjacency(blue_groom8,mode="directed",weighted=TRUE,diag=FALSE)
netb9 = graph.adjacency(blue_groom9,mode="directed",weighted=TRUE,diag=FALSE)

names = c("CAS","CON","ELI","GUL","ISA","IZZ","NIC","PEN","ROC","SAR","TIM","WAN")
strength = strength(netb1, mode = "out")
degree = degree(netb1, mode = "out", normalized = TRUE)
window = rep(1,12)
b1_ind_data = cbind(names,window,strength,degree)

strength = strength(netb2, mode = "out")
degree = degree(netb2, mode = "out", normalized = TRUE)
window = rep(2,12)
b2_ind_data = cbind(names,window,strength,degree)

strength = strength(netb3, mode = "out")
degree = degree(netb3, mode = "out", normalized = TRUE)
window = rep(3,12)
b3_ind_data = cbind(names,window,strength,degree)

strength = strength(netb4, mode = "out")
degree = degree(netb4, mode = "out", normalized = TRUE)
window = rep(4,12)
b4_ind_data = cbind(names,window,strength,degree)

strength = strength(netb5, mode = "out")
degree = degree(netb5, mode = "out", normalized = TRUE)
window = rep(5,12)
b5_ind_data = cbind(names,window,strength,degree)

strength = strength(netb6, mode = "out")
degree = degree(netb6, mode = "out", normalized = TRUE)
window = rep(6,12)
b6_ind_data = cbind(names,window,strength,degree)

strength = strength(netb6, mode = "out")
degree = degree(netb6, mode = "out", normalized = TRUE)
window = rep(6,12)
b6_ind_data = cbind(names,window,strength,degree)

strength = strength(netb7, mode = "out")
degree = degree(netb7, mode = "out", normalized = TRUE)
window = rep(7,12)
b7_ind_data = cbind(names,window,strength,degree)

strength = strength(netb8, mode = "out")
degree = degree(netb8, mode = "out", normalized = TRUE)
window = rep(8,12)
b8_ind_data = cbind(names,window,strength,degree)

strength = strength(netb9, mode = "out")
degree = degree(netb9, mode = "out", normalized = TRUE)
window = rep(9,12)
b9_ind_data = cbind(names,window,strength,degree)

blue_data = rbind(b1_ind_data,b2_ind_data,b3_ind_data,b4_ind_data,b5_ind_data,b6_ind_data,
                  b7_ind_data,b8_ind_data,b9_ind_data)
rm(b1_ind_data,b2_ind_data,b3_ind_data,b4_ind_data,b5_ind_data,b6_ind_data,
   b7_ind_data,b8_ind_data,b9_ind_data)

e1 = edge_density(netb1)
a1 = assortment.continuous(blue_groom1,b1_4,weighted=TRUE)$r
e2 = edge_density(netb2)
a2 = assortment.continuous(blue_groom1,b1_4,weighted=TRUE)$r
e3 = edge_density(netb3, loops = FALSE)
a3 = assortment.continuous(blue_groom1,b1_4,weighted=TRUE)$r
e4 = edge_density(netb4, loops = FALSE)
a4 = assortment.continuous(blue_groom1,b1_4,weighted=TRUE)$r
e5 = edge_density(netb5, loops = FALSE)
a5 = assortment.continuous(blue_groom1,b4_6,weighted=TRUE)$r
e6 = edge_density(netb6, loops = FALSE)
a6 = assortment.continuous(blue_groom1,b4_6,weighted=TRUE)$r
e7 = edge_density(netb7, loops = FALSE)
a7 = assortment.continuous(blue_groom1,b7_9,weighted=TRUE)$r
e8 = edge_density(netb8, loops = FALSE)
a8 = assortment.continuous(blue_groom1,b7_9,weighted=TRUE)$r
e9 = edge_density(netb9, loops = FALSE)
a9 = assortment.continuous(blue_groom1,b7_9,weighted=TRUE)$r

assorts = c(a1, a2, a3, a4, a5, a6, a7, a8, a9)
rm(a1, a2, a3, a4, a5, a6, a7, a8, a9)

edges = c(e1, e2, e3, e4, e5, e6, e7, e8, e9)
rm(e1, e2, e3, e4, e5, e6, e7, e8, e9)

# green group

netg1 = graph.adjacency(green_groom1,mode="directed",weighted=TRUE,diag=FALSE)
netg2 = graph.adjacency(green_groom2,mode="directed",weighted=TRUE,diag=FALSE)
netg3 = graph.adjacency(green_groom3,mode="directed",weighted=TRUE,diag=FALSE)
netg4 = graph.adjacency(green_groom4,mode="directed",weighted=TRUE,diag=FALSE)
netg5 = graph.adjacency(green_groom5,mode="directed",weighted=TRUE,diag=FALSE)
netg6 = graph.adjacency(green_groom6,mode="directed",weighted=TRUE,diag=FALSE)
netg7 = graph.adjacency(green_groom7,mode="directed",weighted=TRUE,diag=FALSE)
netg8 = graph.adjacency(green_groom8,mode="directed",weighted=TRUE,diag=FALSE)
netg9 = graph.adjacency(green_groom9,mode="directed",weighted=TRUE,diag=FALSE)

names = c("ANN","ART","DAK","DAN","GEO","HEL","JOA","KER","KRI","LEW","MAC","NOD",
          "OZZ","REB","SIM")

strength = strength(netg1, mode = "out")
degree = degree(netg1, mode = "out", normalized = TRUE)
window = rep(1,15)
g1_ind_data = cbind(names,window,strength,degree)

strength = strength(netg2, mode = "out")
degree = degree(netg2, mode = "out", normalized = TRUE)
window = rep(2,15)
g2_ind_data = cbind(names,window,strength,degree)

strength = strength(netg3, mode = "out")
degree = degree(netg3, mode = "out", normalized = TRUE)
window = rep(3,15)
g3_ind_data = cbind(names,window,strength,degree)

strength = strength(netg4, mode = "out")
degree = degree(netg4, mode = "out", normalized = TRUE)
window = rep(4,15)
g4_ind_data = cbind(names,window,strength,degree)

strength = strength(netg5, mode = "out")
degree = degree(netg5, mode = "out", normalized = TRUE)
window = rep(5,15)
g5_ind_data = cbind(names,window,strength,degree)

strength = strength(netg6, mode = "out")
degree = degree(netg6, mode = "out", normalized = TRUE)
window = rep(6,15)
g6_ind_data = cbind(names,window,strength,degree)

strength = strength(netg7, mode = "out")
degree = degree(netg7, mode = "out", normalized = TRUE)
window = rep(7,15)
g7_ind_data = cbind(names,window,strength,degree)

strength = strength(netg8, mode = "out")
degree = degree(netg8, mode = "out", normalized = TRUE)
window = rep(8,15)
g8_ind_data = cbind(names,window,strength,degree)

strength = strength(netg9, mode = "out")
degree = degree(netg9, mode = "out", normalized = TRUE)
window = rep(9,15)
g9_ind_data = cbind(names,window,strength,degree)

green_data = rbind(g1_ind_data,g2_ind_data,g3_ind_data,g4_ind_data,g5_ind_data,g6_ind_data,
                   g7_ind_data,g8_ind_data,g9_ind_data)
rm(g1_ind_data,g2_ind_data,g3_ind_data,g4_ind_data,g5_ind_data,g6_ind_data,
   g7_ind_data,g8_ind_data,g9_ind_data)

e1 = edge_density(netg1, loops = FALSE)
a1 = assortment.continuous(green_groom1, g1_4,weighted=TRUE)$r
e2 = edge_density(netg2, loops = FALSE)
a2 = assortment.continuous(green_groom2, g1_4,weighted=TRUE)$r
e3 = edge_density(netg3, loops = FALSE)
a3 = assortment.continuous(green_groom3, g1_4,weighted=TRUE)$r
e4 = edge_density(netg4, loops = FALSE)
a4 = assortment.continuous(green_groom4, g1_4,weighted=TRUE)$r
e5 = edge_density(netg5, loops = FALSE)
a5 = assortment.continuous(green_groom5, g4_6,weighted=TRUE)$r
e6 = edge_density(netg6, loops = FALSE)
a6 = assortment.continuous(green_groom6, g4_6,weighted=TRUE)$r
e7 = edge_density(netg7, loops = FALSE)
a7 = assortment.continuous(green_groom7, g7_9,weighted=TRUE)$r
e8 = edge_density(netg8, loops = FALSE)
a8 = assortment.continuous(green_groom8,g7_9,weighted=TRUE)$r
e9 = edge_density(netg9, loops = FALSE)
a9 = assortment.continuous(green_groom9,g7_9,weighted=TRUE)$r

gassorts = c(a1, a2, a3, a4, a5, a6, a7, a8, a9)
rm(a1, a2, a3, a4, a5, a6, a7, a8, a9)

gedges = c(e1, e2, e3, e4, e5, e6, e7, e8, e9)
rm(e1, e2, e3, e4, e5, e6, e7, e8, e9)

asso = c(assorts, gassorts)
dens = c(edges, gedges)

test = as.data.frame(glob_mets)
test$assortivity = asso
test$density = dens

glob_mets = test

rm(test, zz, asso, assorts, degree, dens, edges, gassorts, gedges, names, strength, window)

# 16 Network plots ####


V(netb1)$bs = b1_4
V(netb1)$size <- V(netb1)$bs*5
V(netb1)$label <- NA
E(netb1)$width <- E(netb1)$weight
E(netb1)$arrow.size <- .05
plot(netb1, edge.curved=.1,vertex.color="lightsteelblue1", edge.color="dodgerblue4")

V(netb2)$bs = b1_4
V(netb2)$size <- V(netb2)$bs*5
V(netb2)$label <- NA
E(netb2)$width <- E(netb2)$weight
E(netb2)$arrow.size <- .05
plot(netb2, edge.curved=.1,vertex.color="lightsteelblue1", edge.color="dodgerblue4")

V(netb3)$bs = b1_4
V(netb3)$size <- V(netb3)$bs*5
V(netb3)$label <- NA
E(netb3)$width <- E(netb3)$weight
E(netb3)$arrow.size <- .05
plot(netb3, edge.curved=.1,vertex.color="lightsteelblue1", edge.color="dodgerblue4")

V(netb4)$bs = b1_4
V(netb4)$size <- V(netb4)$bs*5
V(netb4)$label <- NA
E(netb4)$width <- E(netb4)$weight
E(netb4)$arrow.size <- .05
plot(netb4, edge.curved=.1,vertex.color="lightsteelblue1", edge.color="dodgerblue4")


V(netb5)$bs = b4_6
V(netb5)$size <- V(netb5)$bs*5
V(netb5)$label <- NA
E(netb5)$width <- E(netb5)$weight
E(netb5)$arrow.size <- .05
plot(netb5, edge.curved=.1,vertex.color="lightsteelblue1", edge.color="dodgerblue4")

V(netb6)$bs = b4_6
V(netb6)$size <- V(netb6)$bs*5
V(netb6)$label <- NA
E(netb6)$width <- E(netb6)$weight
E(netb6)$arrow.size <- .05
plot(netb6, edge.curved=.1,vertex.color="lightsteelblue1", edge.color="dodgerblue4")


V(netb7)$bs = b7_9
V(netb7)$size <- V(netb7)$bs*5
V(netb7)$label <- NA
E(netb7)$width <- E(netb7)$weight
E(netb7)$arrow.size <- .05
plot(netb7, edge.curved=.1,vertex.color="lightsteelblue1", edge.color="dodgerblue4")

V(netb8)$bs = b7_9
V(netb8)$size <- V(netb8)$bs*5
V(netb8)$label <- NA
E(netb8)$width <- E(netb8)$weight
E(netb8)$arrow.size <- .05
plot(netb8, edge.curved=.1,vertex.color="lightsteelblue1", edge.color="dodgerblue4")

V(netb9)$bs = b7_9
V(netb9)$size <- V(netb9)$bs*5
V(netb9)$label <- NA
E(netb9)$width <- E(netb9)$weight
E(netb9)$arrow.size <- .05
plot(netb9, edge.curved=.1,vertex.color="lightsteelblue1", edge.color="dodgerblue4")

# green

V(netg1)$bs = g1_4
V(netg1)$size <- V(netg1)$bs*5
V(netg1)$label <- NA
E(netg1)$width <- E(netg1)$weight
E(netg1)$arrow.size <- .05
plot(netg1, edge.curved=.1,vertex.color="#c5f1c5", edge.color="green4")

V(netg2)$bs = g1_4
V(netg2)$size <- V(netg2)$bs*5
V(netg2)$label <- NA
E(netg2)$width <- E(netg2)$weight
E(netg2)$arrow.size <- .05
plot(netg2, edge.curved=.1,vertex.color="#c5f1c5", edge.color="green4")

V(netg3)$bs = g1_4
V(netg3)$size <- V(netg3)$bs*5
V(netg3)$label <- NA
E(netg3)$width <- E(netg3)$weight
E(netg3)$arrow.size <- .05
plot(netg3, edge.curved=.1,vertex.color="#c5f1c5", edge.color="green4")

V(netg4)$bs = g1_4
V(netg4)$size <- V(netg4)$bs*5
V(netg4)$label <- NA
E(netg4)$width <- E(netg4)$weight
E(netg4)$arrow.size <- .05
plot(netg4, edge.curved=.1,vertex.color="#c5f1c5", edge.color="green4")


V(netg5)$bs = g4_6
V(netg5)$size <- V(netg5)$bs*5
V(netg5)$label <- NA
E(netg5)$width <- E(netg5)$weight
E(netg5)$arrow.size <- .05
plot(netg5, edge.curved=.1,vertex.color="#c5f1c5", edge.color="green4")

V(netg6)$bs = g4_6
V(netg6)$size <- V(netg6)$bs*5
V(netg6)$label <- NA
E(netg6)$width <- E(netg6)$weight
E(netg6)$arrow.size <- .05
plot(netg6, edge.curved=.1,vertex.color="#c5f1c5", edge.color="green4")


V(netg7)$bs = g7_9
V(netg7)$size <- V(netg7)$bs*5
V(netg7)$label <- NA
E(netg7)$width <- E(netg7)$weight
E(netg7)$arrow.size <- .05
plot(netg7, edge.curved=.1,vertex.color="#c5f1c5", edge.color="green4")

V(netg8)$bs = g7_9
V(netg8)$size <- V(netg8)$bs*5
V(netg8)$label <- NA
E(netg8)$width <- E(netg8)$weight
E(netg8)$arrow.size <- .05
plot(netg8, edge.curved=.1,vertex.color="#c5f1c5", edge.color="green4")

V(netg9)$bs = g7_9
V(netg9)$size <- V(netg9)$bs*5
V(netg9)$label <- NA
E(netg9)$width <- E(netg9)$weight
E(netg9)$arrow.size <- .05
plot(netg9, edge.curved=.1,vertex.color="#c5f1c5", edge.color="green4")

# 17 Add metrics to window data ####

library(dplyr)

# combine network data

blue_data = as.data.frame(blue_data)
blue_data$Group = "B"
green_data = as.data.frame(green_data)
green_data$Group = "G"
network_d = rbind(blue_data, green_data)

names(network_d)
names(ind_win_d)
network_d = network_d[order(network_d$window, network_d$Group, network_d$names),]

test = network_d %>%
  mutate(Subject = as.character(names),
         time_window = as.character(window),
         strength = as.numeric(as.character(strength)),
         degree = as.numeric(as.character(degree)),
         Group = as.character(Group))
network_d = test

test = ind_win_d %>%
  mutate(Subject = as.character(Subject),
         Group = as.character(Group),
         time_window = as.character(time_window),
         win_rank = as.numeric(as.character(win_rank)),
         win_exc = as.numeric(as.character(win_exc)),
         time_win_min_temp = as.numeric(as.character(time_win_min_temp)),
         time_win_fruit = as.numeric(as.character(time_win_fruit)),
         time_win_anthro = as.numeric(as.character(time_win_anthro)),
         season = as.character(season))
ind_win_d = test
rm(test, blue_data, green_data)

# let's try a left join to add on the network metrics

names(network_d)
network_d = network_d[,c(3:7)]

test = ind_win_d %>%
  left_join(network_d, by = c("Subject", "time_window"))
test$Group.y = NULL
test$group = test$Group.x
test$Group.x = NULL

names(test)

ind_win_d = test[,c(1,11,2,3:10)]

ind_win_d$seasonx = ifelse(ind_win_d$season == 1, 1, 0)
ind_win_d$seasony = ifelse(ind_win_d$seasonx == 1, "Mating", "Non-mating")

# last thing to add is the sex of the individual

names(scan_d)
xx = as.data.frame(table(scan_d$Subject, scan_d$Sex))
xx = filter(xx, !Freq == 0)
xx$Freq = NULL
colnames(xx) = c("Subject","sex")

test = ind_win_d %>%
  left_join(xx, by = "Subject")

ind_win_d = test

rm(test,xx)

# 18 Finalise dataset and variables ####

pd2 = ind_win_d
names(pd)
names(pd2)

test = pd2 %>%
  mutate(subject_ID = as.character(Subject),
         group = as.character(group),
         sex = as.character(sex),
         seasonx = as.numeric(as.character(seasonx)),
         rank = as.numeric(as.character(win_rank)),
         behavioural_syndrome = as.numeric(as.character(win_exc)),
         temp = as.numeric(as.character(time_win_min_temp)),
         food = as.numeric(as.character(time_win_fruit)),
         anthro = as.numeric(as.character(time_win_anthro)))

names(test)
names(pd)
test = test[,c(15,10,11,3,2,14,12,13,20,18,19,16,17)]
pd2 = test
rm(test)

pd2$syndrome_z <- (pd2$behavioural_syndrome - mean(pd2$behavioural_syndrome))/sd(pd2$behavioural_syndrome)
range(pd2$syndrome_z)

range(pd2$seasonx)

pd2$anthro_z <- (pd2$anthro - mean(pd2$anthro))/sd(pd2$anthro)
range(pd2$anthro_z)

pd2$food_z <- (pd2$food - mean(pd2$food))/sd(pd2$food)
range(pd2$food_z)

pd2$temp_z <- (pd2$temp - mean(pd2$temp))/sd(pd2$temp)
range(pd2$temp_z)

pd2$rank_z <- (pd2$rank - mean(pd2$rank))/sd(pd2$rank)
range(pd2$rank_z)

# remove time_window 4 for Blue group, rows 82 to 93

test = pd2[-c(82:93),]
pd2 = test

# 19 Create pairs plots to check predictor variables ####

names(pd2)

covees = pd2[,c(5,6,7,14:18)]
colnames(covees) = c("group", "sex", "seasonx","syndrome_z","anthro_z",           
                     "food_z","temp_z","rank_z")

library(car)

mod = lm(degree ~ seasonx + syndrome_z + anthro_z + food_z + temp_z +  
           rank_z + sex + group, data = pd2)

vif(mod)

' seasonx syndrome_z   anthro_z     food_z     temp_z     rank_z        sex      group 
  2.008551   2.832602  14.226437   3.694265   1.253264   2.969974   4.396870  16.418126 '

mod = lm(degree ~ seasonx + syndrome_z + anthro_z + food_z + temp_z +  
           rank_z + sex, data = pd2)

vif(mod)

'   seasonx syndrome_z   anthro_z     food_z     temp_z     rank_z        sex 
  1.713562   2.814005   2.817399   3.135801   1.247714   2.969237   4.378225'

mod = lm(strength ~ seasonx + syndrome_z + anthro_z + food_z + temp_z +  
           rank_z + sex, data = pd2)

vif(mod)

'   seasonx syndrome_z   anthro_z     food_z     temp_z     rank_z        sex 
  1.713562   2.814005   2.817399   3.135801   1.247714   2.969237   4.378225 '

library(ggplot2)
library(GGally)
windows()
ggpairs(covees) # all <=0.7

rm(covees)

# Check for model complexity

library(lme4)

names(pd2)

modc = lmer(degree ~ syndrome_z*seasonx +  syndrome_z*anthro_z + 
              syndrome_z*food_z + syndrome_z*temp_z +  
              rank_z + sex +
              (seasonx+anthro_z+food_z+temp_z|subject_ID), data = pd2)

round(length(residuals(modc))/(length(fixef(modc))+length(summary(modc)$varcor)+1),1)
#16.5 so no issue with complexity

# 20 Degree model ####

library(brms)

hist(pd2$degree)
range(pd2$degree)
pd2$degree = pd2$degree + 0.001
range(pd2$degree)

hist(pd2$strength)
range(pd2$strength)
pd2$strength = pd2$strength + 0.001
range(pd2$degree)

mprior3 = get_prior(degree ~ seasonx*syndrome_z + anthro_z*syndrome_z + 
                      food_z*syndrome_z + temp_z*syndrome_z + rank_z + sex +
                      (seasonx+anthro_z+food_z+temp_z|subject_ID), 
                    data = pd2, family = Beta(link = "logit", link_phi = "log"))

mprior3$prior[2:12] <- "normal(0,1)"

make_stancode(degree ~ seasonx*syndrome_z + anthro_z*syndrome_z + 
                food_z*syndrome_z + temp_z*syndrome_z + rank_z + sex +
                (seasonx+anthro_z+food_z+temp_z|subject_ID), 
              data = pd2, family = Beta(link = "logit", link_phi = "log"), prior = mprior3)

mod2 = brm(degree ~ seasonx*syndrome_z + anthro_z*syndrome_z + 
             food_z*syndrome_z + temp_z*syndrome_z + rank_z + sex +
             (seasonx+anthro_z+food_z+temp_z|subject_ID), 
           data = pd2, family = Beta(link = "logit", link_phi = "log"), prior = mprior3, 
           chains = 3, cores = 3, iter=5000, warmup = 2000, thin = 2, 
           control = list(adapt_delta = 0.95), sample_prior = "yes")

summary(mod2)
plot(mod2)

# all diagnostics ok - remove 'nonsign' interactions

mprior4 = get_prior(degree ~ temp_z*syndrome_z + seasonx*syndrome_z + anthro_z*syndrome_z +
                      food_z + rank_z + sex +
                      (seasonx+anthro_z+food_z+temp_z|subject_ID), 
                    data = pd2, family = Beta(link = "logit", link_phi = "log"))

mprior4$prior[2:11] <- "normal(0,1)"

make_stancode(degree ~ temp_z*syndrome_z + seasonx*syndrome_z + anthro_z*syndrome_z +
                food_z + rank_z + sex +
                (seasonx+anthro_z+food_z+temp_z|subject_ID), 
              data = pd2, family = Beta(link = "logit", link_phi = "log"), prior = mprior4)

red_mod2 = brm(degree ~ temp_z*syndrome_z + seasonx*syndrome_z + anthro_z*syndrome_z +
                 food_z + rank_z + sex +
                 (seasonx+anthro_z+food_z+temp_z|subject_ID), 
               data = pd2, family = Beta(link = "logit", link_phi = "log"), prior = mprior4, 
               chains = 3, cores = 3, iter=5000, warmup = 2000, thin = 2, 
               control = list(adapt_delta = 0.95), sample_prior = "yes")


summary(red_mod2)
plot(red_mod2)

# 21 PPC of degree model ####

#save.image("C:/Users/Megaport/OneDrive/Project_data/social_comp/analysis/results.Rdata")

library(ggplot2)

pp_check(red_mod2)

length(pd2$syndrome_z[pd2$syndrome_z <= -1]) #29
length(pd2$syndrome_z[pd2$syndrome_z > -1 & pd2$syndrome_z <1]) #161
length(pd2$syndrome_z[pd2$syndrome_z >=1]) #41

# so we need a ppc for relationship between:
# seasonality
# sex
# excitable*temp

pred_vals2 = predict(red_mod2, summary = FALSE)
predictions2 = as.data.frame(t(pred_vals2))

mean = rowMeans(predictions2)

fit_mod2 <- 
  fitted(red_mod2) %>%
  as_tibble() %>%
  bind_cols(pd2)

fit_mod2$predicted = as.numeric(mean)

# create new plotting object with 100 draws from posterior

# now proper plot by randomly selecting 100 draws from the posterior

xx = floor(runif(100, 1, 3000))

# extract those from our predictions object

zz = predictions2[,xx]

colnames(zz) = c("pred1", "pred2", "pred3", "pred4", "pred5", "pred6", "pred7", "pred8", "pred9", "pred10",
                 "pred11", "pred12", "pred13", "pred14", "pred15", "pred16", "pred17", "pred18", "pred19", "pred20",
                 "pred21", "pred22", "pred23", "pred24", "pred25", "pred26", "pred27", "pred28", "pred29", "pred30",
                 "pred31", "pred32", "pred33", "pred34", "pred35", "pred36", "pred37", "pred38", "pred39", "pred40",
                 "pred41", "pred42", "pred43", "pred44", "pred45", "pred46", "pred47", "pred48", "pred49", "pred50",
                 "pred51", "pred52", "pred53", "pred54", "pred55", "pred56", "pred57", "pred58", "pred59", "pred60",
                 "pred61", "pred62", "pred63", "pred64", "pred65", "pred66", "pred67", "pred68", "pred69", "pred70",
                 "pred71", "pred72", "pred73", "pred74", "pred75", "pred76", "pred77", "pred78", "pred79", "pred80",
                 "pred81", "pred82", "pred83", "pred84", "pred85", "pred86", "pred87", "pred88", "pred89", "pred90",
                 "pred91", "pred92", "pred93", "pred94", "pred95", "pred96", "pred97", "pred98", "pred99", "pred100")

# create new plotting object

names(fit_mod2)
plot_data = fit_mod2[, c(1,3,4,10,11,13,18,14,23,7)]
names(plot_data)
colnames(plot_data) = c("Estimate", "lower", "upper", "sex", "season","anthro", 
                        "behavioural_syndrome","temp","mean_prediction", "observed")
plot_data$season2 = ifelse(plot_data$season == 1, "Mating", "Non-mating")

plot_data = cbind(plot_data, zz)

# for the season and sex plots might need a sub plot_data that summarises the result
# there's an interaction with boldness so need to take that into account too

round(mean(plot_data$behavioural_syndrome),3)
length(plot_data$behavioural_syndrome[plot_data$behavioural_syndrome <=0])
length(plot_data$behavioural_syndrome[plot_data$behavioural_syndrome >0])

bold1 = filter(plot_data, behavioural_syndrome <=0)
bold1$behavioural_syndrome3 = "Shy"
bold2 = filter(plot_data, behavioural_syndrome >0)
bold2$behavioural_syndrome3 = "Bold"
test=rbind(bold1,bold2)
plot_data = test
rm(bold1, bold2, test)

names(plot_data)
levels(as.factor(plot_data$behavioural_syndrome))
bold = filter(plot_data, behavioural_syndrome3 == "Bold")
shy = filter(plot_data, behavioural_syndrome3 == "Shy")


bold_mean_season1 = mean(bold$Estimate[bold$season2 == "Mating"])
bold_sd_season1 = sd(bold$Estimate[bold$season2 == "Mating"])
bold_lower_season1 = mean(bold$lower[bold$season2 == "Mating"])
bold_upper_season1 = mean(bold$upper[bold$season2 == "Mating"])
bold_mean_season2 = mean(bold$Estimate[bold$season2 == "Non-mating"])
bold_sd_season2 = sd(bold$Estimate[bold$season2 == "Non-mating"])
bold_lower_season2 = mean(bold$lower[bold$season2 == "Non-mating"])
bold_upper_season2 = mean(bold$upper[bold$season2 == "Non-mating"])

season = c("Mating", "Non-mating")
means = c(bold_mean_season1, bold_mean_season2)
lowers = c(bold_lower_season1, bold_lower_season2)
uppers = c(bold_upper_season1, bold_upper_season2)
sds = c(bold_sd_season1,bold_sd_season2)

bold_data = cbind(means,lowers,uppers,sds)
bold_data = as.data.frame(bold_data)
bold_data$season = season

shy_mean_season1 = mean(shy$Estimate[shy$season2 == "Mating"])
shy_sd_season1 = sd(shy$Estimate[shy$season2 == "Mating"])
shy_lower_season1 = mean(shy$lower[shy$season2 == "Mating"])
shy_upper_season1 = mean(shy$upper[shy$season2 == "Mating"])
shy_mean_season2 = mean(shy$Estimate[shy$season2 == "Non-mating"])
shy_sd_season2 = sd(shy$Estimate[shy$season2 == "Non-mating"])
shy_lower_season2 = mean(shy$lower[shy$season2 == "Non-mating"])
shy_upper_season2 = mean(shy$upper[shy$season2 == "Non-mating"])

season = c("Mating", "Non-mating")
means = c(shy_mean_season1, shy_mean_season2)
lowers = c(shy_lower_season1, shy_lower_season2)
uppers = c(shy_upper_season1, shy_upper_season2)
sds = c(shy_sd_season1,shy_sd_season2)

shy_data = cbind(means,lowers,uppers,sds)
shy_data = as.data.frame(shy_data)
shy_data$season = season

boldseason_ppc = ggplot(bold)+
  geom_jitter(aes(x=season2, y=observed+0.025),height=0.05, alpha=1/10)+
  geom_jitter(aes(x=season2, y=mean_prediction), height=0.05, alpha=1/10, colour = "red")+
  geom_errorbar(data=bold_data, aes(x=season,ymin=means-sds, ymax=means+sds), width=0.5)+
  geom_point(data=bold_data, aes(x=season,y=means))+
  labs(x = "Season", y = "Network degree (normalised)")
boldseason_ppc = boldseason_ppc+coord_cartesian(ylim = c(0, 1))
boldseason_ppc = boldseason_ppc+theme_classic()
#season_ppc = season_ppc+theme(legend.position="none")
boldseason_ppc = boldseason_ppc+ggtitle("Bold")
boldseason_ppc

shyseason_ppc = ggplot(shy)+
  geom_jitter(aes(x=season2, y=observed+0.025),height=0.05, alpha=1/10)+
  geom_jitter(aes(x=season2, y=mean_prediction), height=0.05, alpha=1/10, colour = "red")+
  geom_errorbar(data=shy_data, aes(x=season,ymin=means-sds, ymax=means+sds), width=0.5)+
  geom_point(data=shy_data, aes(x=season,y=means))+
  labs(x = "Season", y = "")
shyseason_ppc = shyseason_ppc+coord_cartesian(ylim = c(0, 1))
shyseason_ppc = shyseason_ppc+theme_classic()
#season_ppc = season_ppc+theme(legend.position="none")
shyseason_ppc = shyseason_ppc+ggtitle("Shy")
shyseason_ppc

degseason_ppc = boldseason_ppc+shyseason_ppc
degseason_ppc

# sex plot

mean_male = mean(plot_data$Estimate[plot_data$sex == "M"])
sd_male = sd(plot_data$Estimate[plot_data$sex == "M"])
lower_male = mean(plot_data$lower[plot_data$sex == "M"])
upper_male = mean(plot_data$upper[plot_data$sex == "M"])
mean_fem = mean(plot_data$Estimate[plot_data$sex == "F"])
sd_fem = sd(plot_data$Estimate[plot_data$sex == "F"])
lower_fem = mean(plot_data$lower[plot_data$sex == "F"])
upper_fem = mean(plot_data$upper[plot_data$sex == "F"])

sex = c("M", "F")
means2 = c(mean_male, mean_fem)
lowers2 = c(lower_male, lower_fem)
uppers2 = c(upper_male, upper_fem)
sds2 = c(sd_male,sd_fem)

plot_data3 = cbind(means2,lowers2,uppers2,sds2)
plot_data3 = as.data.frame(plot_data3)
plot_data3$sex = sex

sex_ppc = ggplot(plot_data)+
  geom_jitter(aes(x=sex, y=observed+0.025),height=0.05, alpha=1/10)+
  geom_jitter(aes(x=sex, y=mean_prediction), height=0.05, alpha=1/10, colour = "red")+
  geom_errorbar(data=plot_data2, aes(x=sex,ymin=means2-sds2, ymax=means2+sds2), width=0.5)+
  geom_point(data=plot_data2, aes(x=sex,y=means2))+
  labs(x = "Sex", y = "Network degree (normalised)")
sex_ppc = sex_ppc+theme_classic()
#season_ppc = season_ppc+theme(legend.position="none")
sex_ppc = sex_ppc+ggtitle("(a)")
sex_ppc

# now for the temperature interaction plot

names(plot_data)
temp_ppc = ggplot(plot_data)+
  geom_jitter(aes(x=temp, y=observed+0.025),height=0.05, alpha=1/10)+
  geom_jitter(aes(x=temp, y=mean_prediction), height=0.05, alpha=1/10, colour = "red")+
  geom_smooth(aes(x=temp, y=observed), method = "glm", 
              method.args = list(family = binomial(link='logit')),
              se = TRUE)+
  geom_smooth(aes(x=temp, y=mean_prediction, colour = "red"), method = "glm", 
              method.args = list(family = binomial(link='logit')), 
              se = TRUE)+
  facet_wrap(~behavioural_syndrome3)+labs(x = "Mean minimum temperature (C)",
                                          y = "Network degree (normalised)")
temp_ppc = temp_ppc+theme_classic()
temp_ppc = temp_ppc+theme(legend.position="none")
temp_ppc = temp_ppc+ggtitle("(b)")
temp_ppc

# anthro ppc

names(plot_data)
anthro_ppc = ggplot(plot_data)+
  geom_jitter(aes(x=anthro, y=observed+0.025),height=0.05, alpha=1/10)+
  geom_jitter(aes(x=anthro, y=mean_prediction), height=0.05, alpha=1/10, colour = "red")+
  geom_smooth(aes(x=anthro, y=observed), method = "glm", 
              method.args = list(family = binomial(link='logit')),
              se = TRUE)+
  geom_smooth(aes(x=anthro, y=mean_prediction, colour = "red"), method = "glm", 
              method.args = list(family = binomial(link='logit')), 
              se = TRUE)+
  facet_wrap(~behavioural_syndrome3)+labs(x = "Mean anthropogenic index",
                                          y = "Network degree (normalised)")
anthro_ppc = anthro_ppc+theme_classic()
anthro_ppc = anthro_ppc+theme(legend.position="none")
anthro_ppc = anthro_ppc+ggtitle("(d)")
anthro_ppc

rm(bs_ppc, bs_ppc2)

rm(plot_data2, plot_data3, season_ppc, sex_ppc, time_ppc)
rm(lower_fem, lower_male, lower_season1, lower_season2, lowers, lowers2)
rm(upper_fem, upper_male, upper_season1, upper_season2, uppers, uppers2)
rm(mean_fem, mean_male, mean_season1, mean_season2, means, means2)
rm(i, medians, sd_fem, sd_male, sd_season1, sd_season2, sds, sds2, season, sex, xx)
rm(bold, bold_data, boldseason_ppc, degseason_ppc, shy, shy_data, shyseason_ppc)
rm(shy_lower_season1,shy_lower_season2, 
   shy_mean_season1,shy_mean_season2,shy_sd_season1,shy_sd_season2,shy_upper_season1, 
   shy_upper_season2)
rm(bold_lower_season1,
   bold_lower_season2,bold_mean_season1,bold_mean_season2,bold_sd_season1,bold_sd_season2,   
   bold_upper_season1,bold_upper_season2)

# 22 Model Estimates and Table 2 ####

#save.image("C:/Users/Megaport/OneDrive/Project_data/social_comp/analysis/results.Rdata")

round(posterior_summary(red_mod2, probs = c(0.05, 0.95)),3)[c(1:16),]

m2.post <- posterior_samples(red_mod2)
names(m2.post)

sum(m2.post$b_syndrome_z < 0) / length(m2.post$b_syndrome_z)
sum(m2.post$b_seasonx > 0) / length(m2.post$b_seasonx)
sum(m2.post$b_anthro_z > 0) / length(m2.post$b_anthro_z)
sum(m2.post$b_food_z > 0) / length(m2.post$b_food_z)
sum(m2.post$b_temp_z < 0) / length(m2.post$b_temp_z)
sum(m2.post$b_rank_z < 0) / length(m2.post$b_rank_z)
sum(m2.post$b_sexM < 0) / length(m2.post$b_sexM)
sum(m2.post$'b_syndrome_z:seasonx' > 0) / length(m2.post$'b_syndrome_z:seasonx')
sum(m2.post$'b_syndrome_z:anthro_z' > 0) / length(m2.post$'b_syndrome_z:anthro_z')
sum(m2.post$'b_temp_z:syndrome_z' < 0) / length(m2.post$'b_temp_z:syndrome_z')

# also the original model estimates for the supplementary

round(posterior_summary(mod2, probs = c(0.05, 0.95)),3)[c(1:16),]

m2.postb <- posterior_samples(mod2)
names(m2.postb)

sum(m2.postb$b_syndrome_z < 0) / length(m2.postb$b_syndrome_z)
sum(m2.postb$b_seasonx > 0) / length(m2.postb$b_seasonx)
sum(m2.postb$b_anthro_z > 0) / length(m2.postb$b_anthro_z)
sum(m2.postb$b_food_z > 0) / length(m2.postb$b_food_z)
sum(m2.postb$b_temp_z < 0) / length(m2.postb$b_temp_z)
sum(m2.postb$b_rank_z < 0) / length(m2.postb$b_rank_z)
sum(m2.postb$b_sexM < 0) / length(m2.postb$b_sexM)
sum(m2.postb$'b_seasonx:syndrome_z' > 0) / length(m2.postb$'b_seasonx:syndrome_z')
sum(m2.postb$'b_syndrome_z:anthro_z' > 0) / length(m2.postb$'b_syndrome_z:anthro_z')
sum(m2.postb$'b_syndrome_z:food_z' > 0) / length(m2.postb$'b_syndrome_z:food_z')
sum(m2.postb$'b_syndrome_z:temp_z' < 0) / length(m2.postb$'b_syndrome_z:temp_z')

# 23 Plot Degree Result ####

# so now we want to use the draws to plot between bold and shy

# probably best to split the data, make two plots and put them together

shy = filter(plot_data, behavioural_syndrome3 == "Shy")
bold = filter(plot_data, behavioural_syndrome3 == "Bold")

# season*excitability

mean_season1 = mean(bold$Estimate[bold$season2 == "Mating"])
sd_season1 = sd(bold$Estimate[bold$season2 == "Mating"])
lower_season1 = mean(bold$lower[bold$season2 == "Mating"])
upper_season1 = mean(bold$upper[bold$season2 == "Mating"])
mean_season2 = mean(bold$Estimate[bold$season2 == "Non-mating"])
sd_season2 = sd(bold$Estimate[bold$season2 == "Non-mating"])
lower_season2 = mean(bold$lower[bold$season2 == "Non-mating"])
upper_season2 = mean(bold$upper[bold$season2 == "Non-mating"])

season = c("Mating", "Non-mating")
means = c(mean_season1, mean_season2)
lowers = c(lower_season1, lower_season2)
uppers = c(upper_season1, upper_season2)
sds = c(sd_season1,sd_season2)

bold2 = cbind(means,lowers,uppers,sds)
bold2 = as.data.frame(bold2)
bold2$season = season

mean_season1 = mean(shy$Estimate[shy$season2 == "Mating"])
sd_season1 = sd(shy$Estimate[shy$season2 == "Mating"])
lower_season1 = mean(shy$lower[shy$season2 == "Mating"])
upper_season1 = mean(shy$upper[shy$season2 == "Mating"])
mean_season2 = mean(shy$Estimate[shy$season2 == "Non-mating"])
sd_season2 = sd(shy$Estimate[shy$season2 == "Non-mating"])
lower_season2 = mean(shy$lower[shy$season2 == "Non-mating"])
upper_season2 = mean(shy$upper[shy$season2 == "Non-mating"])

season = c("Mating", "Non-mating")
means = c(mean_season1, mean_season2)
lowers = c(lower_season1, lower_season2)
uppers = c(upper_season1, upper_season2)
sds = c(sd_season1,sd_season2)

shy2 = cbind(means,lowers,uppers,sds)
shy2 = as.data.frame(shy2)
shy2$season = season

season2 = c(1,0)
bold2$season2=season2
shy2$season2=season2


str(bold)
fig2a = ggplot(bold)+
  geom_line(colour = "black", aes(x=season, y=observed),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred2),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred3),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred4),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred5),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred6),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred7),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred8),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred9),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred10),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred11),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred12),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred13),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred14),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred15),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred16),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred17),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred18),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred19),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred20),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred21),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred22),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred23),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred24),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred25),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred26),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred27),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred28),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred29),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred30),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred31),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred32),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred33),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred34),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred35),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred36),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred37),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred38),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred39),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred40),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred41),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred42),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred43),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred44),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred45),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred46),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred47),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred48),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred49),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred50),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred51),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred52),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred53),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred54),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred55),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred56),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred57),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred58),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred59),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred60),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred61),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred63),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred64),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred65),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred66),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred67),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred68),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred69),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred70),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred71),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred72),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred73),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred74),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred75),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred76),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred77),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred78),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred79),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred80),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred81),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred82),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred83),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred84),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred85),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred86),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred87),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred88),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred89),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred90),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred91),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred92),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred93),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred94),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred95),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred96),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred97),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred98),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred99),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred100),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_errorbar(data=bold2, aes(x=season2,ymin=means-sds, ymax=means+sds), width=0.5)+
  geom_point(data=bold2, aes(x=season2,y=means))
fig2a = fig2a + coord_cartesian(ylim = c(0, 0.5))
fig2a = fig2a+theme_classic()
fig2a = fig2a + labs(x = "",
                     y = "Network degree (normalised")
#
fig2a = fig2a+ggtitle("Bold")   
fig2a  

names(shy)
fig2b = ggplot(shy)+
  geom_line(colour = "black", aes(x=season, y=observed),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred2),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred3),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred4),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred5),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred6),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred7),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred8),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred9),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred10),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred11),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred12),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred13),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred14),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred15),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred16),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred17),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred18),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred19),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred20),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred21),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred22),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred23),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred24),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred25),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred26),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred27),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred28),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred29),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred30),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred31),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred32),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred33),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred34),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred35),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred36),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred37),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred38),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred39),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred40),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred41),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred42),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred43),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred44),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred45),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred46),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred47),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred48),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred49),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred50),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred51),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred52),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred53),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred54),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred55),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred56),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred57),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred58),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred59),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred60),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred61),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred63),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred64),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred65),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred66),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred67),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred68),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred69),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred70),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred71),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred72),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred73),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred74),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred75),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred76),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred77),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred78),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred79),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred80),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred81),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred82),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred83),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred84),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred85),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred86),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred87),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred88),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred89),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred90),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred91),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred92),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred93),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred94),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred95),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred96),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred97),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred98),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred99),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred100),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_errorbar(data=shy2, aes(x=season2,ymin=means-sds, ymax=means+sds), width=0.5)+
  geom_point(data=shy2, aes(x=season2,y=means))
fig2b = fig2b + coord_cartesian(ylim = c(0, 0.5))
fig2b = fig2b + labs(x = "",
                     y = "")
fig2b = fig2b+theme_classic()
fig2b = fig2b+ggtitle("Shy")   
fig2b

library(patchwork)

fig2 = fig2a+fig2b

# fig2b temperature

names(bold)
fig2c = ggplot(bold)+
  geom_line(colour = "black", aes(x=temp, y=observed),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred2),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred3),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred4),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred5),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred6),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred7),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred8),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred9),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred10),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred11),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred12),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred13),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred14),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred15),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred16),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred17),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred18),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred19),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred20),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred21),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred22),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred23),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred24),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred25),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred26),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred27),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred28),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred29),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred30),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred31),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred32),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred33),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred34),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred35),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred36),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred37),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred38),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred39),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred40),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred41),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred42),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred43),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred44),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred45),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred46),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred47),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred48),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred49),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred50),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred51),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred52),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred53),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred54),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred55),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred56),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred57),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred58),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred59),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred60),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred61),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred63),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred64),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred65),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred66),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred67),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred68),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred69),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred70),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred71),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred72),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred73),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred74),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred75),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred76),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred77),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred78),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred79),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred80),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred81),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred82),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred83),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred84),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred85),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred86),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred87),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred88),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred89),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred90),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred91),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred92),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred93),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred94),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred95),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred96),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred97),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred98),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred99),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred100),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)
fig2c = fig2c + coord_cartesian(xlim = c(1.7, 12.5), ylim = c(0, 0.5))
fig2c = fig2c + labs(x = "",
                     y = "Network degree (normalised)")
fig2c = fig2c+theme_classic()
fig2c = fig2c+ggtitle("Bold")   
fig2c   

fig2d = ggplot(shy)+
  geom_line(colour = "black", aes(x=temp, y=observed),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred2),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred3),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred4),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred5),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred6),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred7),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred8),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred9),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred10),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred11),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred12),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred13),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred14),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred15),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred16),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred17),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred18),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred19),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred20),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred21),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred22),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred23),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred24),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred25),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred26),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred27),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred28),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred29),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred30),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred31),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred32),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred33),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred34),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred35),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred36),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred37),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred38),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred39),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred40),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred41),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred42),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred43),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred44),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred45),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred46),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred47),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred48),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred49),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred50),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred51),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred52),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred53),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred54),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred55),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred56),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred57),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred58),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred59),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred60),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred61),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred63),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred64),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred65),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred66),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred67),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred68),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred69),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred70),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred71),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred72),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred73),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred74),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred75),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred76),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred77),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred78),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred79),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred80),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred81),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred82),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred83),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred84),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred85),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred86),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred87),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred88),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred89),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred90),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred91),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred92),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred93),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred94),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred95),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred96),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred97),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred98),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred99),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=temp, y=pred100),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)
fig2d = fig2d + coord_cartesian(xlim = c(1.7, 12.5), ylim = c(0, 0.5))
fig2d = fig2d + labs(x = "",
                     y = "")
fig2d = fig2d+theme_classic()
fig2d = fig2d+ggtitle("Shy")   
fig2d   

range(plot_data$temp)
range(plot_data$mean_prediction)

library(patchwork)

figure2b = fig2c+fig2d

# anthro plot

names(bold)
fig2e = ggplot(bold)+
  geom_line(colour = "black", aes(x=anthro, y=observed),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred2),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred3),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred4),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred5),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred6),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred7),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred8),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred9),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred10),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred11),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred12),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred13),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred14),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred15),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred16),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred17),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred18),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred19),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred20),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred21),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred22),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred23),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred24),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred25),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred26),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred27),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred28),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred29),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred30),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred31),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred32),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred33),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred34),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred35),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred36),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred37),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred38),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred39),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred40),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred41),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred42),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred43),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred44),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred45),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred46),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred47),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred48),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred49),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred50),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred51),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred52),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred53),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred54),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred55),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred56),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred57),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred58),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred59),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred60),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred61),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred63),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred64),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred65),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred66),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred67),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred68),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred69),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred70),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred71),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred72),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred73),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred74),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred75),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred76),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred77),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred78),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred79),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred80),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred81),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred82),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred83),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred84),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred85),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred86),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred87),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred88),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred89),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred90),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred91),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred92),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred93),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred94),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred95),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred96),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred97),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred98),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred99),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred100),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)
fig2e = fig2e + coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.5))
fig2e = fig2e + labs(x = "",
                     y = "Network degree (normalised)")
fig2e = fig2e+theme_classic()
fig2e = fig2e+ggtitle("Bold")   
fig2e   

fig2f = ggplot(shy)+
  geom_line(colour = "black", aes(x=anthro, y=observed),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred2),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred3),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred4),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred5),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred6),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred7),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred8),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred9),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred10),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred11),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred12),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred13),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred14),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred15),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred16),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred17),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred18),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred19),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred20),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred21),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred22),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred23),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred24),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred25),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred26),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred27),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred28),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred29),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred30),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred31),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred32),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred33),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred34),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred35),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred36),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred37),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred38),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred39),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred40),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred41),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred42),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred43),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred44),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred45),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred46),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred47),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred48),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred49),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred50),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred51),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred52),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred53),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred54),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred55),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred56),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred57),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred58),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred59),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred60),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred61),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred63),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred64),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred65),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred66),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred67),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred68),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred69),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred70),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred71),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred72),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred73),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred74),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred75),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred76),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred77),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred78),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred79),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred80),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred81),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred82),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred83),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred84),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred85),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred86),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred87),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred88),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred89),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred90),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred91),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred92),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred93),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred94),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred95),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred96),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred97),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred98),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred99),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=anthro, y=pred100),
            stat = "smooth", method = "glm", method.args = list(family = binomial(link='logit')),
            se = FALSE, size=1, alpha=0.1)
fig2f = fig2f + coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.5))
fig2f = fig2f + labs(x = "",
                     y = "")
fig2f = fig2f+theme_classic()
fig2f = fig2f+ggtitle("Shy")   
fig2f   

range(plot_data$temp)
range(plot_data$mean_prediction)

library(patchwork)

figure2c = fig2e+fig2f

rm(fig2,fig2a,fig2b,bold,mod,modc,shy,smode,zz,mean)
rm(anthro_ppc, bold2, fig2c, fig2d, fig2e, fig2f, figure2b, figure2c, shy2)
rm(lower_season1, lower_season2, lowers, mean_season1, mean_season2, means)
rm(sd_season1, sd_season2, sds, season, season2, upper_season1, upper_season2)
rm(uppers)

# 24 Strength model ####

#save.image("C:/Users/Megaport/OneDrive/Project_data/social_comp/analysis/results.Rdata")

hist(pd2$strength)
range(pd2$strength)
pd2$logS = log10(pd2$strength+1)
hist(pd2$logS)

mprior5 = get_prior(strength ~ seasonx*syndrome_z + anthro_z*syndrome_z + 
                      food_z*syndrome_z + temp_z*syndrome_z + rank_z + sex +
                      (seasonx+anthro_z+food_z+temp_z|subject_ID), 
                    data = pd2, family = weibull(link="log",link_shape="log"))

mprior5$prior[2:12] <- "normal(0,1)"

make_stancode(strength ~ seasonx*syndrome_z + anthro_z*syndrome_z + 
                food_z*syndrome_z + temp_z*syndrome_z + rank_z + sex +
                (seasonx+anthro_z+food_z+temp_z|subject_ID), 
              data = pd2, family = weibull(link="log",link_shape="log"), prior = mprior5)

mod3 = brm(strength ~ seasonx*syndrome_z + anthro_z*syndrome_z + 
             food_z*syndrome_z + temp_z*syndrome_z + rank_z + sex +
             (seasonx+anthro_z+food_z+temp_z|subject_ID), 
           data = pd2, family = weibull(link="log",link_shape="log"), prior = mprior5, 
           chains = 3, cores = 3, iter=5000, warmup = 2000, thin = 2, 
           control = list(adapt_delta = 0.95), sample_prior = "yes")

summary(mod3)
plot(mod3)

round(posterior_summary(mod3, probs = c(0.05, 0.95)),3)

# all diagnostics ok - remove 'nonsig' interactions

mprior6 = get_prior(strength ~ temp_z + syndrome_z + 
                      seasonx + anthro_z + food_z + rank_z + sex +
                      (seasonx+anthro_z+food_z+temp_z|subject_ID), 
                    data = pd2, family = weibull(link="log",link_shape="log"))

mprior6$prior[2:8] <- "normal(0,1)"

make_stancode(strength ~ temp_z + syndrome_z + 
                seasonx + anthro_z + food_z + rank_z + sex +
                (seasonx+anthro_z+food_z+temp_z|subject_ID), 
              data = pd2, family = weibull(link="log",link_shape="log"), prior = mprior6)

red_mod3 = brm(strength ~ temp_z + syndrome_z + 
                 seasonx + anthro_z + food_z + rank_z + sex +
                 (seasonx+anthro_z+food_z+temp_z|subject_ID), 
               data = pd2, family = weibull(link="log",link_shape="log"), 
               chains = 3, cores = 3, iter=5000, warmup = 2000, thin = 2, 
               control = list(adapt_delta = 0.975), sample_prior = "yes")

summary(red_mod3)
plot(red_mod3)

round(posterior_summary(red_mod3, probs = c(0.05, 0.95)),3)

# 25 PPC of strength model ####

#load("C:/Users/Megaport/OneDrive/Project_data/social_comp/analysis/results.Rdata")

pp_check(red_mod3)

summary(red_mod3)
# so we need a ppc for relationship between:
# seasonality
# sex

pred_vals3 = predict(red_mod3, summary = FALSE)
predictions3 = as.data.frame(t(pred_vals3))

mean = rowMeans(predictions3)

fit_mod3 <- 
  fitted(red_mod3) %>%
  as_tibble() %>%
  bind_cols(pd2)

fit_mod3$predicted = as.numeric(mean)

# create new plotting object with 100 draws from posterior

# now proper plot by randomly selecting 100 draws from the posterior

xx = floor(runif(100, 1, 3000))

# extract those from our predictions object

zz = predictions3[,xx]

colnames(zz) = c("pred1", "pred2", "pred3", "pred4", "pred5", "pred6", "pred7", "pred8", "pred9", "pred10",
                 "pred11", "pred12", "pred13", "pred14", "pred15", "pred16", "pred17", "pred18", "pred19", "pred20",
                 "pred21", "pred22", "pred23", "pred24", "pred25", "pred26", "pred27", "pred28", "pred29", "pred30",
                 "pred31", "pred32", "pred33", "pred34", "pred35", "pred36", "pred37", "pred38", "pred39", "pred40",
                 "pred41", "pred42", "pred43", "pred44", "pred45", "pred46", "pred47", "pred48", "pred49", "pred50",
                 "pred51", "pred52", "pred53", "pred54", "pred55", "pred56", "pred57", "pred58", "pred59", "pred60",
                 "pred61", "pred62", "pred63", "pred64", "pred65", "pred66", "pred67", "pred68", "pred69", "pred70",
                 "pred71", "pred72", "pred73", "pred74", "pred75", "pred76", "pred77", "pred78", "pred79", "pred80",
                 "pred81", "pred82", "pred83", "pred84", "pred85", "pred86", "pred87", "pred88", "pred89", "pred90",
                 "pred91", "pred92", "pred93", "pred94", "pred95", "pred96", "pred97", "pred98", "pred99", "pred100")

# create new plotting object

names(fit_mod3)
plot_data = fit_mod3[, c(1,3,4,6,10,11,23)]
names(plot_data)
colnames(plot_data) = c("Estimate", "lower", "upper", "observed", "sex", 
                        "season","mean_prediction")
plot_data = cbind(plot_data, zz)
plot_data$season2 = ifelse(plot_data$season == 1, "Mating", "Non-mating")


# season first

sort(unique(plot_data$season2))

mean_mating = mean(plot_data$Estimate[plot_data$season2 == "Mating"])
sd_mating = sd(plot_data$Estimate[plot_data$season2 == "Mating"])
lower_mating = mean(plot_data$lower[plot_data$season2 == "Mating"])
upper_mating = mean(plot_data$upper[plot_data$season2 == "Mating"])

mean_non = mean(plot_data$Estimate[plot_data$season2 == "Non-mating"])
sd_non = sd(plot_data$Estimate[plot_data$season2 == "Non-mating"])
lower_non = mean(plot_data$lower[plot_data$season2 == "Non-mating"])
upper_non = mean(plot_data$upper[plot_data$season2 == "Non-mating"])

seasons = c("Mating", "Non-mating")
means2 = c(mean_mating, mean_non)
lowers2 = c(lower_mating, lower_non)
uppers2 = c(upper_mating, upper_non)
sds2 = c(sd_mating,sd_non)

plot_data3 = cbind(means2,lowers2,uppers2,sds2)
plot_data3 = as.data.frame(plot_data3)
plot_data3$season = seasons

season_ppc = ggplot(plot_data)+
  geom_jitter(aes(x=season2, y=observed+0.025),height=0.05, alpha=1/10)+
  geom_jitter(aes(x=season2, y=mean_prediction), height=0.05, alpha=1/10, colour = "red")+
  geom_errorbar(data=plot_data3, aes(x=season,ymin=means2-sds2, ymax=means2+sds2), width=0.5)+
  geom_point(data=plot_data3, aes(x=season,y=means2))+
  labs(x = "Season", y = "Network strength")
season_ppc = season_ppc+theme_classic()
#season_ppc = season_ppc+theme(legend.position="none")
season_ppc = season_ppc+ggtitle("(a)")
season_ppc

# sex next

mean_male = mean(plot_data$Estimate[plot_data$sex == "M"])
sd_male = sd(plot_data$Estimate[plot_data$sex == "M"])
lower_male = mean(plot_data$lower[plot_data$sex == "M"])
upper_male = mean(plot_data$upper[plot_data$sex == "M"])
mean_fem = mean(plot_data$Estimate[plot_data$sex == "F"])
sd_fem = sd(plot_data$Estimate[plot_data$sex == "F"])
lower_fem = mean(plot_data$lower[plot_data$sex == "F"])
upper_fem = mean(plot_data$upper[plot_data$sex == "F"])

sex = c("M", "F")
means2 = c(mean_male, mean_fem)
lowers2 = c(lower_male, lower_fem)
uppers2 = c(upper_male, upper_fem)
sds2 = c(sd_male,sd_fem)

plot_data3 = cbind(means2,lowers2,uppers2,sds2)
plot_data3 = as.data.frame(plot_data3)
plot_data3$sex = sex

sex_ppc = ggplot(plot_data)+
  geom_jitter(aes(x=sex, y=observed+0.025),height=0.05, alpha=1/10)+
  geom_jitter(aes(x=sex, y=mean_prediction), height=0.05, alpha=1/10, colour = "red")+
  geom_errorbar(data=plot_data3, aes(x=sex,ymin=means2-sds2, ymax=means2+sds2), width=0.5)+
  geom_point(data=plot_data3, aes(x=sex,y=means2))+
  labs(x = "Sex", y = "Network strength")
sex_ppc = sex_ppc+theme_classic()
#season_ppc = season_ppc+theme(legend.position="none")
sex_ppc = sex_ppc+ggtitle("(b)")
sex_ppc


# 26 Model Estimates for Table 2 ####

round(posterior_summary(red_mod3, probs = c(0.05, 0.95)),3)

m3.post <- posterior_samples(red_mod3)
names(m3.post)

sum(m3.post$b_syndrome_z > 0) / length(m3.post$b_syndrome_z)
sum(m3.post$b_seasonx > 0) / length(m3.post$b_seasonx)
sum(m3.post$b_anthro_z < 0) / length(m3.post$b_anthro_z)
sum(m3.post$b_food_z < 0) / length(m3.post$b_food_z)
sum(m3.post$b_temp_z < 0) / length(m3.post$b_temp_z)
sum(m3.post$b_rank_z > 0) / length(m3.post$b_rank_z)
sum(m3.post$b_sexM < 0) / length(m3.post$b_sexM)
sum(m3.post$'b_seasonx:syndrome_z' > 0) / length(m3.post$'b_seasonx:syndrome_z')
sum(m3.post$'b_syndrome_z:temp_z' < 0) / length(m3.post$'b_syndrome_z:temp_z')

# and for the supplementary materials

m3.postb <- posterior_samples(mod3)
names(m3.postb)

sum(m3.postb$b_syndrome_z > 0) / length(m3.postb$b_syndrome_z)
sum(m3.postb$b_seasonx > 0) / length(m3.postb$b_seasonx)
sum(m3.postb$b_anthro_z < 0) / length(m3.postb$b_anthro_z)
sum(m3.postb$b_food_z < 0) / length(m3.postb$b_food_z)
sum(m3.postb$b_temp_z < 0) / length(m3.postb$b_temp_z)
sum(m3.postb$b_rank_z > 0) / length(m3.postb$b_rank_z)
sum(m3.postb$b_sexM < 0) / length(m3.postb$b_sexM)
sum(m3.postb$'b_seasonx:syndrome_z' > 0) / length(m3.postb$'b_seasonx:syndrome_z')
sum(m3.postb$'b_syndrome_z:anthro_z' > 0) / length(m3.postb$'b_syndrome_z:anthro_z')
sum(m3.postb$'b_syndrome_z:food_z' > 0) / length(m3.postb$'b_syndrome_z:food_z')
sum(m3.postb$'b_syndrome_z:temp_z' < 0) / length(m3.postb$'b_syndrome_z:temp_z')

# 27 Plot Strength Results ####

# season first

sort(unique(plot_data$season2))

mean_mating = mean(plot_data$Estimate[plot_data$season2 == "Mating"])
sd_mating = sd(plot_data$Estimate[plot_data$season2 == "Mating"])
lower_mating = mean(plot_data$lower[plot_data$season2 == "Mating"])
upper_mating = mean(plot_data$upper[plot_data$season2 == "Mating"])

mean_non = mean(plot_data$Estimate[plot_data$season2 == "Non-mating"])
sd_non = sd(plot_data$Estimate[plot_data$season2 == "Non-mating"])
lower_non = mean(plot_data$lower[plot_data$season2 == "Non-mating"])
upper_non = mean(plot_data$upper[plot_data$season2 == "Non-mating"])

seasons = c("Mating", "Non-mating")
means2 = c(mean_mating, mean_non)
lowers2 = c(lower_mating, lower_non)
uppers2 = c(upper_mating, upper_non)
sds2 = c(sd_mating,sd_non)

plot_data2 = cbind(means2,lowers2,uppers2,sds2)
plot_data2 = as.data.frame(plot_data2)
plot_data2$season = seasons

str(plot_data)
str(plot_data2)

plot_data2$season2 = ifelse(plot_data2$season == "Mating", 1, 0)

fig3a = ggplot(plot_data)+
  geom_line(colour = "black", aes(x=season, y=observed),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred2),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred3),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred4),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred5),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred6),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred7),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred8),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred9),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred10),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred11),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred12),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred13),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred14),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred15),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred16),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred17),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred18),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred19),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred20),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred21),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred22),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred23),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred24),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred25),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred26),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred27),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred28),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred29),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred30),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred31),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred32),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred33),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred34),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred35),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred36),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred37),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred38),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred39),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred40),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred41),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred42),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred43),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred44),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred45),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred46),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred47),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred48),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred49),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred50),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred51),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred52),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred53),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred54),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred55),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred56),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred57),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred58),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred59),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred60),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred61),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred63),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred64),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred65),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred66),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred67),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred68),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred69),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred70),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred71),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred72),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred73),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred74),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred75),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred76),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred77),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred78),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred79),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred80),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred81),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred82),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred83),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred84),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred85),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred86),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred87),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred88),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred89),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred90),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred91),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred92),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred93),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred94),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred95),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred96),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred97),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred98),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred99),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred100),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_errorbar(data=plot_data2, aes(x=season2,ymin=means2-sds2, ymax=means2+sds2), width=0.5)
fig3a = fig3a + coord_cartesian(ylim = c(0, 20))
fig3a = fig3a+theme_classic()
fig3a = fig3a + labs(x = "",
                     y = "Network strength")
#
fig3a = fig3a+ggtitle("(a)")   
fig3a  


# sex next

mean_male = mean(plot_data$Estimate[plot_data$sex == "M"])
sd_male = sd(plot_data$Estimate[plot_data$sex == "M"])
lower_male = mean(plot_data$lower[plot_data$sex == "M"])
upper_male = mean(plot_data$upper[plot_data$sex == "M"])
mean_fem = mean(plot_data$Estimate[plot_data$sex == "F"])
sd_fem = sd(plot_data$Estimate[plot_data$sex == "F"])
lower_fem = mean(plot_data$lower[plot_data$sex == "F"])
upper_fem = mean(plot_data$upper[plot_data$sex == "F"])

sex = c("M", "F")
means2 = c(mean_male, mean_fem)
lowers2 = c(lower_male, lower_fem)
uppers2 = c(upper_male, upper_fem)
sds2 = c(sd_male,sd_fem)

plot_data3 = cbind(means2,lowers2,uppers2,sds2)
plot_data3 = as.data.frame(plot_data3)
plot_data3$sex = sex

sex_ppc = ggplot(plot_data)+
  geom_jitter(aes(x=sex, y=observed+0.025),height=0.05, alpha=1/10)+
  geom_jitter(aes(x=sex, y=mean_prediction), height=0.05, alpha=1/10, colour = "red")+
  geom_errorbar(data=plot_data3, aes(x=sex,ymin=means2-sds2, ymax=means2+sds2), width=0.5)+
  geom_point(data=plot_data3, aes(x=sex,y=means2))+
  labs(x = "Sex", y = "Network strength")
sex_ppc = sex_ppc+theme_classic()
#season_ppc = season_ppc+theme(legend.position="none")
sex_ppc = sex_ppc+ggtitle("(b)")
sex_ppc


# figure 4

mean_season1 = mean(bold$Estimate[bold$season2 == "Mating"])
sd_season1 = sd(bold$Estimate[bold$season2 == "Mating"])
lower_season1 = mean(bold$lower[bold$season2 == "Mating"])
upper_season1 = mean(bold$upper[bold$season2 == "Mating"])
mean_season2 = mean(bold$Estimate[bold$season2 == "Non-mating"])
sd_season2 = sd(bold$Estimate[bold$season2 == "Non-mating"])
lower_season2 = mean(bold$lower[bold$season2 == "Non-mating"])
upper_season2 = mean(bold$upper[bold$season2 == "Non-mating"])

season = c("Mating", "Non-mating")
means = c(mean_season1, mean_season2)
lowers = c(lower_season1, lower_season2)
uppers = c(upper_season1, upper_season2)
sds = c(sd_season1,sd_season2)

bold2 = cbind(means,lowers,uppers,sds)
bold2 = as.data.frame(bold2)
bold2$season = season

mean_season1 = mean(shy$Estimate[shy$season2 == "Mating"])
sd_season1 = sd(shy$Estimate[shy$season2 == "Mating"])
lower_season1 = mean(shy$lower[shy$season2 == "Mating"])
upper_season1 = mean(shy$upper[shy$season2 == "Mating"])
mean_season2 = mean(shy$Estimate[shy$season2 == "Non-mating"])
sd_season2 = sd(shy$Estimate[shy$season2 == "Non-mating"])
lower_season2 = mean(shy$lower[shy$season2 == "Non-mating"])
upper_season2 = mean(shy$upper[shy$season2 == "Non-mating"])

season = c("Mating", "Non-mating")
means = c(mean_season1, mean_season2)
lowers = c(lower_season1, lower_season2)
uppers = c(upper_season1, upper_season2)
sds = c(sd_season1,sd_season2)

shy2 = cbind(means,lowers,uppers,sds)
shy2 = as.data.frame(shy2)
shy2$season = season


geom_errorbar(data=bold2, aes(x=season,ymin=means-sds, ymax=means+sds), width=0.5)+
  geom_point(data=bold2, aes(x=season,y=means))+
  labs(x = "", y = "Network strength (log transformed)")
bseason_ppc = bseason_ppc + coord_cartesian(ylim = c(0, 1.2))
bseason_ppc = bseason_ppc+theme_classic()
#season_ppc = season_ppc+theme(legend.position="none")
bseason_ppc = bseason_ppc+ggtitle("Bold")



str(bold)
fig4a = ggplot(bold)+
  geom_line(colour = "black", aes(x=season, y=observed),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred2),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred3),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred4),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred5),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred1),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred6),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred7),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred8),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred9),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred10),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred11),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred12),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred13),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred14),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred15),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred16),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred17),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred18),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred19),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred20),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred21),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred22),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred23),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred24),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred25),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred26),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred27),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred28),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred29),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred30),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred31),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred32),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred33),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred34),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred35),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred36),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred37),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred38),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred39),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred40),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred41),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred42),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred43),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred44),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred45),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred46),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred47),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred48),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred49),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred50),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred51),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred52),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred53),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred54),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred55),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred56),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred57),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred58),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred59),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred60),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred61),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred62),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred63),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred64),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred65),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred66),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred67),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred68),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred69),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred70),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred71),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred72),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred73),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred74),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred75),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred76),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred77),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred78),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred79),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred80),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred81),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred82),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred83),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred84),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred85),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred86),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred87),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred88),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred89),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred90),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred91),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred92),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred93),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred94),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred95),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred96),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred97),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred98),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred99),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_line(colour = "#A40122", aes(x=season, y=pred100),
            stat = "smooth", method = "glm", method.args = list(family = gaussian(link='identity')),
            se = FALSE, size=1, alpha=0.1)+
  geom_errorbar(data=bold2, aes(x=season2,ymin=means-sds, ymax=means+sds), width=0.5)+
  geom_point(data=bold2, aes(x=season2,y=means))
fig4a = fig4a + coord_cartesian(ylim = c(0, 0.85))
fig4a = fig4a+theme_classic()
fig4a = fig4a + labs(x = "",
                     y = "Network strength (log transformed)")
#
fig4a = fig4a+ggtitle("Bold")   
fig4a  



# 28 ESS tables for supplementary ####

effective_sample(red_mod1, effects = "all")
effective_sample(red_mod2, effects = "all")
effective_sample(red_mod3, effects = "all")

# 31 Save ####

save.image("results.Rdata")

