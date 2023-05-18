#### Canada Goose Single Choice Analysis 

rm(list = ls ()) # clears R's memory

library(lme4)
library(tidyverse)
library(afex)
library(psych)
library(ggplot2)
library(stargazer)
library(ggeffects)
library(sjPlot)
library(emmeans)
library(outliers)
library(Rmisc)
library(sjstats)
library(car)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(PCAmixdata)
library(fastICA)
library(easystats)
library(MuMIn)
library(merTools)
library(interactions)
library(ggpubr)
library(jtools)
library(pavo)
library(latex2exp)
library(wavecolour)
getwd()
setwd("/Users/Ryan/Desktop/Single_Choice_MS_Final")

###########################
####Loading Input Files####

##Load object of interest/background reflectance/radiance info
OBJ<-read.csv("Clear.12pm_LEDs.csv", header=TRUE)
OBJECTS<-as.rspec(OBJ)

OBJ_cloud<-read.csv("Cloudy.12pm_LEDs.csv", header=TRUE)
OBJECTS_cloud<-as.rspec(OBJ_cloud)
OBJECTS_cloud
head(OBJ_cloud)

OBJECTS_cloud<-cbind(OBJECTS_cloud,OBJECTS[2:203])


##Load ambient light/irradiance spectrum file
irr_clear<-read.csv("Clear.12pm Irradiance.csv", header=TRUE)# opening the ambient light/irradiance spectrum file
irradiance_clear<-as.rspec(irr_clear)  # converts to rspec object for pavo
head(irr_clear)

irr_cloud<-read.csv("Cloudy.12pm Irradiance.csv", header=TRUE)# opening the ambient light/irradiance spectrum file
irradiance_cloud<-as.rspec(irr_cloud)  # converts to rspec object for pavo
head(irr_cloud)

##Ocular Media Transmittance curve
cago_om<-read.csv("Ocular_Media_CAGO.csv", header=T)
head(cago_om)
cago_OMT<-as.vector(cago_om$x)

##Sensitivity Curves Calculation CAGO Bird##
CAGO_CCpeaksense <- c(409, 458, 509, 580) #vector of lambda maxes for UV/VS, SWS, MWS, AND LWS
CAGO_CCLcut <- c(200, 459, 506, 559) #vector of lambda cuts for oil droplets T, C, Y, R
CAGO_CCBM <- c(1, 0.031, 0.042, 0.030) #vector of Bmids for T, C, Y, R
CAGO_CC_sens <- sensmodel(CAGO_CCpeaksense, lambdacut = CAGO_CCLcut, Bmid = CAGO_CCBM, om = cago_OMT, integrate=TRUE, range = c(300, 700), beta = TRUE)
plot(CAGO_CC_sens)


##Calculate Chromatic Contrast CAGO Bird    
CAGO.vis.mod_iter<- vismodel(OBJECTS, qcatch="Qi", visual=CAGO_CC_sens, illum= "ideal", relative=FALSE) #incorporating visual system information from a bird
summary(CAGO.vis.mod_iter) #a summary of what was performed in the calculation
#the order of relative densities used in the argument below is the following CAGO, SWS, MWS, LWS
CAGO.cc.col.dis<- coldist(CAGO.vis.mod_iter, noise="neural", achro=FALSE, n = c(1, 10.93, 12.46, 16.93), weber = 0.1)
CAGO.Bird.CC<-CAGO.cc.col.dis[1:201,2:3]
colnames(CAGO.Bird.CC)<-c("LEDs", "CAGO.CC")
CAGO.Bird.CC$wl<-seq(300,700,by=2)

ggplot(CAGO.Bird.CC, aes(x = wl, y = CAGO.CC, colour=wl)) + 
  geom_line(lwd=3)+
  scale_color_wavelength() + 
  ylab("Chromatic Contrast (JND)")+
  xlab("Wavelength (nm)")+
  theme_classic(base_size = 18)+
  theme(legend.position="none")
head(CAGO.Bird.CC)
str(CAGO.Bird.CC)

#Theoretical Contrast for Blue
CAGO.Bird.CC[which(CAGO.Bird.CC$wl==482),]$CAGO.CC
CAGO.Bird.CC[which(CAGO.Bird.CC$wl==484),]$CAGO.CC
##Theoretical Contrast for Red
CAGO.Bird.CC[which(CAGO.Bird.CC$wl==630),]$CAGO.CC
CAGO.Bird.CC[which(CAGO.Bird.CC$wl==632),]$CAGO.CC

##Calculate Chromatic Contrast CAGO Bird    
CAGO.vis.mod_iter_clear<- vismodel(OBJECTS, qcatch="Qi", visual=CAGO_CC_sens, illum= irradiance_clear$X1200_Clear, relative=FALSE) #incorporating visual system information from a bird
summary(CAGO.vis.mod_iter_clear) #a summary of what was performed in the calculation
#the order of relative densities used in the argument below is the following CAGO, SWS, MWS, LWS
CAGO.cc.col.dis_clear<- coldist(CAGO.vis.mod_iter_clear, noise="neural", achro=FALSE, n = c(1, 10.93, 12.46, 16.93), weber = 0.1)
CAGO.Bird.CC_clear<-CAGO.cc.col.dis_clear[1:201,2:3]
colnames(CAGO.Bird.CC_clear)<-c("LEDs", "CAGO.CC")
CAGO.Bird.CC_clear$wl<-seq(300,700,by=2)

##Calculate Chromatic Contrast CAGO Bird    
CAGO.vis.mod_iter_cloud<- vismodel(OBJECTS_cloud, qcatch="Qi", visual=CAGO_CC_sens, illum= irradiance_cloud$X1200_Cloudy, relative=FALSE) #incorporating visual system information from a bird
summary(CAGO.vis.mod_iter_cloud) #a summary of what was performed in the calculation
#the order of relative densities used in the argument below is the following CAGO, SWS, MWS, LWS
CAGO.cc.col.dis_cloud<- coldist(CAGO.vis.mod_iter_cloud, noise="neural", achro=FALSE, n = c(1, 10.93, 12.46, 16.93), weber = 0.1)
CAGO.Bird.CC_cloud<-CAGO.cc.col.dis_cloud[1:201,2:3]
colnames(CAGO.Bird.CC_cloud)<-c("LEDs", "CAGO.CC")
CAGO.Bird.CC_cloud$wl<-seq(300,700,by=2)

#Sunny
#Theoretical Contrast for Blue
Blue_JND_Clear<-mean(CAGO.Bird.CC_clear[which(CAGO.Bird.CC_clear$wl==482),]$CAGO.CC,CAGO.Bird.CC_clear[which(CAGO.Bird.CC_clear$wl==484),]$CAGO.CC)

##Theoretical Contrast for Red
Red_JND_Clear<-mean(CAGO.Bird.CC_clear[which(CAGO.Bird.CC_clear$wl==630),]$CAGO.CC,CAGO.Bird.CC_clear[which(CAGO.Bird.CC_clear$wl==632),]$CAGO.CC)


#Cloudy
#Theoretical Contrast for Blue
Blue_JND_Cloud<-mean(CAGO.Bird.CC_cloud[which(CAGO.Bird.CC_cloud$wl==482),]$CAGO.CC,CAGO.Bird.CC_cloud[which(CAGO.Bird.CC_cloud$wl==484),]$CAGO.CC)

##Theoretical Contrast for Red
Red_JND_Cloud<-mean(CAGO.Bird.CC_cloud[which(CAGO.Bird.CC_cloud$wl==630),]$CAGO.CC,CAGO.Bird.CC_cloud[which(CAGO.Bird.CC_cloud$wl==632),]$CAGO.CC)


ggplot(irradiance_clear, aes(x = wl, y = X1200_Clear, colour=wl)) + 
  geom_line(lwd=2)+
  scale_color_wavelength() + 
  ylab(expression(mu~"mol/m"^2~"/s"))+
  xlab("Wavelength (nm)")+
  theme_classic(base_size = 24)+
  theme(legend.position="none")

head(irr_cloud)
ggplot(irradiance_cloud, aes(x = wl, y = X1200_Cloudy, colour=wl)) + 
  geom_line(lwd=2)+
  scale_color_wavelength() + 
  ylab(expression(mu~"mol/m"^2~"/s"))+
  ylim(0,3)+
  xlab("Wavelength (nm)")+
  theme_classic(base_size = 24)+
  theme(legend.position="none")

wl<-seq(300,700,by=1)
ggplot(data=NULL, aes(x = wl , y = OBJ$x524, colour=wl)) + 
  geom_line(lwd=2)+
  ylim(c(0,4000))+
  scale_color_wavelength() + 
  ylab("Photon(counts)")+
  xlab("Wavelength (nm)")+
  theme_classic(base_size = 24)+
  theme(legend.position="none")


#Spectra Intensity Selection

getwd()
setwd("/Users/Ryan/Desktop/Single_Choice_MS")
spec_data<-read.csv("Light_Spectra.csv")


#Blue Light
sum(spec_data$Blue_20cd)
sum(spec_data$Blue_40cd)
sum(spec_data$Blue_80cd)
sum(spec_data$Blue_120cd)

#Red Light
sum(spec_data$Red_40cd)
sum(spec_data$Red_80cd)
sum(spec_data$Red_120cd)
sum(spec_data$Red_240cd)


#Blue Light
max(spec_data$Blue_20cd)
max(spec_data$Blue_40cd)
spec_data[which.max(spec_data$Blue_80cd),]
max(spec_data$Blue_120cd)

#Red Light
max(spec_data$Red_40cd)
max(spec_data$Red_80cd)
max(spec_data$Red_120cd)
max(spec_data$Red_240cd)
spec_data[which.max(spec_data$Red_120cd),]

#Differences in photon emitted at peak wavelength
abs(max(spec_data$Red_40cd)-max(spec_data$Blue_20cd))
abs(max(spec_data$Red_80cd)-max(spec_data$Blue_40cd))
abs(max(spec_data$Red_120cd)-max(spec_data$Blue_80cd))
abs(max(spec_data$Red_240cd)-max(spec_data$Blue_120cd))


#Goose ID generation

a<-(replicate(n=23,sample(1:6, 4, replace=TRUE),simplify = T))
df<-t(a)
df
df1<-replace(df, df==1,"W")
df2<-replace(df1, df1==2,"Y")
df3<-replace(df2, df2==3,"G")
df4<-replace(df3, df3==4,"B")
df5<-replace(df4, df4==5,"R")
df_<-replace(df5, df5==6,"O")
df_

colnames(df_)<-c("one","two","three","four")
df<-data.frame(df_)
df
df$ID<-paste0(df$one,df$two,df$three,df$four, sep="")
df

write.csv(df, file="BIRD_ID.csv")

###Manipulative Experiment 
man_data<-read.csv("Side_Bias_Test_Data.csv")
head(man_data)
man_data$Notes<-NULL

#Separating the data into three distinct days 
man_data_1 <-man_data[which(man_data$Date=="12/6/20"),]
man_data_2<-man_data[which(man_data$Date=="12/9/20"),]
man_data_3<-man_data[which(man_data$Date=="12/11/20"),]

#Testing to see if the differences in left versus right choices within the arena were significantly different for each day of the trial 
binom.test(sum(man_data_1$Choice..Binary.),nrow(man_data_1),.5)
binom.test(sum(man_data_2$Choice..Binary.),nrow(man_data_2),.5)
binom.test(sum(man_data_3$Choice..Binary.),nrow(man_data_3),.5)

#Testing to see if the differences in left versus right choices within the arena were significantly different across all three days of the trial
binom.test(sum(man_data$Choice..Binary.),nrow(man_data),.5)

#Reading the dataset for the single choice test
bird_data <- read.csv("Light_Stimuli_Experiment_Data.csv", na.strings = c("","NA"), header=TRUE) #you can save a normal Excel file as a .csv using "Save As" in excel

# important for afex
set_sum_contrasts()

# summary of the types of variables and their values 
str(bird_data) 
head(bird_data)

str(bird_data)
# explicitly assigning different variables as factors
bird_data$bird_id <- as.factor(bird_data$bird_id) 
bird_data$color <- as.factor(bird_data$color) 
bird_data$frequency <- as.factor(bird_data$frequency) 
bird_data$light_position <- as.factor(bird_data$light_position) 
bird_data$trial_order <- as.factor(bird_data$trial_order)

# explicitly assigning binary as a numeric veriable
bird_data$preference_binary <- as.numeric(bird_data$preference_binary)

#checking changes were correctly executed
str(bird_data)

# checking whether there is no missing trial
table(bird_data$bird_id, bird_data$trial_order, bird_data$color)

# changing order of levels within frequency treatment
bird_data$frequency <- relevel(bird_data$frequency, "Steady") 

# turning light position into a binary variable 
bird_data$light_position_binary <-ifelse(bird_data$light_position=="Left",0,1) 

# checking the correlation between light at position and light NOT at position
# grouping all independent continuous factors into a single object to facilitate running the correlations
ind.cont <- bird_data[c("light_int_at_position", "light_int_not_at_position")]

# runs the pairwise correlations between the independent continuous factors. 
corr.test(ind.cont, use = "pairwise", method = "pearson", adjust = "none") 
pairs(~bird_data$light_int_at_position + bird_data$light_int_not_at_position, pch = 19, lower.panel = NULL)


#### creating the dependent variables of head movement and body movement rate by dividing head and body movement frequency by video time. 
# estimating head movement rate as an index of scanning 
bird_data <- mutate(bird_data, head_move_rate = (head_move_freq/video_new))

# estimating body movement rate as an index of spatial displacement
bird_data <- mutate(bird_data, body_move_rate = (body_move_freq/video_new))


# checking the correlation between head and body movement rate
# grouping all independent continuous factors into a single object to facilitate running the correlations
ind.cont <- bird_data[c("head_move_rate", "body_move_rate")] 
# runs the pairwise correlations between the independent continuous factors. 
corr.test(ind.cont, use = "pairwise", method = "pearson", adjust = "none") 
pairs(~bird_data$head_move_rate + bird_data$body_move_rate, pch = 19, lower.panel = NULL)

####Summarizing ambient light intensity on both sides of the experimental arena with a PCA 

#Creating a data frame for light intensity at the position of the light source and not at the position of the light 
data.pca1 <- data.frame(bird_data$light_int_at_position, bird_data$light_int_not_at_position)

#######Running PCA
res.pca1 <- PCA (data.pca1, graph = T)

#The eigenvalues for the PCA with light intensity at and not at position
eig.val1 <- get_eigenvalue((res.pca1))
eig.val1
#Only PCA 1 has an eigenvalue above 1

#ploting dimensions 1 and 2 for PCA
fviz_pca_var (res.pca1, col.var = "blue")

#extracting the coordinates for light intensity at position for PCA 
ind1 <- get_pca_ind(res.pca1)

#Putting the individual coordinates into a dataframe
df<-data.frame(ind1$coord)
colnames(df)<-c("PCA_1_light_int","PCA_2_light_int")

#numbering the rows for each coordinate that correspond to the rows of the original dataset
df<-cbind(row=rownames(df),df)

#creating and numbering the rows in the orginal dataset
bird_data<-cbind(row=rownames(bird_data),bird_data)

#merging the rows based on row name so the coordinates match the corresponding row of each observation in the orginal dataset
bird_data<-merge(bird_data,df, by="row")

# We will use in the rest of the analysis PCA_1_light_int as an indicator of light intensity in the enclosure

### PCA ends

### Checking associations of confounding factors with PCA factor
# grouping all independent continuous factors into a single object to facilitate running the correlations
ind.cont <- bird_data[c("PCA_1_light_int", "time_decimal", "temp")] 
# runs the pairwise correlations between the independent continuous factors. 
corr.test(ind.cont, use = "pairwise", method = "pearson", adjust = "none") 
pairs(~bird_data$PCA_1_light_int + bird_data$time_decimal + bird_data$temp, pch = 19, lower.panel = NULL)

# Given the significant associations between PCA_1_light_int (PCA factor) and the other environmental factors, we chose to
# only include in the models the PCA factor. 

###### checking the association between light position vs. color:frequency:trail order

bird_data$light_position_binary <-ifelse(bird_data$light_position=="Right",1,0)

# model 
vscftrialorder <- glm(data = bird_data, light_position_binary ~ color*frequency*trial_order , family = binomial(link = "logit")) 

# summary of the model
summary(vscftrialorder)

# looking at whether each independent factor is significantly explaining the variation in the dependent vars 
anova(vscftrialorder, test="Chisq")

# Given the significant association generating a multicollinearity issue, we chose not to include light position as an independent factor

#####################################3 LATENCY Analysis

# Phase 1: Model with only the interaction between color and frequency 
latency_m1 <- mixed(latency_new ~ color + frequency + trial_order + PCA_1_light_int + color:frequency + (1|bird_id), data = bird_data, method = "KR",  
                    control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (latency_m1)

# removing non significant interactions: color:frequency


# Phase 2: Exploring interactions with just trial_order between color and frequency
latency_m2 <- mixed(latency_new ~ color + frequency + trial_order + PCA_1_light_int + color:trial_order + frequency:trial_order + (1|bird_id), data = bird_data, method = "KR",  
                    control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (latency_m2)

# removing non significant interactions: color:trial_order and frequency:trial_order 


#Phase 3: Exploring interactions with just light intensity between color and frequency
latency_m3 <- mixed(latency_new ~ color + frequency + trial_order + PCA_1_light_int + color:PCA_1_light_int + frequency:PCA_1_light_int + (1|bird_id), data = bird_data, method = "KR",  
                    control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (latency_m3)

# we remove non significant interactions: color:PCA_1_light_int and frequency:PCA_1_light_int 

# Phase 4: final model removing non-significant interactions in previous phases
latency_m4 <- mixed(latency_new ~ color + frequency + trial_order + PCA_1_light_int + (1|bird_id), data = bird_data, method = "KR",  
                    control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova(latency_m4)
summary(latency_m4)

# checking assumptions
# homogeneity of variances in two ways
# way 1:
plot(latency_m4$full_model)

#way 2:
boxplot(residuals(latency_m4$full_model) ~ bird_data$color + bird_data$frequency)

# normality of the residuals
qqnorm(residuals(latency_m4$full_model))

# there is some deviation of the residuals from homogeneity of variances and normality of residuals

# log-transforming latency_new_frames
bird_data <- mutate(bird_data, latency_new_log = log (latency_new))


# new model with log transformation 
latency_m4_t <- mixed(latency_new_log ~ color + frequency + trial_order + PCA_1_light_int + (1|bird_id), data = bird_data,  method = "KR",
                      control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (latency_m4_t)

# checking assumptions again 

# way 1:
plot(latency_m4_t$full_model)

#way 2:
boxplot(residuals(latency_m4_t$full_model) ~ bird_data$color + bird_data$frequency)

# normality of the residuals
qqnorm(residuals(latency_m4_t$full_model))

# a log transformation improved the fit of the model to the assumptions, however, to avoid interpretation issues with the 
# parameter estimates after a log transformation, we keep the untrasformed model

# interpreting fixed effects
anova (latency_m4)

# calculating means
emm_options(lmer.df = "kenward-roger") 
emm_color <- emmeans(latency_m4, "color", model = "multivariate")
emm_color

# calculating means
emm_options(lmer.df = "kenward-roger") 
emm_freq <- emmeans(latency_m4, "frequency", model = "multivariate")
emm_freq

# calculating means
emm_options(lmer.df = "kenward-roger") 
emm_contrast <- emmeans(latency_m4, "trial_order", model = "multivariate")
emm_contrast


# interpreting random effects
summary (latency_m4$full_model)
vcov(latency_m4$full_model)

#Estmating the marginal and conditional R-squared value
r.squaredGLMM(latency_m4$full_model)

# R2marginal (associated with fixed effects): 0.05298755
# R2conditional (associated with both fixed and random effects): 0.1897191

summary (latency_m4$full_model)

# repeatability
987.8 / (987.8 + 5854.4)


# 0.1443688 repeatability means that about 14% of the population behavioral variance is associated with
# between-individual differences in behavior

# confirming some results using other packages

# code for interpreting random effects 
report_performance(latency_m4$full_model)

# R2marginal (associated with fixed effects): 0.05
# R2conditional (associated with both fixed and random effects): 0.19
get_variance(latency_m4$full_model)

# bird_id var: 987.8194
# residual var: 5854.391

# repeatability
987.8194 / (987.8194 + 5854.391)


# 0.144382 repeatability means that about 14% of the population behavioral variance is associated with
# between-individual differences in behavior


# estimating the uncertainty in repeatability for the latency
set.seed(1)

simulated <- sim (latency_m4$full_model, n.sim = 1000)

posterior_bird_id <- apply(simulated@ranef$"bird_id"[ , , 1],1,var)
posterior_residual <- simulated@sigma^2

quantile(posterior_bird_id /
           (posterior_bird_id + posterior_residual),
         prob=c(0.025, 0.5, 0.975))

# After controlling for the fixed factors, about 14% of the remaining variance in latency
# is associated with differences between individuals, with a range between 9% and 21%
# this means that some Canada geese have really high latency vs others with really low latency 
# that is not accounted for any of the fixed factors included in the model. 


# We can characterize the latency responses of individual animals that show high vs. low values 
randomSims <- REsim(latency_m4$full_model, n.sims = 1000)

head(randomSims[randomSims$groupFctr=="bird_id",])

randomSims$groupID <- factor(randomSims$groupID,
                             levels = randomSims$groupID[order(randomSims$mean)])

randomSims$mean <- randomSims$mean +
  fixef(latency_m4$full_model)["(Intercept)"]

#Supplementary Figure

afex_plot(latency_m4_t, x = "trial_order", id = "bird_id", dodge = 0.4, point_arg = list(size = 4), factor_levels = list(trial_order = c("1th", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th"))) + labs(y = "Latency (s)", x = "Trial order") + theme_pubr(20) 


ggplot()+
  geom_errorbar(data = randomSims,
                aes(x = groupID, ymin = mean-sd,
                    ymax = mean+sd))+
  geom_point(data = randomSims,
             aes(x = groupID, y = mean), shape = 19, size = 3)+
  theme_classic(base_size = 20) + ylab("Latency to choose a side (sec)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())




###################### HEAD MOVEMENT RATE Analysis 

# Phase 1: Model with only the interaction between color and frequency 
hmr_m1 <- mixed(head_move_rate ~ color + frequency + trial_order + PCA_1_light_int + color:frequency + (1|bird_id), data = bird_data, method = "KR",  
                control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (hmr_m1)

# removing color:frequency 

# Phase 2: Exploring interactions with just trial_order between color and frequency
hmr_m2 <- mixed(head_move_rate ~ color + frequency + trial_order + PCA_1_light_int + color:trial_order + frequency:trial_order + (1|bird_id), data = bird_data, method = "KR",  
                control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (hmr_m2)
# removing non-significant interactions: color:trial_order and frequency:trial_order

#Phase 3: Exploring interactions with just light intensity between color and frequency
hmr_m3 <- mixed(head_move_rate ~ color + frequency + trial_order + PCA_1_light_int + color:PCA_1_light_int + frequency:PCA_1_light_int + (1|bird_id), data = bird_data, method = "KR",  
                control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (hmr_m3)
# we remove non significant interactions: color:trial_order and frequency:color:PCA_1_light_int and frequency:PCA_1_light_int 


# Phase 4: final model removing non-significant interactions in previous phases
hmr_m4 <- mixed(head_move_rate ~ color + frequency + trial_order + PCA_1_light_int + (1|bird_id), data = bird_data, method = "KR",  
                control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (hmr_m4)


# checking assumptions
# homogeneity of variances in two ways
# way 1:
plot(hmr_m4$full_model)

#way 2:
boxplot(residuals(hmr_m4$full_model) ~ bird_data$color + bird_data$frequency)

# normality of the hmr_m4
qqnorm(residuals(hmr_m4$full_model))

# homogeneity of variances and normality of residuals appear to be met

# fixed effects are not significant but we get the means of color and frequency for reporting purposes
# calculating means
emm_options(lmer.df = "kenward-roger") 
emm_color <- emmeans(hmr_m4$full_model, "color", model = "multivariate")
emm_color

# calculating means
emm_options(lmer.df = "kenward-roger") 
emm_freq <- emmeans(hmr_m4$full_model, "frequency", model = "multivariate")
emm_freq

# interpreting random effects
summary (hmr_m4$full_model)
vcov(hmr_m4$full_model)
r.squaredGLMM(hmr_m4$full_model)

# R2marginal (associated with fixed effects): 0.04587605 
# R2conditional (associated with both fixed and random effects): 0.1730945


summary (hmr_m4$full_model)
# repeatability
0.04946 / (0.04946 + 0.32148)

# 0.1333369 repeatability means that about 13% of the population behavioral variance is associated with
# between-individual differences in behavior

# estimating the uncertainty in repeatability
set.seed(1)

simulated <- sim (hmr_m4$full_model, n.sim = 1000)

posterior_bird_id <- apply(simulated@ranef$"bird_id"[ , , 1],1,var)
posterior_residual <- simulated@sigma^2

quantile(posterior_bird_id /
           (posterior_bird_id + posterior_residual),
         prob=c(0.025, 0.5, 0.975))

# After controlling for the fixed factors, about 13% of the remaining variance in head movement rate
# is associated with differences between individuals, with a range between 7% and 20%
# this means that some Canada geese have really relatively low variation 
# that is not accounted for any of the fixed factors included in the model. 


# We can characterize the latency responses of individual animals that show high vs. low values 
randomSims <- REsim(hmr_m4$full_model, n.sims = 1000)

head(randomSims[randomSims$groupFctr=="bird_id",])

randomSims$groupID <- factor(randomSims$groupID,
                             levels = randomSims$groupID[order(randomSims$mean)])

randomSims$mean <- randomSims$mean +
  fixef(hmr_m4$full_model)["(Intercept)"]

#Supplementary Figure
afex_plot(hmr_m4, x = "trial_order", id = "bird_id", dodge = 0.4, point_arg = list(size = 4), factor_levels = list(trial_order = c("1th", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th"))) + labs(y = "Head Movement Rate (events/second)", x = "Trial order") + theme_pubr(20) 


ggplot()+
  geom_errorbar(data = randomSims,
                aes(x = groupID, ymin = mean-sd,
                    ymax = mean+sd))+
  geom_point(data = randomSims,
             aes(x = groupID, y = mean), shape = 19, size = 3)+
  theme_classic(base_size = 20) + ylab("Head movement rate (events/sec)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())


###### BODY MOVEMENT RATE Analysis

# Phase 1: Model with only the interaction between color and frequency
bmr_m1 <- mixed(body_move_rate ~ color + frequency + trial_order + PCA_1_light_int + color:frequency + (1|bird_id), data = bird_data, method = "KR",  
                control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (bmr_m1)

# removing the non-significant interaction: color:frequency

# Phase 2: Exploring interactions with just trial_order between color and frequency
bmr_m2 <- mixed(body_move_rate ~ color + frequency + trial_order + PCA_1_light_int + color:trial_order + frequency:trial_order + (1|bird_id), data = bird_data, method = "KR",  
                control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (bmr_m2)

# removing non-significant interactions: color:trial_order and frequency:trial_order

#Phase 3: Exploring interactions with just light intensity between color and frequency
bmr_m3 <- mixed(body_move_rate ~ color + frequency + trial_order + PCA_1_light_int + color:PCA_1_light_int  + frequency:PCA_1_light_int  + (1|bird_id), data = bird_data, method = "KR",  
                control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (bmr_m3)

# we remove non significant interactions: color:PCA_1_light_int  and frequency:PCA_1_light_int


# Phase 4: final model removing non-significant interactions in previous phases
bmr_m4 <- mixed(body_move_rate ~ color + frequency + trial_order + PCA_1_light_int + (1|bird_id), data = bird_data, method = "KR",  
                control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (bmr_m4)


# checking assumptions
# homogeneity of variances in two ways
# way 1:
plot(bmr_m4$full_model)

#way 2:
boxplot(residuals(bmr_m4$full_model) ~ bird_data$color + bird_data$frequency)

# normality of the residuals
qqnorm(residuals(bmr_m4$full_model))

# normality of residuals and homogeneity of variances seem to be met
# interpreting fixed effects
anova (bmr_m4)
summary(bmr_m4)

# calculating means
emm_options(lmer.df = "kenward-roger") 
emm_color <- emmeans(bmr_m4, "color", model = "multivariate")
emm_color


# calculating means
emm_options(lmer.df = "kenward-roger") 
emm_freq <- emmeans(bmr_m4, "frequency", model = "multivariate")
emm_freq

# calculating means
emm_options(lmer.df = "kenward-roger") 
emm_contrast <- emmeans(bmr_m4, "trial_order", model = "multivariate")
emm_contrast

pairs(emm_contrast)

# the pace of body movements does not seem to differ between color and frequency treatments,
# but it changes across increasing exposure to the treatments

# trail order figure
afex_plot(bmr_m4, x = "trial_order", id = "bird_id", dodge = 0.4, point_arg = list(size = 4), factor_levels = list(trial_order = c("1th", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th"))) + labs(y = "Body Movement Rate (events/second)", x = "Trial order") + theme_pubr(20) 


# interpreting random effects

summary (bmr_m4$full_model)

vcov(bmr_m4$full_model)

r.squaredGLMM(bmr_m4$full_model)

# R2marginal (associated with fixed effects): 0.1083641 
# R2conditional (associated with both fixed and random effects): 0.2642684

summary (bmr_m4$full_model)

# repeatability



0.02973 / (0.02973 + 0.14032)

# 0.1748309 repeatability means that about 17% of the population behavioral variance is associated with
# between-individual differences in behavior


# estimating the uncertainty in repeatability
set.seed(1)

simulated <- sim (bmr_m4$full_model, n.sim = 1000)

posterior_bird_id <- apply(simulated@ranef$"bird_id"[ , , 1],1,var)
posterior_residual <- simulated@sigma^2

quantile(posterior_bird_id /
           (posterior_bird_id + posterior_residual),
         prob=c(0.025, 0.5, 0.975))

# After controlling for the fixed factors, about 17% of the remaining variance in body movement rate
# is associated with differences between individuals, with a range between 10% and 25%
# this means that Canada geese have really relatively low variation in body movement
# that is not accounted for any of the fixed factors included in the model. 


# We can characterize the latency responses of individual animals that show high vs. low values 
randomSims <- REsim(bmr_m4$full_model, n.sims = 1000)

head(randomSims[randomSims$groupFctr=="bird_id",])

randomSims$groupID <- factor(randomSims$groupID,
                             levels = randomSims$groupID[order(randomSims$mean)])

randomSims$mean <- randomSims$mean +
  fixef(bmr_m4$full_model)["(Intercept)"]

#Supplementary Figure
ggplot()+
  geom_errorbar(data = randomSims,
                aes(x = groupID, ymin = mean-sd,
                    ymax = mean+sd))+
  geom_point(data = randomSims,
             aes(x = groupID, y = mean), shape = 19, size = 3)+
  theme_classic(base_size=20) + ylab("Body movement rate (events/sec)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())


############### Probability of Avoidance Analysis 

# Phase 1: Model with only the interaction between color and frequency
avoid_m1 <- mixed(preference_binary ~ color + frequency + trial_order + PCA_1_light_int + color:frequency + (1|bird_id), data = bird_data, method = "LRT", family = binomial, 
                  control = glmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (avoid_m1)
# removing non-significant interactions: color:frequency

# Phase 2: Exploring interactions with just trial_order between color and frequency
avoid_m2 <- mixed(preference_binary ~ color + frequency + trial_order + PCA_1_light_int + color:trial_order + frequency:trial_order + (1|bird_id), data = bird_data, method = "LRT", family = binomial, 
                  control = glmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (avoid_m2)

# removing non-significant interaction of frequency:trial_order but retained the significant interaction between color and trial order


#Phase 3: Exploring interactions with just light intensity between color and frequency
avoid_m3 <- mixed(preference_binary ~ color + frequency + trial_order + PCA_1_light_int + color:trial_order + color:PCA_1_light_int + frequency:PCA_1_light_int + (1|bird_id), data = bird_data, method = "LRT", family = binomial, 
                  control = glmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (avoid_m3)

# we removed non-significant interaction of color:PCA_1_light_int but retained the significant interaciton between freqeuncy:PCA_1_light_int

# Phase 4: final model removing non-significant interactions in previous phases
head(bird_data)
avoid_m4 <- mixed(preference_binary ~ color + frequency + trial_order + PCA_1_light_int + color:trial_order + frequency:PCA_1_light_int +  (1|bird_id), data = bird_data, method = "LRT", family = binomial, 
                  control = glmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (avoid_m4)
summary(avoid_m4)


# Estimating the fixed effects of the model 
emmeans(avoid_m4, "color", type = "response")

#Plot between color and bird ID


afex_plot(avoid_m4, x = "color", id = "bird_id", dodge = 0.8,
                mapping = c("color"),
                point_arg = list(size = 4), 
                error_arg = list(linewidth = 1, width = 0))+
  guides(color=guide_legend("Color"))+
  labs(y = "Probability of Avoidance (%)",x="Color") +
  scale_color_manual(values=c("blue","red"))+
  scale_y_continuous(labels=function(x)x*100)+ theme_pubr(20) 


emmeans(avoid_m4, "frequency", type = "response")

emmeans(avoid_m4, "trial_order", type = "response")

#Plot of trial order, color, and bird_ID
afex_plot(avoid_m4, x = "trial_order", trace = "color",id = "bird_id",dodge = 0.6, mapping=list("color","fill"),data_arg = list(
  position = ggplot2::position_jitterdodge(
      jitter.width = 0.3, 
      jitter.height = .015, 
      dodge.width = 0.6),alpha=.001), point_arg = list(size = 4), factor_levels = list(trial_order = c("1", "2", "3", "4", "5", "6", "7", "8"), color = c("Blue", "Red")), legend_title = "Color") + labs(y = "Probability of Avoidance (%)", x = "Trial order")+
  scale_color_manual(values=c("blue","red"))+
  scale_y_continuous(labels=function(x)x*100)+
  theme_pubr(20) 


emm_int <- emmeans(avoid_m4, "color", by = c("trial_order"), type = "response")
emm_int
pairs(emm_int)


# plotting interaction effect between frequency and PCA_1_light_int
interact_plot(avoid_m4$full_model, pred = PCA_1_light_int, modx = frequency, interval = TRUE, int.width = 0.95,
              line.thickness = 1.5, x.label = "Ambient light intensity (PCA1)",
              y.label = "Probability of Avoidance (%)", modx.labels = c("Steady", "Pulsing"), legend.main = "Light Frequency") + 
  theme_apa() + 
  scale_y_continuous(labels=function(x)x*100)+
  theme (axis.title.y = element_text(size=20), axis.text.y = element_text(size=20), 
         axis.title.x = element_text(size=20), axis.text.x = element_text(size=20),
         legend.title = element_text(size=20), legend.text = element_text(size=20)) 

help("interact_plot")
# code for interpreting random effects 

report_performance(avoid_m4$full_model)

# R2marginal (associated with fixed effects): 0.84 
# R2conditional (associated with both fixed and random effects): 0.86

summary(avoid_m4$full_model)

get_variance(avoid_m4$full_model)

# bird_id var: 0.2003978
# residual var: 3.289868


# repeatability
0.0917/ (0.0917 + 3.289868)

# 0.0271176 repeatability means that about 3% of the population behavioral variance is associated with
# between-individual differences in behavior


# We can characterize the responses of individual animals that show high vs. low values, but the output will be in the log odds scale rather than prob scale
randomSims <- REsim(avoid_m4$full_model, n.sims = 1000)

head(randomSims[randomSims$groupFctr=="bird_id",])

randomSims$groupID <- factor(randomSims$groupID,
                             levels = randomSims$groupID[order(randomSims$mean)])



randomSims$mean <- randomSims$mean +
  fixef(avoid_m4$full_model)["(Intercept)"]

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

randomSims$prob<-logit2prob(randomSims$mean)
randomSims$prob_sd<-logit2prob(randomSims$sd)

ggplot()+
  geom_errorbar(data = randomSims,
                aes(x = groupID, ymin = mean-sd,
                    ymax = mean+sd))+
  geom_point(data = randomSims,
             aes(x = groupID, y = mean), shape = 19, size = 3)+
  theme_classic(base_size=20) + ylab("Odds of avoiding light (log odds scale)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())


ggplot()+
  geom_errorbar(data = randomSims,
                aes(x = groupID, ymin = prob-prob_sd,
                    ymax = prob+prob_sd))+
  geom_point(data = randomSims,
             aes(x = groupID, y = prob), shape = 19, size = 3)+
  theme_classic() + ylab("Probability of avoiding light ()")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

# including random slopes for trial_order and comparing them with random intercepts
# this is to see the between individual variation per trial 

avoid_m4s <- mixed(preference_binary ~ color + frequency + light_position + trial_order + PCA_1_light_int + color:frequency + frequency:PCA_1_light_int + color:trial_order + (trial_order||bird_id), data = bird_data, method = "LRT", family = binomial, 
                   control = glmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

anova (avoid_m4s)

AIC (avoid_m4$full_model, avoid_m4s$full_model)


summary (avoid_m4s)

RI<-augment(avoid_m4$full_model) %>%
  dplyr::select(preference_binary, trial_order , bird_id, .fitted)
RIS<-augment(avoid_m4s$full_model) %>%
  dplyr::select(preference_binary, trial_order , bird_id, .fitted)
RI$.fittedRIS<-RIS$.fitted

df <- RI %>%
  dplyr::group_by(bird_id, trial_order) %>%
  dplyr::summarise(RI = mean(.fitted),
                   RIS = mean(.fittedRIS))%>%
  gather(type, Value, `RI`:`RIS`)
df$bird_id <- as.character(df$bird_id)

#Supplementary Figure
ggplot(df[df$type == "RI",], aes(x = trial_order, y = Value, group = bird_id)) +
  geom_line() +
  theme_classic() +
  labs(y="", x="")+
  ggtitle("Random Intercept") + guides(color = guide_legend(nrow = 2, byrow = TRUE))


ggplot(df[df$type == "RIS",], aes(x = trial_order, y = Value, group = bird_id)) +
  geom_line() +
  theme_classic()+
  labs(y="", x="")+
  ggtitle("Random Intercept\n and Slope")

