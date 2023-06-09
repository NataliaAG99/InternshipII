#########################################
#Estimating Genetic Liabilities
#########################################
#### Preparing the environment ####
rm(list=ls(all=TRUE))
set.seed(270922)
#### Install and load packages ####
## First specify the packages of interest
packages = c("lavaan", "MASS", "psych", "dplyr", "polycor")
## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
# 
#### Load the data ####
# Import the data and look at the first six rows
dades1 <- read.csv(file = 'data_cohort1_students.csv', header = T, sep = ";")
dades1 <- subset(dades1, !is.na(dades1$mh_prob_fam_1)) # excluding the individual with NA for family history
dades2 <- read.csv(file = 'data_cohort2_students.csv', header = T, sep = ";")
dades2 <- subset(dades2, !is.na(dades2$mh_prob_fam_1)) # excluding the individual with NA for family history
dades2['X'] <- c(448:896)
dades <- rbind(dades1,dades2)
head(dades)
describe(dades)
summary(dades)

#Subset only the used variables
data_fh <- subset(dades, select=c(X, age, sex, gender,  mh_prob_fam_1, mh_prob_ever_1,
                                  mh_prob_fam_spec_1_1, mh_prob_fam_spec_1_2,
                                  mh_prob_fam_spec_1_3, 
                                  mh_prob_ever_2, mh_prob_fam_spec_2_1, mh_prob_fam_spec_2_2,
                                  mh_prob_fam_spec_2_3,  mh_prob_ever_3,
                                  mh_prob_fam_spec_3_1, mh_prob_fam_spec_3_2,
                                  mh_prob_fam_spec_3_3,  mh_prob_ever_4,
                                  mh_prob_fam_spec_4_1, mh_prob_fam_spec_4_2,
                                  mh_prob_fam_spec_4_3,  mh_prob_ever_5,
                                  mh_prob_fam_spec_5_1, mh_prob_fam_spec_5_2,
                                  mh_prob_fam_spec_5_3,  mh_prob_ever_6,
                                  mh_prob_fam_spec_6_1, mh_prob_fam_spec_6_2,
                                  mh_prob_fam_spec_6_3,  mh_prob_ever_7,
                                  mh_prob_fam_spec_7_1, mh_prob_fam_spec_7_2,
                                  mh_prob_fam_spec_7_3,  mh_prob_ever_8,
                                  mh_prob_fam_spec_8_1, mh_prob_fam_spec_8_2,
                                  mh_prob_fam_spec_8_3, lidas_composite.recipientfirstname)) 

colnames(data_fh) <- c("ID", "age", "sex", "gender", "probdepr_fam", "probdepr_ever", 
                       "probdepr_pare", "probdepr_sibl", "probdepr_gp",
                       "probmbp_ever", "probmbp_pare", 
                       "probmbp_sibl", "probmbp_gp",  
                       "probpanic_ever", "probpanic_pare", "probpanic_sibl", 
                       "probpanic_gp",  "probanx_ever",
                       "probanx_pare", "probanx_sibl", "probanx_gp",
                       "probsubs_ever", "probsubs_pare", 
                       "probsubs_sibl", "probsubs_gp",  
                       "probeat_ever", "probeat_pare", "probeat_sibl", 
                       "probeat_gp",  "probPTSD_ever",
                       "probPTSD_pare", "probPTSD_sibl", "probPTSD_gp",
                       "probADHD_ever", "probADHD_pare", 
                       "probADHD_sibl", "probADHD_gp", "first_name")

head(data_fh)
describe(data_fh)
summary(data_fh)

#Change NAs for 0
data_fh <- data_fh %>% 
  mutate(probdepr_pare = coalesce(probdepr_pare, 0),
         probdepr_sibl = coalesce(probdepr_sibl, 0),
         probdepr_gp = coalesce(probdepr_gp, 0),
         probmbp_pare = coalesce(probmbp_pare, 0),
         probmbp_sibl = coalesce(probmbp_sibl, 0),
         probmbp_gp = coalesce(probmbp_gp, 0),
         probpanic_pare = coalesce(probpanic_pare, 0),
         probpanic_sibl = coalesce(probpanic_sibl, 0),
         probpanic_gp = coalesce(probpanic_gp, 0),
         probanx_pare = coalesce(probanx_pare, 0),
         probanx_sibl = coalesce(probanx_sibl, 0),
         probanx_gp = coalesce(probanx_gp, 0),
         probsubs_pare = coalesce(probsubs_pare, 0),
         probsubs_sibl = coalesce(probsubs_sibl, 0),
         probsubs_gp = coalesce(probsubs_gp, 0),
         probeat_pare = coalesce(probeat_pare, 0),
         probeat_sibl = coalesce(probeat_sibl, 0),
         probeat_gp = coalesce(probeat_gp, 0),
         probPTSD_pare = coalesce(probPTSD_pare, 0),
         probPTSD_sibl = coalesce(probPTSD_sibl, 0),
         probPTSD_gp = coalesce(probPTSD_gp, 0),
         probADHD_pare = coalesce(probADHD_pare, 0),
         probADHD_sibl = coalesce(probADHD_sibl, 0),
         probADHD_gp = coalesce(probADHD_gp, 0))

data_fh2 <- subset(data_fh, select=-c(ID, age, sex, gender, probdepr_fam, first_name)) 

real_observed  <- data_fh2 == 1 # make the true/false real matrix

#Genetic and environmental correlations between relatives
# focal, parent, sib, gp rg within triat:
corA <- matrix(c(1,  .5,     .5,    .25,
                 .5,   1,     .5,   .25,
                 .5,  .5,     1,     .25,
                 .25, .25, .25,   1),4,4)

corE <- diag(4)

# genetic correlations between disorders:
rg <- matrix(c(1,  .3452,     .79,    .52, .4518, .1584, .39, .538,
               .3452,   1,     .1381,   .13, .1798, .1887, .08, .1205,
               .79,  .1381,     1,     .55, .17, -.15, .4619, .1381, 
               .52, .13, .55,   1, .45, .13, .29, .26,
               .4518, .1798, .17, .45, 1, -.0052, .21, .5405,
               .1584, .1887, -.15, .13, -.0052, 1, .02, -.2451,
               .39, .08, .4619, .29, .21, .02, 1, .48,
               .538, .1205, .1381, .26, .5405, -.2451, .48, 1),8,8)

# heritability
h2 <- c(.4, .8, .3, .34, .5, .5, .4, .75)

# genetic covariance:
cov_g <- diag(sqrt(h2)) %*% rg %*% diag(sqrt(h2))

# environmental covariance:
cov_e <- diag(1-diag(cov_g))

simulationsize<-5000000
# liability, note the %x% which multiplies the relation betwene relatives with the rg to make that 32*32 matrix!
liability <- mvrnorm(n=simulationsize,mu=rep(0,32),Sigma = (cov_g%x%corA) + (cov_e%x%corE))

# 8  trait specific thresholds:
thresh <- c(.1828564, 2.126254, 0.1290458, -0.2400805, 1.457098, 0.9038566, 1.260222, 1.123652)
# repeat 4 times each because 4 relatives
thresholds <- rep(thresh,each=4) 

# make empty simulated observed matrix:
mcobserved <- matrix(NA,simulationsize,32)

# make simulated observed
for(i in 1:32){
  mcobserved[,i]  <- liability[,i] > thresholds[i]
}

#Estimate the liability without the individual's report
# mak
out <- matrix(NA,896,8)
for(i in 1:896){
  # here is the trick to omit each 4th column: -c(4*0:9+1)
  ak <- (mcobserved[,-c(4*0:7+1)]-.5) %*% (real_observed[i,-c(4*0:7+1)]-.5)
  
  out[i,1] <- mean(liability[ak==max(ak),1])
  out[i,2] <- mean(liability[ak==max(ak),5])
  out[i,3] <- mean(liability[ak==max(ak),9])
  out[i,4] <- mean(liability[ak==max(ak),13])
  out[i,5] <- mean(liability[ak==max(ak),17])
  out[i,6] <- mean(liability[ak==max(ak),21])
  out[i,7] <- mean(liability[ak==max(ak),25])
  out[i,8] <- mean(liability[ak==max(ak),29])
  
  print(i)
}

colnames(out) <- c("Depr", "BP", "PD", "Anx", "Subs", "Eat", "PTSD", "ADHD")
out <- as.data.frame(out)

#Create the matrix with only the estimates, the real status and the name
data_fhStudents <- subset(data_fh2, select=c(probdepr_ever, probmbp_ever, probpanic_ever, probanx_ever, probsubs_ever, probeat_ever, probPTSD_ever, probADHD_ever)) 
IDs <- subset(data_fh, select = first_name)
correlations_ <- cbind(out, data_fhStudents)
correlations <- cbind(IDs, correlations_)

#Logistic regression

glm.Dep <- glm(probdepr_ever ~ Depr, data = correlations, family = binomial)
summary(glm.Dep)
polyserial(correlations$Depr, correlations$probdepr_ever, std.err=T)
glm.BP <- glm(probmbp_ever ~ BP, data = correlations, family = binomial)
summary(glm.BP)
polyserial(correlations$BP, correlations$probmbp_ever, std.err=T)
glm.PD <- glm(probpanic_ever ~ PD, data = correlations, family = binomial)
summary(glm.PD)
polyserial(correlations$PD, correlations$probpanic_ever, std.err=T)
glm.Anx <- glm(probanx_ever ~ Anx, data = correlations, family = binomial)
summary(glm.Anx)
polyserial(correlations$Anx, correlations$probanx_ever, std.err=T)
glm.Subs <- glm(probsubs_ever ~ Subs, data = correlations, family = binomial)
summary(glm.Subs)
polyserial(correlations$Subs, correlations$probsubs_ever, std.err=T)
glm.Eat <- glm(probeat_ever ~ Eat, data = correlations, family = binomial)
summary(glm.Eat)
polyserial(correlations$Eat, correlations$probeat_ever, std.err=T)
glm.PTSD <- glm(probPTSD_ever ~ PTSD, data = correlations, family = binomial)
summary(glm.PTSD)
polyserial(correlations$PTSD, correlations$probPTSD_ever, std.err=T)
glm.ADHD <- glm(probADHD_ever ~ ADHD, data = correlations, family = binomial)
summary(glm.ADHD)
polyserial(correlations$ADHD, correlations$probADHD_ever, std.err=T)

table(correlations$Depr, correlations$probdepr_ever)
table(correlations$BP, correlations$probmbp_ever)
table(correlations$PD, correlations$probpanic_ever)
table(correlations$Anx, correlations$probanx_ever)
table(correlations$Subs, correlations$probsubs_ever)
table(correlations$Eat, correlations$probeat_ever)
table(correlations$PTSD, correlations$probPTSD_ever)
table(correlations$ADHD, correlations$probADHD_ever)

#Odds ratio
(ORDepr <- exp(cbind(coef(glm.Dep), confint(glm.Dep))))
(ORBP <- exp(cbind(coef(glm.BP), confint(glm.BP))))
(ORPD <- exp(cbind(coef(glm.PD), confint(glm.PD))))
(ORAnx <- exp(cbind(coef(glm.Anx), confint(glm.Anx))))
(ORSubs <- exp(cbind(coef(glm.Subs), confint(glm.Subs))))
(OREat <- exp(cbind(coef(glm.Eat), confint(glm.Eat))))
(ORPTSD <- exp(cbind(coef(glm.PTSD), confint(glm.PTSD))))
(ORADHD <- exp(cbind(coef(glm.ADHD), confint(glm.ADHD))))

#########################################
#PHQ9
#########################################
#### Install and load packages ####
## First specify the packages of interest
packages = c("lavaan", "MASS", "psych", "dplyr", "readr", "ggplot2", "ggpubr", "tidyverse", "magrittr", "tidyr")
#
#
## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
# 
#### Load the data ####
df <- read_csv("Combined_weekend_surveysC1_C2.csv")
View(df)

length(unique(df$user_id)) # users, should be around 835
length(unique(df$day)) # consecutive timepoints, around 19 total if people did all
tokeep<- c("phq_anhedonia", "phq_depressed", "phq_sleeping_little",
           "phq_sleeping_much","phq_tired", "phq_appetite", "phq_overeating",
           "phq_failure", "phq_concentration", "phq_slowly", "phq_restless",
           "phq_suic_harm")
df2<- subset(df, item %in% tokeep)
df2 <- subset(df2, select=-c(...1))

#reshape it in wide format
reshaped_df <- df2 %>%
  group_by(user_id, day, first_name) %>%
  mutate(row = row_number()) %>%
  spread(item, answer_id) %>%
  fill(everything()) %>%
  group_by(user_id, day, first_name) %>%
  slice(which.min(rowSums(is.na(.))))

#subset the used variables
data_ <- subset(reshaped_df, select=c(user_id, day,first_name, phq_anhedonia, phq_depressed, phq_sleeping_little,
                                      phq_sleeping_much,phq_tired, phq_appetite, phq_overeating,
                                      phq_failure, phq_concentration, phq_slowly, phq_restless,
                                      phq_suic_harm))

#combine the PHQ9 variables
data_$phq45_inhyp_w <- pmax(data_$phq_sleeping_little, data_$phq_sleeping_much)
data_$phq78_apea_w <- pmax(data_$phq_appetite, data_$phq_overeating)
data_$phq1112_retag_w <- pmax(data_$phq_slowly, data_$phq_restless)

data <- subset(data_, select=-c(phq_sleeping_little, phq_sleeping_much, phq_appetite, phq_overeating,
                                phq_slowly, phq_restless))

#Delete >70 missingness
users <- as.vector(unique(data$user_id))

dfmiss <- data.frame()
for (i in 1:832) {
  userid <- users[i]
  df1 <- filter(data, user_id == userid) # select user
  timescounter<- length(unique(data$day))
  temp <- matrix(NA, 1, 2)
  temp [,1] <- userid
  temp [,2] <- (sum(is.na(df1)))/(9*timescounter)
  dfmiss <- rbind(dfmiss, temp)
  
}
colnames(dfmiss)<- c("user_id", "Missingness")
miss70 <- subset(dfmiss, Missingness > .70) #subset with individuals with more than 70% of missing values

values_to_delete70 <- miss70$user_id
data70 <- data[!data$user_id %in% values_to_delete70, ] #subset with individuals without more than 70% of missing values

data70$sum <- rowSums(data70[ , c(4:12)], na.rm=TRUE) #Now the ones with NA, are changed to 0
data70["sum"][data70["sum"] == 0] <- NA

#scale the variable to make it comparable
sum_scaled <- scale(data70$sum)
data70 <- cbind(data70, sum_scaled)
colnames(data70) <- c("user_id", "day", "first_name", "phq_anhedonia", "phq_depressed", 
                      "phq_tired", "phq_failure", "phq_concentration", "phq_suic_harm", 
                      "phq45_inhyp_w", "phq78_apea_w", "phq1112_retag_w", "sum", "sum_scaled")
# Calculate mean of sum for each ID
means <- aggregate(sum_scaled ~ first_name, data = data70, FUN = mean)

# Define quantile function
my_quantile <- function(x) {
  quantile(x, probs = 0.90)  # Change the probability (probs) as needed (0.9 or 0.1)
}

# Calculate quantil of sum for each ID
quant <- aggregate(sum_scaled ~ first_name, data = data70, FUN = my_quantile)

# Define variance function
my_variance <- function(x) {
  var(x, na.rm = T)
}

# Calculate var of sum for each ID
varis <- aggregate(sum_scaled ~ first_name, data = data70, FUN = my_variance)

#Merge with the genetic liabilities
new_df1 <- means %>% inner_join(correlations, by = "first_name", multiple = "all")
new_df <- quant %>% inner_join(correlations, by = "first_name", multiple = "all")
new_df2 <- varis %>% inner_join(correlations, by = "first_name", multiple = "all")
colnames(new_df1) <- c("first_name", "mean", "Depr", "BP", "PD", "Anx", "Subs", 
                       "Eat", "PTSD", "ADHD", "probdepr_ever", "probmbp_ever", 
                       "probpanic_ever", "probanx_ever", "probsubs_ever", 
                       "probeat_ever", "probPTSD_ever", "probADHD_ever", "rank")

new_df3 <- varis %>% inner_join(new_df1, by = "first_name", multiple = "all")
colnames(new_df3) <- c("first_name", "variance", "mean", "Depr", "BP", "PD", "Anx", "Subs", 
                       "Eat", "PTSD", "ADHD", "probdepr_ever", "probmbp_ever", 
                       "probpanic_ever", "probanx_ever", "probsubs_ever", 
                       "probeat_ever", "probPTSD_ever", "probADHD_ever", "rank")

####################
#Plots
# Plot quantiles using ggplot2
ggplot(new_df, aes(x = sum_scaled, y = Depr)) +
  geom_point(color = "#5ab4ac", size = 1) +
  labs(title = "90% quantile of sum of PHQ9 by GL", x = "Quantile of PHQ9", y = "GL")

# Plot means using ggplot2
ggplot(new_df1, aes(x = sum_scaled, y = Depr)) +
  geom_point(color = "#5ab4ac", size = 1) +
  labs(title = "Mean of sum of PHQ9 by GL", x = "Mean of PHQ9", y = "GL")

# Plot variance using ggplot2
ggplot(new_df2, aes(x = sum_scaled, y = Depr)) +
  geom_point(color = "#5ab4ac", size = 1) +
  labs(title = "Variance of sum of PHQ9 by GL", x = "Variance of PHQ9", y = "GL")

#######################
#Correlations
# Correlation quantiles
cor.test(new_df$sum_scaled, new_df$Depr, method = "pearson")

# Correlation means
cor.test(new_df1$sum_scaled, new_df1$Depr, method = "pearson")

# Correlation variance
cor.test(new_df2$sum_scaled, new_df2$Depr, method = "pearson")

#######################
#Regressions
# Anova quantiles
ano.quaNA90 <- lm(sum_scaled ~ Depr, data = new_df)
summary(ano.quaNA90)

ano.quaNA10 <- lm(sum_scaled ~ Depr, data = new_df)
summary(ano.quaNA10)

# Anova means
ano.meaNA <- lm(sum_scaled ~ Depr, data = new_df1)
summary(ano.meaNA)

# Anova variance
ano.varNA <- lm(sum_scaled ~ Depr, data = new_df2)
summary(ano.varNA)

#Other regressions
ano.varNABP <- lm(sum_scaled ~ BP, data = new_df2)
summary(ano.varNABP)

ano.varNAAnx <- lm(sum_scaled ~ Anx, data = new_df2)
summary(ano.varNAAnx)

ano.varNAADHD <- lm(sum_scaled ~ ADHD, data = new_df2)
summary(ano.varNAADHD)

#Correlation between the variance and the mean
cor.test(new_df2$sum_scaled, new_df1$sum_scaled, method = "pearson")

ano.varmeanNA <- lm(variance ~ Depr + mean, data = new_df3)
summary(ano.varmeanNA)
ano.varmeanNABP <- lm(variance ~ BP + mean, data = new_df3)
summary(ano.varmeanNABP)
ano.varmeanNAAnx <- lm(variance ~ Anx + mean, data = new_df3)
summary(ano.varmeanNAAnx)
ano.varmeanNAADHD <- lm(variance ~ ADHD + mean, data = new_df3)
summary(ano.varmeanNAADHD)

#########################################
#NA
#########################################
#### Install and load packages ####
## First specify the packages of interest
packages = c("lavaan", "MASS", "psych", "dplyr", "readr", "ggplot2", "ggpubr", "tidyverse", "magrittr")
#
## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
# 
#### Load the data ####

df <- read_csv("C1_C2_daily_wide.csv")
View(df)

length(unique(df$user_id)) # users, should be around 855
length(unique(df$counter)) # consecutive timepoints, around 352 total if people did all

#Create the subset with the used variables
data <- subset(df, select=c(user_id, counter, sad, stressed, overwhelmed, nervous,annoyed, neg_thoughts, tired, first_name))

#Remove the >70 missingness
users <- as.vector(unique(data$user_id))

dfmiss <- data.frame()
for (i in 1:855) {
  userid <- users[i]
  df1 <- filter(data, user_id == userid) # select user
  timescounter<- length(unique(data$counter))
  temp <- matrix(NA, 1, 2)
  temp [,1] <- userid
  temp [,2] <- (sum(is.na(df1)))/(7*timescounter)
  dfmiss <- rbind(dfmiss, temp)
  
}
colnames(dfmiss)<- c("user_id", "Missingness")
miss70 <- subset(dfmiss, Missingness > .70) #subset with individuals with more than 70% of missing values

values_to_delete70 <- miss70$user_id
data70 <- data[!data$user_id %in% values_to_delete70, ] #subset with individuals without more than 70% of missing values

data70$sum <- rowSums(data70[ , c(3:9)], na.rm=TRUE) #Now the ones with NA, are changed to 0
data70["sum"][data70["sum"] == 0] <- NA

#Scale the variable to make it comparable
sum_scaled <- scale(data70$sum)
data70 <- cbind(data70, sum_scaled)

# Calculate mean of sum for each ID
means <- aggregate(sum_scaled ~ first_name, data = data70, FUN = mean)

# Define quantile function
my_quantile <- function(x) {
  quantile(x, probs = 0.90)  # Change the probability (probs) as needed (0.9 or 0.1)
}

# Calculate quantil of sum for each ID
quant <- aggregate(sum_scaled ~ first_name, data = data70, FUN = my_quantile)

# Define variance function
my_variance <- function(x) {
  var(x, na.rm = T)
}

# Calculate var of sum for each ID
varis <- aggregate(sum_scaled ~ first_name, data = data70, FUN = my_variance)

#Merge with genetic liabilities
new_df1 <- means %>% inner_join(correlations, by = "first_name", multiple = "all")
new_df <- quant %>% inner_join(correlations, by = "first_name", multiple = "all")
new_df2 <- varis %>% inner_join(correlations, by = "first_name", multiple = "all")
colnames(new_df1) <- c("first_name", "mean", "Depr", "BP", "PD", "Anx", "Subs", 
                       "Eat", "PTSD", "ADHD", "probdepr_ever", "probmbp_ever", 
                       "probpanic_ever", "probanx_ever", "probsubs_ever", 
                       "probeat_ever", "probPTSD_ever", "probADHD_ever", "rank")
new_df3 <- varis %>% inner_join(new_df1, by = "first_name", multiple = "all")
colnames(new_df3) <- c("first_name", "variance", "mean", "Depr", "BP", "PD", "Anx", "Subs", 
                       "Eat", "PTSD", "ADHD", "probdepr_ever", "probmbp_ever", 
                       "probpanic_ever", "probanx_ever", "probsubs_ever", 
                       "probeat_ever", "probPTSD_ever", "probADHD_ever", "rank")

####################
#Plots
# Plot quantiles using ggplot2
ggplot(new_df, aes(x = sum_scaled, y = Depr)) +
  geom_point(color = "#5ab4ac", size = 1) +
  labs(title = "90% quantile of sum of negative symptoms by GL", x = "Quantile of NS", y = "GL")

# Plot means using ggplot2
ggplot(new_df1, aes(x = sum_scaled, y = Depr)) +
  geom_point(color = "#5ab4ac", size = 1) +
  labs(title = "Mean of sum of negative symptoms by GL", x = "Mean of NS", y = "GL")

# Plot variance using ggplot2
ggplot(new_df2, aes(x = sum_scaled, y = Depr)) +
  geom_point(color = "#5ab4ac", size = 1) +
  labs(title = "Variance of sum of negative symptoms by GL", x = "Variance of NS", y = "GL")

#######################
#Correlations
# Correlation quantiles
cor.test(new_df$sum_scaled, new_df$Depr, method = "pearson")

# Correlation means
cor.test(new_df1$sum_scaled, new_df1$Depr, method = "pearson")

# Correlation variance
cor.test(new_df2$sum_scaled, new_df2$Depr, method = "pearson")

#######################
#Regressions
# Anova quantiles
ano.quaNA90 <- lm(sum_scaled ~ Depr, data = new_df)
summary(ano.quaNA90)

ano.quaNA10 <- lm(sum_scaled ~ Depr, data = new_df)
summary(ano.quaNA10)

# Anova means
ano.meaNA <- lm(sum_scaled ~ Depr, data = new_df1)
summary(ano.meaNA)

# Anova variance
ano.varNA <- lm(sum_scaled ~ Depr, data = new_df2)
summary(ano.varNA)

#Other regressions
ano.varNABP <- lm(sum_scaled ~ BP, data = new_df2)
summary(ano.varNABP)

ano.varNAAnx <- lm(sum_scaled ~ Anx, data = new_df2)
summary(ano.varNAAnx)

ano.varNAADHD <- lm(sum_scaled ~ ADHD, data = new_df2)
summary(ano.varNAADHD)

#Correlation between the variance and the mean
cor.test(new_df2$sum_scaled, new_df1$sum_scaled, method = "pearson")

ano.varmeanNA <- lm(variance ~ Depr + mean, data = new_df3)
summary(ano.varmeanNA)
ano.varmeanNABP <- lm(variance ~ BP + mean, data = new_df3)
summary(ano.varmeanNABP)
ano.varmeanNAAnx <- lm(variance ~ Anx + mean, data = new_df3)
summary(ano.varmeanNAAnx)
ano.varmeanNAADHD <- lm(variance ~ ADHD + mean, data = new_df3)
summary(ano.varmeanNAADHD)

#########################################
#PA
#########################################
#### Install and load packages ####
## First specify the packages of interest
packages = c("lavaan", "MASS", "psych", "dplyr", "readr", "ggplot2", "ggpubr", "tidyverse", "magrittr")
#
## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
# 
#### Load the data ####

df <- read_csv("C1_C2_daily_wide.csv")
View(df)

length(unique(df$user_id)) # users, should be around 855
length(unique(df$counter)) # consecutive timepoints, around 352 total if people did all

#Subset with the used variables
data <- subset(df, select=c(user_id, counter, relaxed, happy, motivated, first_name))

#Remove the >70% missingness
users <- as.vector(unique(data$user_id))

dfmiss <- data.frame()
for (i in 1:855) {
  userid <- users[i]
  df1 <- filter(data, user_id == userid) # select user
  timescounter<- length(unique(data$counter))
  temp <- matrix(NA, 1, 2)
  temp [,1] <- userid
  temp [,2] <- (sum(is.na(df1)))/(3*timescounter)
  dfmiss <- rbind(dfmiss, temp)
  
}
colnames(dfmiss)<- c("user_id", "Missingness")
miss70 <- subset(dfmiss, Missingness > .70) #subset with individuals with more than 70% of missing values

values_to_delete70 <- miss70$user_id
data70 <- data[!data$user_id %in% values_to_delete70, ] #subset with individuals without more than 70% of missing values

data70$sum <- rowSums(data70[ , c(3:5)], na.rm=TRUE) #Now the ones with NA, are changed to 0
data70["sum"][data70["sum"] == 0] <- NA

#Scale the variable to make it comparable
sum_scaled <- scale(data70$sum)
data70 <- cbind(data70, sum_scaled)

# Calculate mean of sum for each ID
means <- aggregate(sum_scaled ~ first_name, data = data70, FUN = mean)

# Define quantile function
my_quantile <- function(x) {
  quantile(x, probs = 0.90)  # Change the probability (probs) as needed (0.9 or 0.1)
}

# Calculate quantil of sum for each ID
quant <- aggregate(sum_scaled ~ first_name, data = data70, FUN = my_quantile)

# Define variance function
my_variance <- function(x) {
  var(x, na.rm = T)
}

# Calculate var of sum for each ID
varis <- aggregate(sum_scaled ~ first_name, data = data70, FUN = my_variance)

#Merge with genetic liabilities
new_df1 <- means %>% inner_join(correlations, by = "first_name", multiple = "all")
new_df <- quant %>% inner_join(correlations, by = "first_name", multiple = "all")
new_df2 <- varis %>% inner_join(correlations, by = "first_name", multiple = "all")
colnames(new_df1) <- c("first_name", "mean", "Depr", "BP", "PD", "Anx", "Subs", 
                       "Eat", "PTSD", "ADHD", "probdepr_ever", "probmbp_ever", 
                       "probpanic_ever", "probanx_ever", "probsubs_ever", 
                       "probeat_ever", "probPTSD_ever", "probADHD_ever", "rank")
new_df3 <- varis %>% inner_join(new_df1, by = "first_name", multiple = "all")
colnames(new_df3) <- c("first_name", "variance", "mean", "Depr", "BP", "PD", "Anx", "Subs", 
                       "Eat", "PTSD", "ADHD", "probdepr_ever", "probmbp_ever", 
                       "probpanic_ever", "probanx_ever", "probsubs_ever", 
                       "probeat_ever", "probPTSD_ever", "probADHD_ever", "rank")


####################
#Plots
# Plot quantiles using ggplot2
ggplot(new_df, aes(x = sum_scaled, y = Depr)) +
  geom_point(color = "#5ab4ac", size = 1) +
  labs(title = "90% quantile of sum of positive symptoms by GL", x = "Quantile of NS", y = "GL")

# Plot means using ggplot2

ggplot(new_df1, aes(x = sum_scaled, y = Depr)) +
  geom_point(color = "#5ab4ac", size = 1) +
  labs(title = "Mean of sum of positive symptoms by GL", x = "Mean of NS", y = "GL")

# Plot variance using ggplot2

ggplot(new_df2, aes(x = sum_scaled, y = Depr)) +
  geom_point(color = "#5ab4ac", size = 1) +
  labs(title = "Variance of sum of positive symptoms by GL", x = "Variance of PA", y = "GL")

#######################
#Correlations
# Correlation quantiles
cor.test(new_df$sum_scaled, new_df$Depr, method = "pearson")

# Correlation means
cor.test(new_df1$sum_scaled, new_df1$Depr, method = "pearson")

# Correlation variance
cor.test(new_df2$sum_scaled, new_df2$Depr, method = "pearson")

#######################
#Regressions
# Anova quantiles
ano.quaPA90 <- lm(sum_scaled ~ Depr, data = new_df)
summary(ano.quaPA90)

ano.quaPA10 <- lm(sum_scaled ~ Depr, data = new_df)
summary(ano.quaPA10)

# Anova means
ano.meaPA <- lm(sum_scaled ~ Depr, data = new_df1)
summary(ano.meaPA)

# Anova variance
ano.varPA <- lm(sum_scaled ~ Depr, data = new_df2)
summary(ano.varPA)

#Other regressions
ano.varPABP <- lm(sum_scaled ~ BP, data = new_df2)
summary(ano.varPABP)

ano.varPAAnx <- lm(sum_scaled ~ Anx, data = new_df2)
summary(ano.varPAAnx)

ano.varPAADHD <- lm(sum_scaled ~ ADHD, data = new_df2)
summary(ano.varPAADHD)

#Correlation between the variance and the mean
cor.test(new_df2$sum_scaled, new_df1$sum_scaled, method = "pearson")

ano.varmeanPA <- lm(variance ~ Depr + mean, data = new_df3)
summary(ano.varmeanPA)
ano.varmeanPABP <- lm(variance ~ BP + mean, data = new_df3)
summary(ano.varmeanPABP)
ano.varmeanPAAnx <- lm(variance ~ Anx + mean, data = new_df3)
summary(ano.varmeanPAAnx)
ano.varmeanPAADHD <- lm(variance ~ ADHD + mean, data = new_df3)
summary(ano.varmeanPAADHD)

#########################################
#Network model
#########################################
#### Install and load packages ####
## First specify the packages of interest
packages = c("graphicalVAR", "readr", "dplyr", "qgraph", "imputeTS", "ggplot2")
#
#
## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
# 
#### Load the data ####

df <- read_csv("C1_C2_daily_wide.csv")
View(df)

length(unique(df$user_id)) # users, should be around 855
length(unique(df$counter)) # consecutive timepoints, around 352 total if people did all

vars <- c("sad", "stressed", "overwhelmed", "nervous","annoyed", "neg_thoughts", "tired") # select some items with variance

#### Missingness ####
data <- subset(df, select=c(user_id, counter, sad, stressed, overwhelmed, nervous,annoyed, neg_thoughts, tired, first_name))

users <- as.vector(unique(data$user_id))

dfmiss <- data.frame()
for (i in 1:855) {
  userid <- users[i]
  df1 <- filter(data, user_id == userid) # select user
  timescounter<- length(unique(data$counter))
  temp <- matrix(NA, 1, 2)
  temp [,1] <- userid
  temp [,2] <- (sum(is.na(df1)))/(7*timescounter)
  dfmiss <- rbind(dfmiss, temp)
  
}
colnames(dfmiss)<- c("user_id", "Missingness")
miss70 <- subset(dfmiss, Missingness > .70) #subset with individuals with more than 70% of missing values

values_to_delete70 <- miss70$user_id
data70 <- data[!data$user_id %in% values_to_delete70, ] #subset with individuals without more than 70% of missing values

#### Idiographic ####
# For one individual
df1 <- filter(data, user_id == 44199) # select user
summary(df1[,c(3:9)]) # check items to select for preliminary analyses

vars <- c("sad", "stressed", "overwhelmed", "nervous","annoyed", "neg_thoughts", "tired") # select some items with variance

nw44199 <- graphicalVAR(df1, 
                        beepvar = "counter",
                        gamma = 0,
                        deleteMissings = FALSE,
                        vars = vars)

plot(nw44199, theme="colorblind")

# Loop for all of them
graphic_list <- list()
unique_ids <- as.character(unique(data70$first_name))

for (i in unique_ids[1:624]) {
  tryCatch({
    print(i)
    df1 <- filter(data70, first_name == i) # select user
    nw <- graphicalVAR(df1, 
                       beepvar = "counter",
                       gamma = 0,
                       deleteMissings = FALSE,
                       vars = vars)
    graphic_list[[i]] <- nw 
  }, error = function(e) {
    print(paste("Error occurred for user", i))
    print(e)
  })
}

graphic_list$`44155` #to visualize one specific user

#### Extract the connectivity parameter ####
# Create an empty list to store the results
Connectivity <- list()

# Loop through each element in graphic_list
for (first_name in names(graphic_list)) {
  # Extract the PCC element for the current user
  PCCMatrix <- graphic_list[[first_name]]$PCC
  # Calculate the sum of the upper triangle of the PCC matrix
  sum_upper <- sum(PCCMatrix[upper.tri(PCCMatrix)])
  # Add the user ID and sum_upper to the results list
  Connectivity[[first_name]] <- list(first_name = first_name, sum_upper = sum_upper)
}

# Convert the results list to a data frame
Connectivity_df <- do.call(rbind, lapply(Connectivity, data.frame))

means <- aggregate(sum_upper ~ first_name, data = Connectivity_df, FUN = mean)

#Add the GL
new_df <- means %>% inner_join(correlations, by = "first_name", multiple = "all")

####################
#Plots
# Plot using ggplot2
ggplot(new_df, aes(x = sum_upper, y = Depr)) +
  geom_point(color = "#5ab4ac", size = 1) +
  labs(title = "Network connectivity by GL", x = "Connectivity", y = "GL")

#######################
#Correlations
# Correlation 
cor.test(new_df$sum_upper, new_df$Depr, method = "pearson")

#######################
#Regressions
# Anova 
ano.test <- lm(sum_upper ~ Depr, data = new_df)
summary(ano.test)
