##########################################################################################################
############################################### Initial An ###############################################
##########################################################################################################

#----------------------------------------- LIBRARIES

library(tidyverse)  # for data manipulation and visualization
library(GGally)
library(ggplot2)
library(gridExtra)

#----------------------------------------- GET DATA

# Read a txt file
data <- read.table(file = "Report1_Dataset.txt") 
colnames(data) = c("id", "length_of_stay", "age", "infection_risk", "routine_culturing_ratio", "routine_chest_xray_ratio", "number_of_beds", "medical_school_affiliation", "region", "average_daily_census", "number_of_nurses", "available_facilities_and_services")
data <- data[,2:12]

# Convert the "school" and "region" variables to factors:
data_subset <- data
data_subset$medical_school_affiliation <- as.numeric(data_subset$medical_school_affiliation)
data_subset$region <- as.numeric(data_subset$region)

#---------------------Correlation
cor(data)
ggpairs(data_subset)

#----------------------------------------- Analysis of Infection Risk

## Histogram and qqnorm To know the distribution of infection_risk
hist(data$infection_risk, xlab = "Infection risk", main = "")
# Verify if is normal
qqnorm(data$infection_risk, main="")
qqline(data$infection_risk)

# média, desvio padrão, mediana e min, max, boxplot
summary(data$infection_risk)
sd(data$infection_risk)
length(data$infection_risk)

# Boxplot - Distribuição de valores e possíveis outliers
boxplot(data$infection_risk)


## Analisys per Region 
#Box plot of infection_risk by region
ggplot(data, aes(x=region, y=infection_risk, fill=region)) +
  geom_boxplot() +
  ggtitle("") +
  facet_wrap(~region,ncol = 4)

# Histogram of risk of infection by region
ggplot(data, aes(x=infection_risk, fill=region)) + 
  geom_histogram(bins=30) +
  ggtitle("Histogram of risk of infection by region")


## Analisys per medical_school_affiliation
#Box plot of infection_risk by medical_school_affiliation
ggplot(data, aes(x=medical_school_affiliation, y=infection_risk, fill=medical_school_affiliation)) +
  geom_boxplot() +
  ggtitle("") +
  facet_wrap(~medical_school_affiliation,ncol = 2)

# Histogram of risk of infection by school affiliation
ggplot(data, aes(x=infection_risk, fill=medical_school_affiliation)) + 
  geom_histogram(bins=30) +
  ggtitle("Histogram of risk of infection by school affiliation")


##########################################################################################################
############################################### MCMC #####################################################
##########################################################################################################

# Convert the "school" and "region" variables to numeric:
data$medical_school_affiliation <- as.numeric(data$medical_school_affiliation)
data$region <- as.numeric(data$region)


#Normalize the data
data[, c("length_of_stay", "age", "infection_risk", "routine_culturing_ratio", "routine_chest_xray_ratio", "number_of_beds", "average_daily_census", "number_of_nurses", "available_facilities_and_services")] <- scale(data[, c("length_of_stay", "age", "infection_risk", "routine_culturing_ratio", "routine_chest_xray_ratio", "number_of_beds", "average_daily_census", "number_of_nurses", "available_facilities_and_services")], center = TRUE, scale = TRUE)


#Intercept
library(fastDummies)
data_1 = data
intercept_1 = rep(c(1), times = 113)
data_1$intercept <- intercept_1
data_1 <- data_1[, c(12, 1:11)]

#Change to dummies
data_1 <- dummy_cols(data_1, select_columns = c("medical_school_affiliation", "region"))
data_1 <- subset(data_1, select = -c(medical_school_affiliation, medical_school_affiliation_2, region, region_4))

#Create a matrix with 1's
#data_1 <- model.matrix(~ length_of_stay + age + infection_risk + routine_culturing_ratio + routine_chest_xray_ratio + number_of_beds + medical_school_affiliation +  region +  average_daily_census + number_of_nurses +  available_facilities_and_services, data = data )


#Library
library(rjags) 
library(coda)

#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- MODEL 1 -----------------------------------------------#
#----------------------------------------------------------------------------------------------------#

Y = data_1$infection_risk
X1 = subset(data_1, select = -c(infection_risk))

# Specify the model
file <- "regModel.txt"

cat("model{
  for (i in 1:n){
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- inprod(beta[],X[i,])
  }
  
  for(j in 1:13){
    beta[j] ~ dnorm(0,0.0001)
  }
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)

}", file = file)


# Compile the model
model_1 <- jags.model(file = file, data = list(Y=Y, n=nrow(X1), X=X1), n.chains = 2)
update(model_1,1000)

# Run the MCMC sampling (posterior)
samples_1 <- coda.samples(model_1, variable.names = c("beta", "tau", "sigma"), n.iter = 10000, thin = 1, n.chains = 2)

# Summarize the results
summary(samples_1)

# Check for convergence
library(mcmcplots)
denplot(samples_1)
gelman.plot(samples_1)
traceplot(samples_1)
gelman.diag(samples_1)
heidel.diag(samples_1)
raftery.diag(samples_1)

#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- MODEL 2 -----------------------------------------------#
#----------------------------------------------------------------------------------------------------#

#Remove average_daily_census
data_2 <- subset(data_1, select = -c(average_daily_census))

Y = data_2$infection_risk
X2 = subset(data_2, select = -c(infection_risk))

# Specify the model
file <- "regModel.txt"

cat("model{
  for (i in 1:n){
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- inprod(beta[],X[i,])
    }

  for(j in 1:12){
    beta[j] ~ dnorm(0,0.0001)
  }
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}", file = file)

# Compile the model
model_2 <- jags.model(file = file, data = list(Y=Y, n=nrow(X2), X=X2), n.chains = 2)
update(model_2,1000)

# Run the MCMC sampling (posterior)
samples_2 <- coda.samples(model_2, variable.names = c("beta", "tau", "sigma"), n.iter = 10000, thin = 1, n.chains = 2)

# Summarize the results
summary(samples_2)

# Check for convergence
denplot(samples_2)
gelman.plot(samples_2)
traceplot(samples_2)
gelman.diag(samples_2)
heidel.diag(samples_2)
raftery.diag(samples_2)

#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- MODEL 3 -----------------------------------------------#
#----------------------------------------------------------------------------------------------------#

#Remove number_of_beds 
data_3 <- subset(data_1, select = -c(number_of_beds))

Y = data_3$infection_risk
X3 = subset(data_3, select = -c(infection_risk))

# Specify the model
file <- "regModel.txt"

cat("model{
  for (i in 1:n){
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- inprod(beta[],X[i,])
    }

  for(j in 1:12){
    beta[j] ~ dnorm(0,0.0001)
  }
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}", file = file)


# Compile the model
model_3 <- jags.model(file = file, data = list(Y=Y, n=nrow(X3), X=X3), n.chains = 2)
update(model_3,1000)

# Run the MCMC sampling (posterior)
samples_3 <- coda.samples(model_3, variable.names = c("beta", "tau", "sigma"), n.iter = 10000, thin = 1, n.chains = 2)

# Summarize the results
summary(samples_3)

# Check for convergence
denplot(samples_3)
gelman.plot(samples_3)
traceplot(samples_3)
gelman.diag(samples_3)
heidel.diag(samples_3)
raftery.diag(samples_3)

#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- MODEL 4 -----------------------------------------------#
#----------------------------------------------------------------------------------------------------#

#Remove number_of_beds and average_daily_census
data_4 <- subset(data_2, select = -c(number_of_beds))

Y = data_4$infection_risk
X4 = subset(data_4, select = -c(infection_risk))

# Specify the model
file <- "regModel.txt"

cat("model{
  for (i in 1:n){
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- inprod(beta[],X[i,])
    }

  for(j in 1:11){
    beta[j] ~ dnorm(0,0.0001)
  }
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}", file = file)


# Compile the model
model_4 <- jags.model(file = file, data = list(Y=Y, n=nrow(X4), X=X4), n.chains = 2)
update(model_4,1000)

# Run the MCMC sampling (posterior)
samples_4 <- coda.samples(model_4, variable.names = c("beta", "tau", "sigma"), n.iter = 10000, thin = 1, n.chains = 2)

# Summarize the results
summary(samples_4)

# Check for convergence
denplot(samples_4)
gelman.plot(samples_4)

par(mfrow=c(3,4))
traceplot(samples_4)
gelman.diag(samples_4)
heidel.diag(samples_4)
raftery.diag(samples_4)

#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- MODEL 5 -----------------------------------------------#
#----------------------------------------------------------------------------------------------------#

#Remove number_of_beds and number_of_nurses
data_5 <- subset(data_3, select = -c(number_of_nurses))

Y = data_5$infection_risk
X5 = subset(data_5, select = -c(infection_risk))

# Specify the model
file <- "regModel.txt"

cat("model{
  for (i in 1:n){
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- inprod(beta[],X[i,])
    }

  for(j in 1:11){
    beta[j] ~ dnorm(0,0.0001)
  }
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}", file = file)


# Compile the model
model_5 <- jags.model(file = file, data = list(Y=Y, n=nrow(X5), X=X5), n.chains = 2)
update(model_5,1000)

# Run the MCMC sampling (posterior)
samples_5 <- coda.samples(model_5, variable.names = c("beta", "tau", "sigma"), n.iter = 10000, thin = 1, n.chains = 2)

# Summarize the results
summary(samples_5)

# Check for convergence
denplot(samples_5)
gelman.plot(samples_5)
traceplot(samples_5)
gelman.diag(samples_5)
heidel.diag(samples_5)
raftery.diag(samples_5)

#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- MODEL 6 -----------------------------------------------#
#----------------------------------------------------------------------------------------------------#

#Remove average_daily_census and number_of_nurses
data_6 <- subset(data_2, select = -c(number_of_nurses))

Y = data_6$infection_risk
X6 = subset(data_6, select = -c(infection_risk))

# Specify the model
file <- "regModel.txt"

cat("model{
  for (i in 1:n){
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- inprod(beta[],X[i,])
    }

  for(j in 1:11){
    beta[j] ~ dnorm(0,0.0001)
  }
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}", file = file)


# Compile the model
model_6 <- jags.model(file = file, data = list(Y=Y, n=nrow(X6), X=X6), n.chains = 2)
update(model_6,1000)

# Run the MCMC sampling (posterior)
samples_6 <- coda.samples(model_6, variable.names = c("beta", "tau", "sigma"), n.iter = 10000, thin = 1, n.chains = 2)

# Summarize the results
summary(samples_6)

# Check for convergence
denplot(samples_6)
gelman.plot(samples_6)
traceplot(samples_6)
gelman.diag(samples_6)
heidel.diag(samples_6)
raftery.diag(samples_6)

#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- COMPARE MODELS ----------------------------------------#
#----------------------------------------------------------------------------------------------------#

# Calculate the DIC for model_1
#a = as.mcmc(model_1)
dic1 <- dic.samples(model_1, n.iter = 10000)
dic1
235.2 + 14.32

# Calculate the DIC for model_2
dic2 <- dic.samples(model_2, n.iter = 10000)
dic2
235.4 + 13.2


# Calculate the DIC for model_3
dic3 <- dic.samples(model_3, n.iter = 10000)
dic3
236 + 13.35

# Calculate the DIC for model_4
dic4 <- dic.samples(model_4, n.iter = 10000)
dic4
234.8 + 12.28

# Calculate the DIC for model_5
dic5 <- dic.samples(model_5, n.iter = 10000)
dic5
235.6 + 12.27

# Calculate the DIC for model_5
dic6 <- dic.samples(model_6, n.iter = 10000)
dic6
236.2 + 12.18

# Compare the DIC values
#if (dic1 < dic2) {
#  print("Model 1 is preferred according to DIC")
#} else {
#  print("Model 2 is preferred according to DIC")
#}

#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- Evaluate Best Model -----------------------------------#
#----------------------------------------------------------------------------------------------------#

#########################
### Model Diagnostics ###
#########################

# Summarize the results
summary(samples_4)

# Densitylot the results
densplot(samples_4)

#Autocorrealtion
acfplot(samples_4)
acfplot(samples_1)
acfplot(samples_2)
acfplot(samples_3)
acfplot(samples_5)
acfplot(samples_6)
 
# Check normality of residuals
library(lattice)

chosen_model <- lm(infection_risk ~ length_of_stay + age + routine_culturing_ratio + routine_chest_xray_ratio + medical_school_affiliation +  region  + number_of_nurses +  available_facilities_and_services, data = data_subset )

residuals <- residuals(chosen_model)

qqmath(resid(chosen_model))
hist(residuals)


# Calculate the HPD intervals for the coefficients
coef1_hpd <- HPDinterval(samples_4, prob = 0.95)

# Plot the HPD intervals for the coefficients

library("bayesplot")
library("mcmcse")

mcmc_intervals(samples_4)
 
 
#######################################################################################################
#######################################################################################################
##########################################NON-NORMALIZED###############################################
#######################################################################################################
#######################################################################################################

# Read a txt file
newdata <- read.table(file = "Report1_Dataset.txt") 
colnames(newdata) = c("id", "length_of_stay", "age", "infection_risk", "routine_culturing_ratio", "routine_chest_xray_ratio", "number_of_beds", "medical_school_affiliation", "region", "average_daily_census", "number_of_nurses", "available_facilities_and_services")
newdata<-newdata[,2:12]

# Convert the "school" and "region" variables to numeric:
newdata$medical_school_affiliation <- as.numeric(newdata$medical_school_affiliation)
newdata$region <- as.numeric(newdata$region)

#summary of data
summary(newdata)

#Intercept
library(fastDummies)
newdata_1 = newdata
intercept_1 = rep(c(1), times = 113)
newdata_1$intercept <- intercept_1
newdata_1 <- newdata_1[, c(12, 1:11)]

#Change to dummies
newdata_1 <- dummy_cols(newdata_1, select_columns = c("medical_school_affiliation", "region"))
newdata_1 <- subset(newdata_1, select = -c(medical_school_affiliation, medical_school_affiliation_2, region, region_4))

#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- MODEL 1 -----------------------------------------------#
#----------------------------------------------------------------------------------------------------#

nY = newdata_1$infection_risk
nX1 = subset(newdata_1, select = -c(infection_risk))

# Specify the model
file <- "regModel.txt"

cat("model{
  for (i in 1:n){
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- inprod(beta[],X[i,])
  }
  
  for(j in 1:13){
    beta[j] ~ dnorm(0,0.0001)
  }
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)

}", file = file)


# Compile the model
newmodel_1 <- jags.model(file = file, data = list(Y=nY, n=nrow(nX1), X=nX1), n.chains = 2)
update(newmodel_1,1000)

# Run the MCMC sampling (posterior)
newsamples_1 <- coda.samples(newmodel_1, variable.names = c("beta", "tau", "sigma"), n.iter = 10000, thin = 1, n.chains = 2)

#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- MODEL 2 -----------------------------------------------#
#----------------------------------------------------------------------------------------------------#

#Remove average_daily_census
newdata_2 <- subset(newdata_1, select = -c(average_daily_census))

nY = newdata_2$infection_risk
nX2 = subset(newdata_2, select = -c(infection_risk))

# Specify the model
file <- "regModel.txt"

cat("model{
  for (i in 1:n){
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- inprod(beta[],X[i,])
    }

  for(j in 1:12){
    beta[j] ~ dnorm(0,0.0001)
  }
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}", file = file)

# Compile the model
newmodel_2 <- jags.model(file = file, data = list(Y=nY, n=nrow(nX2), X=nX2), n.chains = 2)
update(newmodel_2,1000)

# Run the MCMC sampling (posterior)
newsamples_2 <- coda.samples(newmodel_2, variable.names = c("beta", "tau", "sigma"), n.iter = 10000, thin = 1, n.chains = 2)

#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- MODEL 3 -----------------------------------------------#
#----------------------------------------------------------------------------------------------------#

#Remove number_of_beds 
newdata_3 <- subset(newdata_1, select = -c(number_of_beds))

nY = newdata_3$infection_risk
nX3 = subset(newdata_3, select = -c(infection_risk))

# Specify the model
file <- "regModel.txt"

cat("model{
  for (i in 1:n){
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- inprod(beta[],X[i,])
    }

  for(j in 1:12){
    beta[j] ~ dnorm(0,0.0001)
  }
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}", file = file)


# Compile the model
newmodel_3 <- jags.model(file = file, data = list(Y=nY, n=nrow(nX3), X=nX3), n.chains = 2)
update(newmodel_3,1000)

# Run the MCMC sampling (posterior)
newsamples_3 <- coda.samples(newmodel_3, variable.names = c("beta", "tau", "sigma"), n.iter = 10000, thin = 1, n.chains = 2)

#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- MODEL 4 -----------------------------------------------#
#----------------------------------------------------------------------------------------------------#

#Remove number_of_beds and average_daily_census
newdata_4 <- subset(newdata_2, select = -c(number_of_beds))

nY = newdata_4$infection_risk
nX4 = subset(newdata_4, select = -c(infection_risk))

# Specify the model
file <- "regModel.txt"

cat("model{
  for (i in 1:n){
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- inprod(beta[],X[i,])
    }

  for(j in 1:11){
    beta[j] ~ dnorm(0,0.0001)
  }
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}", file = file)


# Compile the model
newmodel_4 <- jags.model(file = file, data = list(Y=nY, n=nrow(nX4), X=nX4), n.chains = 2)
update(newmodel_4,1000)

# Run the MCMC sampling (posterior)
newsamples_4 <- coda.samples(newmodel_4, variable.names = c("beta", "tau", "sigma"), n.iter = 10000, thin = 1, n.chains = 2)

#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- MODEL 5 -----------------------------------------------#
#----------------------------------------------------------------------------------------------------#

#Remove number_of_beds and number_of_nurses
newdata_5 <- subset(newdata_3, select = -c(number_of_nurses))

nY = newdata_5$infection_risk
nX5 = subset(newdata_5, select = -c(infection_risk))

# Specify the model
file <- "regModel.txt"

cat("model{
  for (i in 1:n){
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- inprod(beta[],X[i,])
    }

  for(j in 1:11){
    beta[j] ~ dnorm(0,0.0001)
  }
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}", file = file)


# Compile the model
newmodel_5 <- jags.model(file = file, data = list(Y=nY, n=nrow(nX5), X=nX5), n.chains = 2)
update(model_5,1000)

# Run the MCMC sampling (posterior)
newsamples_5 <- coda.samples(newmodel_5, variable.names = c("beta", "tau", "sigma"), n.iter = 10000, thin = 1, n.chains = 2)

#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- MODEL 6 -----------------------------------------------#
#----------------------------------------------------------------------------------------------------#

#Remove average_daily_census and number_of_nurses
newdata_6 <- subset(newdata_2, select = -c(number_of_nurses))

nY = newdata_6$infection_risk
nX6 = subset(newdata_6, select = -c(infection_risk))

# Specify the model
file <- "regModel.txt"

cat("model{
  for (i in 1:n){
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- inprod(beta[],X[i,])
    }

  for(j in 1:11){
    beta[j] ~ dnorm(0,0.0001)
  }
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}", file = file)


# Compile the model
newmodel_6 <- jags.model(file = file, data = list(Y=nY, n=nrow(nX6), X=nX6), n.chains = 2)
update(model_6,1000)

# Run the MCMC sampling (posterior)
newsamples_6 <- coda.samples(newmodel_6, variable.names = c("beta", "tau", "sigma"), n.iter = 10000, thin = 1, n.chains = 2)


#----------------------------------------------------------------------------------------------------#
#-------------------------------------------- Evaluate Best Model -----------------------------------#
#----------------------------------------------------------------------------------------------------#

########################
### Model prediction ###
########################

#Prediction Model 4

covariates_model4 = nX4

nX4 = covariates_model4
lambda4 <- c()

summary_4 <- summary(newsamples_4)
mean_4 <- summary_4$statistics[1:11]
mean_4

betas_model4 <- c(-0.693607547,  0.258652022,  0.011637218,  0.052943469,  0.011618166,  0.001371684,  0.016968457, -0.599630752,-1.081340028, -0.731781495, -0.779306389)

for (i in 1:nrow(nX4)){
  lambda4[i] <-  betas_model4%*%as.numeric(covariates_model4[i,])
}

summary4 <- summary(lambda4) #For the mean
summary(newdata_1$infection_risk)
summary4

var(lambda4) #Variance
var(newdata_1$infection_risk)

#Prediction Model 5

covariates_model5 = nX5

nX5 = covariates_model5
lambda5 <- c()

summary_5 <- summary(newsamples_5)
mean_5 <- summary_5$statistics[1:11]
mean_5

betas_model5 <- c(-0.6929365451,  0.2430707714,  0.0118189989,  0.0544467904,  0.0120337105,  0.0009621916,  0.0200972269,-0.5918019200, -1.0748214459, -0.7389604371, -0.7821436689)

for (i in 1:nrow(nX5)){
  lambda5[i] <-  betas_model5%*%as.numeric(covariates_model5[i,])
}

summary5 <- summary(lambda5) #Mean
summary(newdata_1$infection_risk)
summary5

var(lambda5) #Variance
var(newdata_1$infection_risk)

#Prediction Model 6

covariates_model6 = nX6

nX6 = covariates_model6
lambda6 <- c()

summary_6 <- summary(newsamples_6)
mean_6 <- summary_6$statistics[1:11]
mean_6

betas_model6 <- c(-0.5025645994,  0.2589172172,  0.0068131002,  0.0534312636,  0.0114127375,  0.0004535454,  0.0213170843,-0.5344021502, -1.0732968685, -0.7491996798, -0.7734478015)

for (i in 1:nrow(nX6)){
  lambda6[i] <-  betas_model6%*%as.numeric(covariates_model6[i,])
}

summary6 <- summary(lambda6)

summary(newdata_1$infection_risk)
summary6

var(lambda6) #Variance
var(newdata_1$infection_risk)


#Prediction Model 1

covariates_model1 = nX1

nX1 = covariates_model1
lambda1 <- c()

summary_1 <- summary(newsamples_1)
mean_1 <- summary_1$statistics[1:13]
mean_1

betas_model1 <- c(-0.841825140,  0.239606828,  0.016257379,  0.054831762,  0.011544994, -0.003654753,  0.004039130,  0.001868768, 0.020426134, -0.664416272, -1.152338247, -0.717172854, -0.783174410)

for (i in 1:nrow(nX1)){
  lambda1[i] <-  betas_model1%*%as.numeric(covariates_model1[i,])
}

summary1 <- summary(lambda1)

summary(newdata_1$infection_risk)
summary1

var(lambda1) #Variance
var(newdata_1$infection_risk)

