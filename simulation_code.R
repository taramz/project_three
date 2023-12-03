library(riskCommunicator)
library(tidyverse)
library(tableone)
library(comprehenr)

# Load in pre-processed data and framingham models
framingham_df_men <- readRDS('framingham_men.rds')
framingham_df_women <- readRDS('framingham_women.rds')

mod_men <- glm(CVD~log(HDLC)+log(TOTCHOL)+log(AGE)+log(SYSBP_UT+1)+
                 log(SYSBP_T+1)+CURSMOKE+DIABETES, 
               data= framingham_df_men, family= "binomial")


mod_women <- glm(CVD~log(HDLC)+log(TOTCHOL)+log(AGE)+log(SYSBP_UT+1)+
                   log(SYSBP_T+1)+CURSMOKE+DIABETES, 
                 data= framingham_df_women, family= "binomial")

nhanes <- readRDS('NHANES.rds')
nhanes_eligible <- nhanes %>%
  filter(AGE >=28 & AGE <= 62)

diabetes_mod <- glm(DIABETES ~ AGE + HDLC + BPMEDS, data = nhanes_eligible, family = 'binomial')

# Code for the calculation of Brier scores and AUC
calculate_weights <- function(data) {
  #' Function to calculate the weights
  #' @param data, the data used for building the model
  #'
  #' @return vector of weights
  
  model <- glm(S ~ ., data = data, family = "binomial")
  S1 <- predict(model, type = "response")
  
  return((1-S1)/S1)
  
}

brier_score <- function(data, sex) {
  #' Function to calculate the Brier Score by sex
  #' @param data, the data to analyze
  #' @param sex, the sex to calculate the Brier Score for
  #' 
  #' @return brier score
  
  wt <- data %>%
    filter(SEX == sex & S == 1) %>%
    pull(weights)
  
  nhanes <- data %>%
    filter(SEX == sex & S == 0) %>%
    nrow()
  
  framingham <- data %>%
    filter(S == 1 & SEX == sex)
  
  if (sex == 1) {
    pred = predict(mod_men, type = 'response')
    return(sum(wt * (framingham$CVD - pred)^2)/nhanes)
  } else {
    pred = predict(mod_women, type = 'response')
    return(sum(wt * (framingham$CVD - pred)^2)/nhanes)
  }
}

# Initialize empty vector for male and female Brier scores for each imputed dataset
scores_male <- c()
scores_female <- c()
auc_male <- c()
auc_female <- c()

auc <- function(data, sex) {
  #' Function to calculate the AUC by sex
  #' @param data, the data to analyze
  #' @param sex, the sex to calculate the AUC for
  #' 
  #' @return AUC value
  
  data <- data %>%
    filter(SEX == sex & S == 1)
  
  if (sex == 1) {
    data$pred <- predict(mod_men, type = 'response')
  } else {
    data$pred <- predict(mod_women, type = 'response')
  }
  
  num <- 0
  denom <- 0

  for (i in 1:nrow(data)) {
    for (j in 1:nrow(data)) {
      
      if (i != j) {
        wt <- ((1 - data$weights[i]) * (1 - data$weights[j]))/(data$weights[i] * data$weights[j])
        
        indicator_n <- ifelse((data$pred[i] > data$pred[j]) & 
          (data$CVD[i] == 1 & data$CVD[j] == 0) &
          (data$S[i] == 1 & data$S[j] == 1), 1, 0)
        
        indicator_d <- ifelse((data$CVD[i] == 1 & data$CVD[j] == 0) &
          (data$S[i] == 1 & data$S[j] == 1), 1, 0)

        
        num <- num + (wt * indicator_n)
        denom <- denom + (wt * indicator_d)
      }
      
    }
  }
  return(num/denom)
}


# Loop through each imputed dataset, calculate the weights, and then the Brier scores
for (i in 1:m) {
  
  # Take the ith imputed dataset
  # Calculate variables so that it can be added to to the framingham data
  data <- complete(mice_data_out, i)
  data <- data %>%
    filter(AGE >= 28 & AGE <= 62) %>%
    dplyr::select(c("SEX", "TOTCHOL", "AGE", "BMI", "SYSBP", "HDLC", 
             "CURSMOKE", "BPMEDS", "DIABETES", "DIABP"))
  data$CVD <- NA
  
  data$SYSBP_UT <- ifelse(data$BPMEDS == 0, data$SYSBP, 0)
  data$SYSBP_T <- ifelse(data$BPMEDS ==1, data$SYSBP, 0)
  
  combined_df <- rbind(framingham_df, data) %>%
    mutate(S = case_when(is.na(CVD) ~ 0,
                         TRUE ~ 1))
  
  # Calculate weights and add it as a column in the data
  weights <- calculate_weights(combined_df %>% dplyr::select(-CVD))
  combined_df <- cbind(combined_df, weights)
  
  # Calculate the brier and AUC scores by gender
  scores_male <- c(scores_male, brier_score(combined_df, 1))
  scores_female <- c(scores_female, brier_score(combined_df, 2))
  
  auc_male <- c(auc_male, auc(combined_df, 1))
  auc_female <- c(auc_female, auc(combined_df, 2))
  
}


# Simulate two datasets - male and female
# TODO: Change Mean, SD values for the genders
# Male
set.seed(1)
n <- nrow(nhanes %>% filter(AGE >= 28 & AGE <= 62 & SEX == 1))
n_sim <- 500
cov <- seq(-5,5,1) # Varied Covariances
thresholds <- seq(120, 140, 5) # SYSBP thresholds
full_grid <- expand.grid(cov = cov, thresh = thresholds, iter = seq(1,n_sim))
full_grid$brier <- NA
full_grid$auc <- NA
for (i in 1:nrow(full_grid)) {
  
  # Generate the cholesterol variables
  # Vary the covariance based on the values in the full_grid dataframe
  chol_var <- mvrnorm(n, 
                       mu = c(46, 186),
                       Sigma = matrix(c(15^2, full_grid$cov[i],
                                        full_grid$cov[i], 43^2),
                                      2, 2))
  
  HDLC <- chol_var[,1]
  TOTCHOL <- chol_var[,2]
  
  # Generate the other continuous variables
  SYSBP <- rnorm(n, 118, 20)
  AGE <- round(runif(n, 28, 62))
  
  # Use the threshold in the full_grid dataframe to determine the value of bp meds
  BP_MEDS <- as.numeric(SYSBP >= full_grid$thresh[i])
  
  # Generate smoking variable using the binomial distribution
  CURSMOKE <- rbinom(n, 1, 0.25)
  SEX <- rep(1, n)
  
  data <- data.frame(HDLC = HDLC, TOTCHOL = TOTCHOL, SYSBP = SYSBP, AGE = AGE, 
                     BPMEDS = BP_MEDS, CURSMOKE = CURSMOKE, SEX = SEX)
  
  # Calculate probability of diabetes using the diabetes model
  x_vars <- model.matrix(CURSMOKE ~ AGE + HDLC + BPMEDS, data = data)
  probs <- x_vars %*% coef(diabetes_mod)
  probs <- exp(probs)/(1 + exp(probs))
  
  # Use probability to model diabetes as a binomial variable
  data$DIABETES <- to_vec(for(i in 1:length(probs)) rbinom(1,1,probs[i]))
  
  # Create variables that are in framingham data & the model
  data$CVD <- NA
  data$SYSBP_UT <- ifelse(data$BPMEDS == 0, data$SYSBP, 0)
  data$SYSBP_T <- ifelse(data$BPMEDS ==1, data$SYSBP, 0)
  
  # Join with the framingham data and calculate the brier score
  framingham <- framingham_df %>%
    dplyr::select(HDLC, TOTCHOL, SYSBP, AGE, BPMEDS, CURSMOKE, SYSBP_UT, SYSBP_T, CVD, DIABETES, SEX)
  combined_df <- rbind(framingham, data) %>%
    mutate(S = case_when(is.na(CVD) ~ 0,
                         TRUE ~ 1))
  weights <- calculate_weights(combined_df %>% dplyr::select(-CVD))
  combined_df <- cbind(combined_df, weights)
  full_grid$brier[i] <- brier_score(combined_df, 1)
}

#store results as RDS file
saveRDS(full_grid, 'sim_results_men.rds')

# Female
# Same code as for men - only differences are the number of values & the mean/SD
n <- nrow(nhanes %>% filter(AGE >= 28 & AGE <= 62 & SEX == 1))
set.seed(1)
n_sim <- 500
cov <- seq(-5,5,1) # Varied Covariances
thresholds <- seq(120, 140, 5) # SYSBP thresholds
full_grid <- expand.grid(cov = cov, thresh = thresholds, iter = seq(1,n_sim))
full_grid$brier <- NA
full_grid$auc <- NA
for (i in 1:nrow(full_grid)) {
  chol_var <- mvrnorm(n, 
                      mu = c(56, 201),
                      Sigma = matrix(c(17^2, full_grid$cov[i],
                                       full_grid$cov[i], 41^2),
                                     2, 2))
  
  HDLC <- chol_var[,1]
  TOTCHOL <- chol_var[,2]
  SYSBP <- rnorm(n, 134, 20)
  AGE <- round(runif(n, 28, 62))
  BP_MEDS <- as.numeric(SYSBP >= full_grid$thresh[i])
  CURSMOKE <- rbinom(n, 1, 0.19)
  SEX <- rep(2, n)
  
  data <- data.frame(HDLC = HDLC, TOTCHOL = TOTCHOL, SYSBP = SYSBP, AGE = AGE, 
                     BPMEDS = BP_MEDS, CURSMOKE = CURSMOKE, SEX = SEX)
  
  x_vars <- model.matrix(CURSMOKE ~ AGE + HDLC + BPMEDS, data = data)
  probs <- x_vars %*% coef(diabetes_mod)
  probs <- exp(probs)/(1 + exp(probs))
  
  data$DIABETES <- to_vec(for(i in 1:length(probs)) rbinom(1,1,probs[i]))
  
  data$CVD <- NA
  
  data$SYSBP_UT <- ifelse(data$BPMEDS == 0, data$SYSBP, 0)
  data$SYSBP_T <- ifelse(data$BPMEDS ==1, data$SYSBP, 0)
  
  framingham <- framingham_df %>%
    dplyr::select(HDLC, TOTCHOL, SYSBP, AGE, BPMEDS, CURSMOKE, SYSBP_UT, SYSBP_T, CVD, DIABETES, SEX)
  combined_df <- rbind(framingham, data) %>%
    mutate(S = case_when(is.na(CVD) ~ 0,
                         TRUE ~ 1))
  weights <- calculate_weights(combined_df %>% dplyr::select(-CVD))
  combined_df <- cbind(combined_df, weights)
  full_grid$brier[i] <- brier_score(combined_df, 2)
}

saveRDS(full_grid, 'sim_results_women.rds')

# Store diabetes model
saveRDS(diabetes_mod, 'diabetes_mod.rds')
