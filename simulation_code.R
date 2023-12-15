# Load in packages
library(riskCommunicator)
library(tidyverse)
library(tableone)
library(kableExtra)
library(comprehenr)

data("framingham")

# The Framingham data has been used to create models for cardiovascular risk.
# The variable selection and model below are designed to mimic the models used
# in the paper General Cardiovascular Risk Profile for Use in Primary Care 
# This paper is available (cvd_risk_profile.pdf) on Canvas.

framingham_df <- framingham %>% dplyr::select(c(CVD, TIMECVD, SEX, TOTCHOL, AGE,
                                         SYSBP, DIABP, CURSMOKE, DIABETES, BPMEDS,
                                         HDLC, BMI))
framingham_df <- na.omit(framingham_df)

#CreateTableOne(data=framingham_df, strata = c("SEX"))

# Get blood pressure based on whether or not on BPMEDS
framingham_df$SYSBP_UT <- ifelse(framingham_df$BPMEDS == 0, 
                                 framingham_df$SYSBP, 0)
framingham_df$SYSBP_T <- ifelse(framingham_df$BPMEDS == 1, 
                                framingham_df$SYSBP, 0)

# Looking at risk within 15 years - remove censored data
#dim(framingham_df)
framingham_df <- framingham_df %>%
  filter(!(CVD == 0 & TIMECVD <= 365*15)) %>%
  dplyr::select(-c(TIMECVD))
#dim(framingham_df)

# Filter to each sex
framingham_df_men <- framingham_df %>% filter(SEX == 1)
framingham_df_women <- framingham_df %>% filter(SEX == 2)

# Fit models with log transforms for all continuous variables
mod_men <- glm(CVD~log(HDLC)+log(TOTCHOL)+log(AGE)+log(SYSBP_UT+1)+
                 log(SYSBP_T+1)+CURSMOKE+DIABETES, 
               data= framingham_df_men, family= "binomial")


mod_women <- glm(CVD~log(HDLC)+log(TOTCHOL)+log(AGE)+log(SYSBP_UT+1)+
                   log(SYSBP_T+1)+CURSMOKE+DIABETES, 
                 data= framingham_df_women, family= "binomial")




library(nhanesA)

# blood pressure, demographic, bmi, smoking, and hypertension info
bpx_2017 <- nhanes("BPX_J") %>% 
  dplyr::select(SEQN, BPXSY1, BPXDI1) %>% 
  rename(SYSBP = BPXSY1, DIABP = BPXDI1)
demo_2017 <- nhanes("DEMO_J") %>% 
  dplyr::select(SEQN, RIAGENDR, RIDAGEYR) %>% 
  rename(SEX = RIAGENDR, AGE = RIDAGEYR)
bmx_2017 <- nhanes("BMX_J") %>% 
  dplyr::select(SEQN, BMXBMI) %>% 
  rename(BMI = BMXBMI)
smq_2017 <- nhanes("SMQ_J") %>%
  mutate(CURSMOKE = case_when(SMQ040 %in% c(1,2) ~ 1,
                              SMQ040 == 3 ~ 0, 
                              SMQ020 == 2 ~ 0)) %>%
  dplyr::select(SEQN, CURSMOKE)
bpq_2017 <- nhanes("BPQ_J") %>% 
  mutate(BPMEDS = case_when(
    BPQ020 == 2 ~ 0,
    BPQ040A == 2 ~ 0,
    BPQ050A == 1 ~ 1,
    TRUE ~ NA )) %>%
  dplyr::select(SEQN, BPMEDS) 
tchol_2017 <- nhanes("TCHOL_J") %>% 
  dplyr::select(SEQN, LBXTC) %>% 
  rename(TOTCHOL = LBXTC)
hdl_2017 <- nhanes("HDL_J") %>% 
  dplyr::select(SEQN, LBDHDD) %>% 
  rename(HDLC = LBDHDD)
diq_2017 <- nhanes("DIQ_J") %>% 
  mutate(DIABETES = case_when(DIQ010 == 1 ~ 1, 
                              DIQ010 %in% c(2,3) ~ 0, 
                              TRUE ~ NA)) %>%
  dplyr::select(SEQN, DIABETES) 

# Join data from different tables
df_2017 <- bpx_2017 %>%
  full_join(demo_2017, by = "SEQN") %>%
  full_join(bmx_2017, by = "SEQN") %>%
  full_join(hdl_2017, by = "SEQN") %>%
  full_join(smq_2017, by = "SEQN") %>%
  full_join(bpq_2017, by = "SEQN") %>%
  full_join(tchol_2017, by = "SEQN") %>%
  full_join(diq_2017, by = "SEQN")

complete_cases <- df_2017[complete.cases(df_2017),]
nhanes_eligible <- complete_cases %>%
  filter(AGE >= 30) %>%
  dplyr::select(c("SEX", "TOTCHOL", "AGE", "BMI", "SYSBP", "HDLC", 
                  "CURSMOKE", "BPMEDS", "DIABETES", "DIABP"))

diabetes_mod <- glm(DIABETES ~ AGE + HDLC + BPMEDS, data = framingham_df, family = 'binomial')

# Code for the calculation of Brier scores and AUC
calculate_weights <- function(data) {
  #' Function to calculate the weights
  #' @param data, the data used for building the model
  #'
  #' @return vector of weights
  
  model <- glm(S ~ ., data = data, family = "binomial") #same formula as mod_men, do a sanity check on simulated datasets
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
    filter(AGE >= 30) %>%
    dplyr::select(c("SEX", "TOTCHOL", "AGE", "BMI", "SYSBP", "HDLC", 
             "CURSMOKE", "BPMEDS", "DIABETES", "DIABP"))
  data$CVD <- NA
  
  data$SYSBP_UT <- ifelse(data$BPMEDS == 0, data$SYSBP, 0)
  data$SYSBP_T <- ifelse(data$BPMEDS ==1, data$SYSBP, 0)
  
  combined_df <- rbind(framingham_df, data) %>%
    mutate(S = case_when(is.na(CVD) ~ 0,
                         TRUE ~ 1))
  
  # Calculate weights and add it as a column in the data
  weights <- calculate_weights(combined_df %>% dplyr::select(-CVD,SEX))
  combined_df <- cbind(combined_df, weights)
  
  # Calculate the brier and AUC scores by gender
  scores_male <- c(scores_male, brier_score(combined_df, 1))
  scores_female <- c(scores_female, brier_score(combined_df, 2))
  
  #auc_male <- c(auc_male, auc(combined_df, 1))
  #auc_female <- c(auc_female, auc(combined_df, 2))
  
}

cov_mat <- function(sd1, sd2, rho) {
  
  cov_mat = matrix(0, nrow = 2, ncol = 2)
  cov_mat[1,1] <- sd1^2
  cov_mat[1,2] <- rho * sd1 * sd2
  cov_mat[2,1] <- rho * sd1 * sd2
  cov_mat[2,2] <- sd2^2 
  
  return(cov_mat)
  
}

library(MASS)
library(compositions)
# Simulate two datasets - male and female
# Male
set.seed(1)
nhanes_men <- nhanes_eligible %>% filter(AGE >= 30 & SEX == 1)
n <- nrow(nhanes_men)
n_sim <- 500
rho <- seq(0,0.90,0.1) # Varied Correlations
thresholds <- seq(120, 140, 5) # SYSBP thresholds
full_grid <- expand.grid(rho = rho, thresh = thresholds, iter = seq(1,n_sim))
full_grid$brier <- NA

mean_vals <- as.data.frame(t(apply(X = nhanes_men, MARGIN = 2, FUN = mean)))
sd_vals <- as.data.frame(t(apply(X = nhanes_men, MARGIN = 2, FUN = sd)))
for (i in 1:nrow(full_grid)) {

  cov_mat_chol <- cov_mat(sd_vals$HDLC, sd_vals$TOTCHOL, full_grid$rho[i])
  
  X <- as.data.frame(mvrnorm(n = n, mu = c(mean_vals$HDLC, mean_vals$TOTCHOL), 
                             Sigma = cov_mat_chol, 2))
  
  # Generate the cholesterol variables
  # Vary the covariance based on the values in the full_grid dataframe
  
  HDLC <- X[,1]
  TOTCHOL <- X[,2]
  
  # Generate the other continuous variables
  SYSBP <- rnorm(n, mean_vals$SYSBP, sd_vals$SYSBP)
  AGE <- round(runif(n, 30, 80))
  
  # Use the threshold in the full_grid dataframe to determine the value of bpmeds
  BP_MEDS <- as.numeric(SYSBP >= full_grid$thresh[i])
  
  # Generate smoking variable using the binomial distribution
  CURSMOKE <- rbinom(n, 1, mean_vals$CURSMOKE)
  SEX <- rep(1, n)
  
  data <- data.frame(HDLC = HDLC, TOTCHOL = TOTCHOL, SYSBP = SYSBP, AGE = AGE, 
                     BPMEDS = BP_MEDS, CURSMOKE = CURSMOKE, SEX = SEX)
  
  # Calculate probability of diabetes using the diabetes model
  # Have to use cursmoke as y variable to get the model matrix
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
    dplyr::select(HDLC, TOTCHOL, SYSBP, AGE, BPMEDS, CURSMOKE, 
                  SYSBP_UT, SYSBP_T, CVD, DIABETES, SEX)
  combined_df <- rbind(framingham, data) %>%
    mutate(S = case_when(is.na(CVD) ~ 0,
                         TRUE ~ 1))
  weights <- calculate_weights(combined_df %>% dplyr::select(-CVD, SEX))
  combined_df <- cbind(combined_df, weights)
  full_grid$brier[i] <- brier_score(combined_df, 1)
}

full_grid %>% group_by(rho, thresh) %>% summarise(bias = mean(brier - 0.0963176)) %>% arrange(abs(bias))

full_grid %>% filter(rho == 0.4 & thresh == 140) %>% summarise(mean(brier))


#store results as RDS file
saveRDS(full_grid, 'sim_results_men_1.rds')

# Female
# Same code as for men - only differences are the number of values & the mean/SD
set.seed(1)
nhanes_women <- nhanes_eligible %>% filter(AGE >= 30 & SEX == 2)
n <- nrow(nhanes_women)
n_sim <- 500
rho <- seq(0,0.90,0.1) # Varied Correlations
thresholds <- seq(120, 140, 5) # SYSBP thresholds
full_grid <- expand.grid(rho = rho, thresh = thresholds, iter = seq(1,n_sim))
full_grid$brier <- NA
sds <- c(sd(nhanes_women$HDLC, na.rm=TRUE), sd(nhanes_women$TOTCHOL, na.rm = TRUE))
#full_grid$auc <- NA
for (i in 1:nrow(full_grid)) {
  print(i)
  
  mu <- c(mean(nhanes_women$HDLC, na.rm = TRUE),
          mean(nhanes_women$TOTCHOL, na.rm = TRUE))
  
  cov_mat_chol <- cov_mat(sds[1], sds[2], full_grid$rho[i])
  
  X <- as.data.frame(mvrnorm(n = n, mu = mu, 
                             Sigma = cov_mat_chol, 2))
  
  # Generate the cholesterol variables
  # Vary the covariance based on the values in the full_grid dataframe
  
  HDLC <- X[,1]
  TOTCHOL <- X[,2]
  
  # Generate the other continuous variables
  SYSBP <- rnorm(n, mean(nhanes_women$SYSBP, na.rm = TRUE), sd(nhanes_women$SYSBP, na.rm = TRUE))
  AGE <- round(runif(n, 30, 80))
  
  # Use the threshold in the full_grid dataframe to determine the value of bp meds
  BP_MEDS <- as.numeric(SYSBP >= full_grid$thresh[i])
  
  # Generate smoking variable using the binomial distribution
  CURSMOKE <- rbinom(n, 1, mean(nhanes_women$CURSMOKE, na.rm = TRUE))
  SEX <- rep(2, n)
  
  data <- data.frame(HDLC = HDLC, TOTCHOL = TOTCHOL, SYSBP = SYSBP, AGE = AGE, 
                     BPMEDS = BP_MEDS, CURSMOKE = CURSMOKE, SEX = SEX)
  
  # Calculate probability of diabetes using the diabetes model
  # Have to use cursmoke as y variable to get the model matrix
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
  weights <- calculate_weights(combined_df %>% dplyr::select(-CVD, SEX))
  combined_df <- cbind(combined_df, weights)
  full_grid$brier[i] <- brier_score(combined_df, 2)
}

full_grid %>% group_by(rho, thresh) %>% summarise(bias = mean(brier - 0.0620515)) %>% arrange(abs(bias))

full_grid %>% filter(rho == 0.4 & thresh == 140)

saveRDS(full_grid, 'sim_results_women_1.rds')

# Store diabetes model
saveRDS(diabetes_mod, 'diabetes_mod.rds')
