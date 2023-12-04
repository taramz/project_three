library(riskCommunicator)
library(tidyverse)
library(tableone)
library(ggplot2)
library(tidyverse)

data("framingham")

# The Framingham data has been used to create models for cardiovascular risk.
# The variable selection and model below are designed to mimic the models used
# in the paper General Cardiovascular Risk Profile for Use in Primary Care 
# This paper is available (cvd_risk_profile.pdf) on Canvas.

framingham_df <- framingham %>% select(c(CVD, TIMECVD, SEX, TOTCHOL, AGE,
                                      SYSBP, DIABP, CURSMOKE, DIABETES, BPMEDS,
                                      HDLC, BMI))
framingham_df <- na.omit(framingham_df)

CreateTableOne(data=framingham_df, strata = c("SEX"))

# Get blood pressure based on whether or not on BPMEDS
framingham_df$SYSBP_UT <- ifelse(framingham_df$BPMEDS == 0, 
                                 framingham_df$SYSBP, 0)
framingham_df$SYSBP_T <- ifelse(framingham_df$BPMEDS == 1, 
                                framingham_df$SYSBP, 0)

# Looking at risk within 15 years - remove censored data
dim(framingham_df)
framingham_df <- framingham_df %>%
  filter(!(CVD == 0 & TIMECVD <= 365*15)) %>%
  select(-c(TIMECVD))
dim(framingham_df)

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


# The NHANES data here finds the same covariates among this national survey data
library(nhanesA)

# blood pressure, demographic, bmi, smoking, and hypertension info
bpx_2017 <- nhanes("BPX_J") %>% 
  select(SEQN, BPXSY1, BPXDI1) %>% 
  rename(SYSBP = BPXSY1, DIABP = BPXDI1)
demo_2017 <- nhanes("DEMO_J") %>% 
  select(SEQN, RIAGENDR, RIDAGEYR) %>% 
  rename(SEX = RIAGENDR, AGE = RIDAGEYR)
bmx_2017 <- nhanes("BMX_J") %>% 
  select(SEQN, BMXBMI) %>% 
  rename(BMI = BMXBMI)
smq_2017 <- nhanes("SMQ_J") %>%
  mutate(CURSMOKE = case_when(SMQ040 %in% c(1,2) ~ 1,
                              SMQ040 == 3 ~ 0, 
                              SMQ020 == 2 ~ 0)) %>%
  select(SEQN, CURSMOKE)
bpq_2017 <- nhanes("BPQ_J") %>% 
  mutate(BPMEDS = ifelse(BPQ050A == 1, 1, 0)) %>%
  select(SEQN, BPMEDS) 
tchol_2017 <- nhanes("TCHOL_J") %>% 
  select(SEQN, LBXTC) %>% 
  rename(TOTCHOL = LBXTC)
hdl_2017 <- nhanes("HDL_J") %>% 
  select(SEQN, LBDHDD) %>% 
  rename(HDLC = LBDHDD)
diq_2017 <- nhanes("DIQ_J") %>% 
  mutate(DIABETES = case_when(DIQ010 == 1 ~ 1, 
                              DIQ010 %in% c(2,3) ~ 0, 
                              TRUE ~ NA)) %>%
  select(SEQN, DIABETES) 

# Join data from different tables
df_2017 <- bpx_2017 %>%
  full_join(demo_2017, by = "SEQN") %>%
  full_join(bmx_2017, by = "SEQN") %>%
  full_join(hdl_2017, by = "SEQN") %>%
  full_join(smq_2017, by = "SEQN") %>%
  full_join(bpq_2017, by = "SEQN") %>%
  full_join(tchol_2017, by = "SEQN") %>%
  full_join(diq_2017, by = "SEQN")

CreateTableOne(data = df_2017, strata = c("SEX"))

# Correlation plot
library(corrplot)
head(df_2017)

numeric_nhanes <- df_2017 %>%
  select(SYSBP, DIABP, AGE, BMI, HDLC, TOTCHOL)
numeric_nhanes <- numeric_nhanes[complete.cases(numeric_nhanes),]
corrplot(cor(numeric_nhanes),
         type = 'lower',
         method = 'square',
         tl.col = 'black')

# Patterns of Missingness
# SYSBP and DIABP are missing together,
# HDLC and TOTCHOL are missing together
apply(X = is.na(df_2017), MARGIN = 2, FUN = sum)

sum(is.na(df_2017[is.na(df_2017$HDLC),]$TOTCHOL))

dim(df_2017)
dim(df_2017[complete.cases(df_2017),])

df_2017 %>%
  filter(AGE >= 28 & AGE <= 62) %>%
  group_by(SEX) %>%
  summarise(count = n())

complete_cases <- df_2017[complete.cases(df_2017),]
eligible_2017 <- complete_cases %>%
  filter(AGE >= 28 & AGE <= 62) %>%
  dplyr::select(c("SEX", "TOTCHOL", "AGE", "BMI", "SYSBP", "HDLC", 
                  "CURSMOKE", "BPMEDS", "DIABETES", "DIABP"))

# How to handle missing data
#6000 people who answered questionnaire who answered the questionnaire themself
#3000 people who did not complete questionnaire

# Histograms of Continuous Variables to show distributions
ggplot(eligible_2017) +
  geom_histogram(aes(x = HDLC), color = "black", fill = "white") +
  labs(title = 'Distribution of HDLC',
       x = 'HDLC',
       y = 'Count')

ggplot(eligible_2017) +
  geom_histogram(aes(x = TOTCHOL), color = "black", fill = "white") +
  labs(title = 'Distribution of Cholesterol',
       x = 'TOTCHOL',
       y = 'Count')

ggplot(eligible_2017) +
  geom_histogram(aes(x = AGE), color = "black", fill = "white") +
  labs(title = 'Distribution of AGE',
       x = 'AGE',
       y = 'Count')

ggplot(eligible_2017) +
  geom_histogram(aes(x = SYSBP), color = "black", fill = "white") +
  labs(title = 'Distribution of SYSBP',
       x = 'SYSBP',
       y = 'Count')






