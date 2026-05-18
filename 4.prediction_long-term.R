# Load libraries
library(dplyr)        # for data wrangling
library(tidyr)        # for data wrangling
library(here)         # for relative pathways
library(forcats)      # for wrangling factors
library(readr)
library(lme4)
library(lcmm)
library(ggplot2)
library(gtsummary)
library(yardstick)
library(pROC)
library(gt)
library(boot)
library(nnet)
library(fastDummies)

# Load data
df <- read_rds(here("1.data", "2.processed", "1.cleaned", "merged.rds"))
load(here("1.data", "3.final", "lcga_pf_spl.RData"))

df2 <- df %>% filter(is.na(preop_pf))
length(unique(df3$respondent_id))

df3 <- df2 %>%
  group_by(respondent_id) %>%
  filter(sum(promis_time >= 0 & promis_time <= 365) >= 2) %>%
  ungroup()

# ==============================================================================
# DATA WRANGLING
# ==============================================================================

# Create preoperative promis scores
df <- df %>%
  group_by(respondent_id) %>%
  mutate(
    preop_pf = mean(tscore[promis_time >= -365 & promis_time <= -21 & instrument == "physical_function"], na.rm = TRUE),
    preop_pi = mean(tscore[promis_time >= -365 & promis_time <= -21 & instrument == "pain_interference"], na.rm = TRUE),
    preop_dp = mean(tscore[promis_time >= -365 & promis_time <= -12 & instrument == "depression"], na.rm = TRUE)
  )

# Extract class membership and the max posterior probability
# (respondent ids are identical across lcga objects and also identical in their order)
classes <- data.frame(
  respondent_id = lcga1_pf_spl$pprob$respondent_id,
  lcga1_pf = factor(lcga1_pf_spl$pprob$class),
  lcga2_pf = factor(lcga2_pf_spl$pprob$class),
  lcga3_pf = factor(lcga3_pf_spl$pprob$class),
  lcga4_pf = factor(lcga4_pf_spl$pprob$class),
  lcga5_pf = factor(lcga5_pf_spl$pprob$class))

# Add cluster membership and posptob to dataframe
model_df <- left_join(df %>% filter(promis_time >= 0 & promis_time <=365 & instrument == "physical_function"), 
                     classes, 
                     by = "respondent_id")

# Wrangle lcga factor variables
model_df$lcga2_pf <- fct_recode(
  model_df$lcga2_pf,
  "high recovery" = "1",
  "low recovery" = "2") %>%
  fct_relevel("low recovery", "high recovery")

model_df$lcga3_pf <- fct_recode(
  model_df$lcga3_pf,
  "low recovery" = "2",
  "moderate recovery" = "1",
  "high recovery" = "3") %>%
  fct_relevel("low recovery", "moderate recovery", "high recovery")
 
model_df$lcga4_pf <- fct_recode(
  model_df$lcga4_pf,
  "low recovery" = "3",
  "low-moderate recovery" = "2",
  "high-moderate recovery" = "1",
  "high recovery" = "4") %>%
  fct_relevel("low recovery", "low-moderate recovery", "high-moderate recovery", "high recovery")

model_df$lcga5_pf <- fct_recode(
  model_df$lcga5_pf,
  "low recovery" = "5",
  "low-moderate recovery" = "4",
  "steep growth recovery" = "1",
  "high-moderate recovery" = "2",
  "high recovery" = "3") %>%
  fct_relevel("low recovery", "low-moderate recovery", "steep growth recovery", "high-moderate recovery", "high recovery")

# Set reference values for clusters  
model_df$lcga2_pf <- relevel(model_df$lcga2_pf, ref = "high recovery")
model_df$lcga3_pf <- relevel(model_df$lcga3_pf, ref = "moderate recovery")
model_df$lcga4_pf <- relevel(model_df$lcga4_pf, ref = "high-moderate recovery")
model_df$lcga5_pf <- relevel(model_df$lcga5_pf, ref = "high-moderate recovery")

# Set reference values for predictors
model_df$smoking_status <- relevel(model_df$smoking_status, ref = "Never")
model_df$age_cat <- relevel(model_df$age_cat, ref = "< 55")
model_df$bmi_cat <- factor(
  model_df$bmi_cat, 
  levels = c(
    "Normal or Underweight: < 25", "Overweight: 25 - 29.9", "Class I: 30 - 34.9",
    "Class II: 35 - 39.9", "Class III (Extreme Obesity): 40 or higher"))
model_df$bmi_cat <- relevel(model_df$bmi_cat, ref = "Normal or Underweight: < 25")

model_df$smoking_status[model_df$smoking_status == "Unknown"] <- NA
model_df$smoking_status <- droplevels(model_df$smoking_status)

# Collapse comorbidity levels into 5 classes and set a reference value (none)
model_df <- model_df %>%
  mutate(
    comorbidity = fct_na_value_to_level(comorbidity, "none or missing")
  )

model_df <- model_df %>%
  mutate(
    comorbidity = fct_collapse(comorbidity,
                               "high blood pressure" = "elix: hypertension, uncomplicated",
                               "cardiovascular disease" = c(
                                 "cci: chf (congestive heart)", 
                                 "cci: pvd (peripheral vascular)",
                                 "elix: congestive heart failure",
                                 "elix: peripheral vascular disorders"
                               ),
                               "diabetes" = c(
                                 "elix: diabetes, uncomplicated",
                                 "cci: diabetes",
                                 "elix: diabetes, complicated",
                                 "cci: diabetes + complications"
                               ),
                               "copd" = "cci: copd (chronic obstructive pulmonary)",
                               "other" = c(
                                 "elix: depression",
                                 "elix: chronic pulmonary disease",
                                 "cci: rheumatoid disease",
                                 "elix: renal failure",
                                 "elix: obesity",
                                 "elix: rheumatoid arthritis/collagen vascular",
                                 "cci: rd (renal)",
                                 "cci: mild ld (liver)",
                                 "elix: liver disease",
                                 "elix: alcohol abuse"
                               )))
model_df$comorbidity <- relevel(model_df$comorbidity, ref = "none or missing")


# Visualize clusters
colours <- c(
  "low recovery" = "#E7B800",
  "moderate recovery" = "#00AFBB",
  "high recovery" = "#FC4E07"
)

model_df %>%
  group_by(respondent_id) %>%
  ggplot(aes(x = promis_time, y = tscore, colour = comorbidity)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  theme_bw()

model_df %>%
  filter(!is.na(lcga3_pf)) %>%
  group_by(respondent_id) %>%
  ggplot(aes(x = promis_time, y = tscore, colour = lcga3_pf)) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 2) +
  geom_point(alpha = 0.1) +
  scale_colour_manual(values = colours) +
  coord_cartesian(ylim = c(17, 76)) +
  theme_bw() +
  theme(legend.position = "none")


# Create a wide dataframe. We can just slice, as our predictors and outcomes, lcga clusters, are static
model_df <- model_df %>%
  group_by(respondent_id) %>%
  slice(1) %>%
  ungroup()



# ==============================================================================
# MODELLING 
# ==============================================================================

# 2-cluster model
fit1 <- glm(
  lcga2_pf ~ age_cat + gender + bmi_cat + SV5_sum + relationship + smoking_status + comorbidity + procedure_loc + preop_pf + preop_pi, 
  data = model_df, family = "binomial", model = TRUE)

# 3-cluster model
fit2 <- multinom(
  lcga3_pf ~ age_cat + gender + bmi_cat + SV1_SES + SV2_household + SV3_minority + SV4_housing + relationship + 
  smoking_status + comorbidity + procedure_loc + preop_pf + preop_pi,
  data = model_df, model = TRUE)


save(fit2, file = here("1.data", "3.final", "fit2.RData"))

# 4-cluster model
fit3 <- multinom(
  lcga4_pf ~ age_cat + gender + bmi_cat + SV1_SES + SV2_household + SV3_minority + SV4_housing + relationship + 
  smoking_status + comorbidity + procedure_loc + preop_pf + preop_pi, 
  data = model_df, model = TRUE)


# 5-cluster model
fit4 <- multinom(
  lcga5_pf ~ age_cat + gender + bmi_cat + SV5_sum + relationship + smoking_status + comorbidity + procedure_loc + preop_pf + preop_pi, 
  data = model_df, model = TRUE)

# Diagnostics
probs_4 <- predict(fit4, type = "probs")
pred_class <- predict(fit4)
obs_4 <- fit4$model$lcga5_pf

table(Predicted = pred_class, Obs = obs_4)

df_eval <- data.frame(obs = obs_4, pred = pred_class)
classes <- levels(df_eval$obs)

accuracy_4 <- accuracy(df_eval, truth = obs, estimate = pred)
bal_accuracy_4 <- bal_accuracy(df_eval, truth = obs, estimate = pred)

sens_per_class_4 <- sapply(classes, function(cl) {
  mean(pred_class[obs_4 == cl] == cl)
})
spec_per_class_4 <- sapply(classes, function(cl) {
  mean(pred_class[obs_4 != cl] != cl)
})

auc_4 <- multiclass.roc(obs_4, probs_4)

obs_4_index <- as.numeric(obs_4)
p_true <- probs_4[cbind(1:nrow(probs_4), obs_4_index)]
log_loss_4 <- -mean(log(p_true))

df_onehot <- dummy_cols(data.frame(obs = obs_4), select_columns = "obs")
Y <- as.matrix(df_onehot[, grep("obs_", names(df_onehot))])
brier_4 <- mean(rowSums((Y - probs_4)^2))

tbl_regression(
  fit2,
  intercept = T,
  exponentiate = T,
  pvalue_fun = ~style_sigfig(., digits = 3)) %>%
  bold_labels()

save(fit2_tbl, file = here("1.data", "3.final", "fit2_tbl.RData"))

# Save model dataframe
write_rds(model_df, file = here("1.data", "3.final", "model_df.rds"))

# ==============================================================================
# INTERNAL VALIDATION
# ==============================================================================
set.seed(2026)


# 1 — Patient-level bootstrap

boot_ids <- sample(
  unique(model_df$respondent_id),
  size = length(unique(model_df$respondent_id)),
  replace = TRUE
)

# Assign new unique IDs for duplicates
boot_ids_df <- tibble(
  respondent_id = boot_ids,
  new_id = seq_along(boot_ids)
)

# Rebuild longitudinal data
boot_data <- boot_ids_df %>%
  left_join(model_df, by = "respondent_id", relationship = "many-to-many") %>%
  mutate(respondent_id2 = new_id) %>%  # LCGA subject IDs
  select(-new_id)


# 2 — Run LCGA on bootstrap

lcga <- hlme(
  tscore ~ promis_time,
  subject = "respondent_id2",
  ng = 3,
  data = boot_data,
  mixture = ~ promis_time,
  B = lcga3_pf$best
)

# Extract predicted classes
class_df <- lcga$pprob %>%
  select(respondent_id2, class)


# 3 — Merge class back

boot_data <- boot_data %>%
  left_join(class_df, by = "respondent_id2") %>%
  mutate(class = factor(class, levels = c(1,2,3)))


# 4 — Collapse to one row per patient

boot_data_wide <- boot_data %>%
  group_by(respondent_id2) %>%
  slice(1) %>%   # first row per patient
  ungroup()


# 5 — Fit multinomial model

m <- multinom(
  class ~ preop_pf + preop_pi + SV5_sum + comorbidity + 
    smoking_status + age_cat + bmi_cat + procedure_loc + 
    gender + relationship, 
  data = boot_data_wide, 
  model = TRUE,
  trace = FALSE
)


# 6 — Train AUC

probs_train <- predict(m, type = "probs")
obs_train <- m$model$class

# Ensure factor levels match
obs_train <- factor(obs_train, levels = c(1,2,3))
probs_train <- as.matrix(probs_train)

auc_train <- as.numeric(multiclass.roc(obs_train, probs_train)$auc)


# 7 — Test AUC on original dataset

# Collapse original dataset to one row per patient for prediction
model_df_wide <- model_df %>%
  group_by(respondent_id) %>%
  slice(1) %>%
  ungroup()

obs_test <- model_df_wide$lcga3_pf
valid <- !is.na(obs_test)
obs_test <- factor(obs_test[valid], levels = c(1,2,3))

probs_test <- predict(m, newdata = model_df_wide[valid, ], type = "probs")
probs_test <- as.matrix(probs_test)

auc_test <- as.numeric(multiclass.roc(obs_test, probs_test)$auc)


# 8 — Optimism estimate

optimism <- auc_train - auc_test

optimism




# ------------------------------------------------------------------
# Function to run one bootstrap iteration
# ------------------------------------------------------------------
auc_boot <- function(ids, index) {
  
  # -------------------------
  # 1 — bootstrap patients
  # -------------------------
  boot_ids <- ids[index]
  
  # Assign new unique IDs for duplicates
  boot_ids_df <- tibble(
    respondent_id = boot_ids,
    new_id = seq_along(boot_ids)
  )
  
  # Rebuild longitudinal dataset
  boot_data <- boot_ids_df %>%
    left_join(model_df, by = "respondent_id", relationship = "many-to-many") %>%
    mutate(respondent_id2 = new_id) %>%
    select(-new_id)
  
  # -------------------------
  # 2 — LCGA on bootstrap
  # -------------------------
  lcga <- tryCatch({
    hlme(
      tscore ~ ns(promis_time, df = 3),
      subject = "respondent_id2",
      ng = 3,
      data = boot_data,
      mixture = ~ ns(promis_time, df = 3),
      B = lcga3_pf_spl$best
    )
  }, error = function(e) return(NULL))
  
  if(is.null(lcga)) return(NA)  # skip iteration if LCGA fails
  
  class_df <- lcga$pprob %>%
    select(respondent_id2, class)
  
  boot_data <- boot_data %>%
    left_join(class_df, by = "respondent_id2") %>%
    mutate(class = factor(class, levels = c(1,2,3)))
  
  # -------------------------
  # 3 — Collapse to one row per patient
  # -------------------------
  boot_data_wide <- boot_data %>%
    group_by(respondent_id2) %>%
    slice(1) %>%
    ungroup()
  
  # Ensure at least 2 classes in sample
  if(length(unique(boot_data_wide$class)) < 2) return(NA)
  
  # -------------------------
  # 4 — Fit multinomial model
  # -------------------------
  m <- tryCatch({
    multinom(
      class ~ preop_pf + SV5_sum + comorbidity + 
        smoking_status + age_cat + bmi_cat + procedure_loc + 
        gender + relationship, 
      data = boot_data_wide,
      model = TRUE,
      trace = FALSE
    )
  }, error = function(e) return(NULL))
  
  if(is.null(m)) return(NA)
  
  # -------------------------
  # 5 — Train AUC
  # -------------------------
  probs_train <- predict(m, type = "probs")
  obs_train <- factor(m$model$class, levels = c(1,2,3))
  probs_train <- as.matrix(probs_train)
  
  if(length(unique(obs_train)) < 2) return(NA)
  
  auc_train <- as.numeric(multiclass.roc(obs_train, probs_train)$auc)
  
  # -------------------------
  # 6 — Test AUC on original dataset
  # -------------------------
  model_df_wide <- model_df %>%
    group_by(respondent_id) %>%
    slice(1) %>%
    ungroup()
  
  obs_test <- factor(model_df_wide$lcga3_pf, levels = c(1,2,3))
  valid <- !is.na(obs_test)
  obs_test <- obs_test[valid]
  
  probs_test <- predict(m, newdata = model_df_wide[valid, ], type = "probs")
  probs_test <- as.matrix(probs_test)
  
  auc_test <- as.numeric(multiclass.roc(obs_test, probs_test)$auc)
  
  # -------------------------
  # 7 — Optimism
  # -------------------------
  c(train = auc_train, test = auc_test)
}


# ------------------------------------------------------------------
# Bootstrap IDs
# ------------------------------------------------------------------
ids <- unique(model_df$respondent_id)

# Run bootstrap
set.seed(2026)
auc_out <- boot(
  data = ids,
  statistic = auc_boot,
  R = 500  
)

# Inspect results
mean(auc_out$t[, 1], na.rm = TRUE)
mean(auc_out$t[, 2],  na.rm = TRUE)
mean(auc_out$t[, 1], na.rm = TRUE) - mean(auc_out$t[, 2],  na.rm = TRUE)

# Save results
save(auc_out, file = here("1.data", "3.final", "auc_out.RData"))
write_rds(model_df, file = here("1.data", "3.final", "model_df.rds"))




# ------------------------------------------------------------------
# BLRT - BLRT vectors were calculated in MyDRE
# ------------------------------------------------------------------
load(here("1.data", "3.final", "BLRT_1v2.RData"))
load(here("1.data", "3.final", "BLRT_2v3.RData"))
load(here("1.data", "3.final", "BLRT_3v4.RData"))
load(here("1.data", "3.final", "BLRT_4v5.RData"))

# Calculate Likelihood Ratio Tests to contrast with the BLRT
LR_2class <- -2 * (lcga1_pf_spl$loglik - lcga2_pf_spl$loglik)
LR_3class <- -2 * (lcga2_pf_spl$loglik - lcga3_pf_spl$loglik)
LR_4class <- -2 * (lcga3_pf_spl$loglik - lcga4_pf_spl$loglik)
LR_5class <- -2 * (lcga4_pf_spl$loglik - lcga5_pf_spl$loglik)

# Compute the average times the BLRT were greater than the LRT
mean(BLRT_1v2_lcga >= LR_2class)
mean(BLRT_2v3_lcga >= LR_3class)
mean(BLRT_3v4_lcga >= LR_4class)
mean(BLRT_4v5_lcga >= LR_5class)





