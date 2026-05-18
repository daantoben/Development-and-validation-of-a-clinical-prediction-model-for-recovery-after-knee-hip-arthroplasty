# Load libraries
library(gtsummary)
library(forcats)
library(dplyr)
library(here)

# Load data
model_df <- read_rds(here("1.data", "3.final", "model_df.rds"))

# Filter NA values
model_df <- model_df %>%
  filter(is.na(preop_pf) == FALSE & is.na(preop_pi) == FALSE & is.na(lcga3_pf) == FALSE &
           is.na(SV5_sum) == FALSE)

# Relabel last age category
model_df$age_cat <- fct_recode(model_df$age_cat, "80-88" = ">=80")

# Relabel lcga category
model_df$lcga3_pf <- fct_recode(
  model_df$lcga3_pf,
  "low recovery" = "1",
  "moderate recovery" = "3",
  "high recovery" = "2") %>%
  fct_relevel("low recovery", "moderate recovery", "high recovery")

# Construct table
tbl1 <- model_df %>% 
  tbl_summary(
    include = c(
      age_cat, gender, relationship, SV1_SES, SV2_household, SV3_minority, SV4_housing, 
      smoking_status, bmi_cat, comorbidity, procedure_loc, procedure_year, preop_pf, preop_pi,
    ), missing = "ifany",
    label = list(
      gender ~ "Gender",
      age_cat ~ "Age",
      relationship ~ "Relationship status",
      SV1_SES ~ "Social vulnerability - SES",
      SV2_household ~ "Social vulnerability - household",
      SV3_minority ~ "Social vulnerability - minority",
      SV4_housing ~ "Social vulnerability - housing",
      smoking_status ~ "Smoking habit",
      bmi_cat ~ "bmi",
      comorbidity ~ "Comorbidity",
      procedure_loc ~ "Procedure location",
      procedure_year ~ "Procedure year",
      preop_pf ~ "Preoperative physical functioning",
      preop_pi ~ "Preoperative pain interference"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_categorical() ~ 1)
  ) %>%
  bold_labels()

tbl1
summary(model_df$relationship)
