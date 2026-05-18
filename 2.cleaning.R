# SET-UP =======================================================================
# Load libraries
library(dplyr)         # for data wrangling
library(tidyr)         # for data wrangling
library(here)          # for relative pathways
library(readxl)        # for writing Excel files
library(readr)         # writing rds files
library(writexl)       # writing xlsx files
library(forcats)       # for wrangling factors
library(vroom)         # for efficient reading of tabular data
library(fuzzyjoin)     # for joining comorbidities to PROMIS
library(data.table)    # for more memory-efficient data handling

# Load data
comorbidities <- vroom(here("1.data", "1.original", "comorbidities_v2.csv"))
complications <- vroom(here("1.data", "1.original", "complications_v2.csv"))
patient_data <- read_xlsx(here("1.data", "1.original", "patient-data_v2.xlsx")) # CSV's round fractional numbers to 0 and 1's
promis <- vroom(here("1.data", "1.original", "promis_v2.csv"))




# BASIC WRANGLING ==============================================================
# Change variables to lowercase
comorbidities <- comorbidities %>% rename_with(tolower)
complications <- complications %>% rename_with(tolower)
patient_data <- patient_data %>% rename_with(tolower)
promis <- promis %>% rename_with(tolower)

# Rename variables
comorbidities <- comorbidities %>% 
  rename(
    respondent_id = `claim id (randomized)`,
    days_from_procedure = `days from procedure (<0 is preop)`
  )

complications <- complications %>% 
  rename(
    respondent_id = `claim id (randomized)`,
    days_from_procedure = `days from procedure (<0 is preop)`
  )

patient_data <- patient_data %>% 
  rename(
    respondent_id = `claim id (randomized)`,
    SV1_SES = `socioeconomic status – rpl_theme1`,
    SV2_household = `household characteristics – rpl_theme2`,
    SV3_minority = `racial & ethnic minority status – rpl_theme3`,
    SV4_housing = `housing type & transportation – rpl_theme4`,
    SV5_sum = `overall summary ranking - rpl_theme1-4`,
    smoking_status = `smoking status`,
    age_cat = `age category`,
    bmi_cat = `bmi category`,
    procedure_loc = `procedure location`,
    procedure_year = `year`,
    gender = `gender`,
    relationship = `in relationship`
  )

promis <- promis %>% 
  rename(
    respondent_id = `claim id (randomized)`,
    instrument = `promis survey name`,
    tscore = `t-score`,
    days_from_procedure = `days from procedure (<0 is preop)`
  )

# Factorize PROMIS instruments
promis <- promis %>%
  mutate(
    instrument = fct_recode(
      as_factor(instrument),
      "anxiety" = "PROMIS item bank - emotional distress - anxiety - version 1.0",
      "depression" = "PROMIS item bank - emotional distress - depression - version 1.0",
      "pain_interference" = "PROMIS item bank - pain interference - version 1.1",
      "physical_function" = "PROMIS item bank - physical function - version 2.0",
      "self-efficacy" = "PROMIS item bank - self-efficacy for managing symptoms - version 1.0"
    )
  )

# For how many patients do we have PROMIS data?
length(unique(promis$respondent_id))
length(unique(patient_data$respondent_id))
length(unique(comorbidities$respondent_id))
length(unique(complications$respondent_id))




# COMPLICATIONS ================================================================
# Add complications to patient_data
# First we have to create a wide complications dataframe, if we want the wide structure
# of patient_data to stay intact
complication_flags <- complications %>%
  filter(days_from_procedure > 0) %>%
  select(respondent_id, complication) %>%
  mutate(value = "yes") %>%
  distinct() %>%                           # In case some complication occurs more than once
  pivot_wider(
    names_from = complication,
    values_from = value,
    values_fill = "no"
  )

# Merge the complication set with the patient_data
patient_data <- patient_data %>%
  left_join(complication_flags, by = "respondent_id") %>%
  rename(
    compl_other_or_unspecified = "Other or Unspecified",
    compl_mechanical = "Mechanical",
    compl_break_or_fracture = "Break or Fracture",
    compl_embolism = "Embolism",
    compl_infection = "Infection",
    compl_ssi = "Surgical Site Infection"
  )

# Change any remaining NA values to "no" as these represent patients that do not
# have a record in the complications dataframe
patient_data <- patient_data %>%
  mutate(across(
    c(compl_other_or_unspecified, compl_mechanical, compl_break_or_fracture, compl_embolism, compl_infection, compl_ssi), 
    ~ factor(replace_na(.x, "no"), levels = c("no", "yes"))
  ))

rm(complication_flags)

# Factorize remaining categorical/dichotomous variables in patient data
patient_data <- patient_data %>%
  mutate(
    smoking_status = factor(smoking_status),
    age_cat = factor(age_cat),
    bmi_cat = factor(bmi_cat),
    procedure_loc = factor(procedure_loc),
    procedure_year = factor(procedure_year),
    gender = factor(gender),
    relationship = factor(relationship)
  )

# Add patient_data to PROMIS
promis <- promis %>%
  left_join(patient_data, by = "respondent_id")




# COMORBIDITIES ================================================================
# Add comorbidities to PROMIS by matching them in time; find the nearest time measurement
# for each PROMIS measure in comorbidities and attach the comorbidity label

# For this, we will use the data.table package, which uses C instead of pure R. It also has
# built-in support for non-equitable joins like we have with our unmatching dates

# If we would do this in dplyr, R would crash: it would make promis rows x comorbidities rows comparisons.
# That's 615300 x 390687 = 240389711100 (240 billion)
# That would take terabytes of RAM...

# First, cut down the comorbidities df, which right now lists absent as well as present
# comorbidities. 
comorbidities_long <- comorbidities %>%
  pivot_longer(
    cols = -c(respondent_id, days_from_procedure),
    names_to = "comorbidity",
    values_to = "status"
  )

comorbidities <- comorbidities_long %>%
  filter(status == "Present") %>%
  select(respondent_id, days_from_procedure, comorbidity)

rm(comorbidities_long)

# Convert to data table for more efficient handling
promis_dt <- as.data.table(promis)
comorbidities_dt <- as.data.table(comorbidities)

# Rename shared time variable to avoid confusion
setnames(promis_dt, "days_from_procedure", "promis_time")
setnames(comorbidities_dt, "days_from_procedure", "comorbidities_time")

# Sort comorbidities table and define intervals
setorder(comorbidities_dt, respondent_id, comorbidities_time)
comorbidities_dt[, end := shift(comorbidities_time, type = "lead", fill = Inf), by = respondent_id]
comorbidities_dt[, start := comorbidities_time]

# Create a point interval in promis_dt
promis_dt[, `:=`(start = promis_time, end = promis_time)]

# Set keys
setkey(comorbidities_dt, respondent_id, start, end)
setkey(promis_dt, respondent_id, start, end)

# Join overlaps
merged_dt <- foverlaps(promis_dt, comorbidities_dt, type = "within", nomatch = NA)

# Create a variable denoting the time difference between matched PROMIS and comorbidities measurements
merged_dt[, time_difference := promis_time - comorbidities_time]

# For each respondent + PROMIS measurement + instrument, keep only the closest comorbidity match
# Create a temporary column for absolute time difference
merged_dt[, abs_time_diff := abs(time_difference)]

# Order by respondent, promis_time, instrument, then absolute time difference
setorder(merged_dt, respondent_id, promis_time, instrument, abs_time_diff)

# Keep only the first row per group (closest comorbidity)
merged_dt_filt <- merged_dt[, .SD[1], by = .(respondent_id, promis_time, instrument)]
removed_duplicates <- merged_dt[
  !merged_dt_filt, 
  on = .(respondent_id, promis_time, instrument, comorbidity, tscore)]

# Drop the temporary column
merged_dt_filt[, abs_time_diff := NULL]

# Make it a dataframe
merged_df <- as.data.frame(merged_dt_filt)

# Remove obsolete data objects
rm(merged_dt, promis_dt, comorbidities_dt)
# Factorize comorbidity
merged_df$comorbidity <- as_factor(merged_df$comorbidity)
# Remove obsolete variables
merged_df <- merged_df %>%
  select(-c(start, end, i.start, i.end))
# Reorder columns
merged_df <- merged_df %>%
  select(respondent_id, promis_time, comorbidities_time, time_difference, everything())

# The now filtered dataframe has the same row numbers as the promis dataframe so that seems
# alright. but let's make sure.
merged_df_sorted <- merged_df %>%
  arrange(respondent_id, promis_time, instrument)

promis_sorted <- promis %>%
  arrange(respondent_id, days_from_procedure, instrument)

all(merged_df_sorted$tscore == promis_sorted$tscore) # = TRUE! SO IT WORKS

# Lastly, some individuals have 0 on their PROMIS tscores. This is definitely an error; 
# there's no such thing. Let's make these NA
merged_df$tscore[merged_df$tscore == 0] <- NA




# SAVE =========================================================================
write.csv(comorbidities, here("1.data", "2.processed", "1.cleaned", "comorbidities_v2.csv"), row.names = F)
write.csv(complications, here("1.data", "2.processed", "1.cleaned", "complications_v2.csv"), row.names = F)
write_xlsx(patient_data, here("1.data", "2.processed", "1.cleaned", "patient-data_v2.xlsx"))
write.csv(promis, here("1.data", "2.processed", "1.cleaned", "promis_v2.csv"), row.names = F)
write_rds(merged_df, here("1.data", "2.processed", "1.cleaned", "merged.rds"))




