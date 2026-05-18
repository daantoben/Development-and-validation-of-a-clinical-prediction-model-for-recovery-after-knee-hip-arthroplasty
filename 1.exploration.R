# SET-UP =======================================================================
# Load libraries
library(dplyr)     # for data wrangling
library(tidyr)     # for data wrangling
library(here)      # for relative pathways
library(naniar)    # for missing data exploration
library(ggplot2)   # for data visualization
library(finalfit)  # for missing data exploration
library(knitr)     # for formatting tables
library(readxl)    # for writing Excel files
library(forcats)   # for wrangling factors
library(vroom)     # for efficient reading of tabular data
library(gtsummary) # for summary tables
library(patchwork) # for combining tables
library(stringr)   # for string manipulation

# Load data
comorbidities <- vroom(here("1.data", "1.original", "comorbidities.csv"))
complications <- vroom(here("1.data", "1.original", "complications.csv"))
patient_data <- vroom(here("1.data", "1.original", "patient-data.csv"))
promis <- vroom(here("1.data", "1.original", "promis.csv"))

# Change variables to lowercase
comorbidities <- comorbidities %>% rename_with(tolower)
complications <- complications %>% rename_with(tolower)
patient_data <- patient_data %>% rename_with(tolower)
promis <- promis %>% rename_with(tolower)

# Rename variables
comorbidities <- comorbidities %>% 
  rename(
    respondent_id = `claim id (randomized)`
    )

complications <- complications %>% 
  rename(
    respondent_id = `claim id (randomized)`,
    days_from_procedure = `days from procedure (<0 is preop)`
    )

patient_data <- patient_data %>% 
  rename(
    respondent_id = `claim id (randomized)`,
    zipcode = `zip code (first 3 digits)`,
    smoking_status = `smoking status`,
    age_cat = `age category`,
    bmi_cat = `bmi category`
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




# DATA WRANGLING ===============================================================
# For how many patients do we have PROMIS data?
length(unique(promis$respondent_id))
length(unique(patient_data$respondent_id))
length(unique(comorbidities$respondent_id))
length(unique(complications$respondent_id))

# How many PROMIS domains do we have?
table(promis$instrument)

# How many repeated promis measures per respondent?
desc <- promis %>%
  filter(!is.na(tscore)) %>%
  group_by(respondent_id, instrument) %>%
  summarise(n_tscore = n(), .groups = "drop")

desc %>% 
  ggplot(aes(x = n_tscore)) +
  geom_histogram(binwidth = 1, fill = "lightblue", colour = "black") +
  facet_wrap(~instrument, scales = "free_y") +
  xlab("Number of data points per respondent") +
  ylab("Count") +
  theme_bw()

# What is the length of follow-up?
promis %>%
  ggplot(aes(x = days_from_procedure, y = instrument)) +
  geom_boxplot() +
  xlab("Days from procedure") +
  ylab("PROMIS domain") +
  theme_bw()

promis %>%
  filter(days_from_procedure > 0 & instrument == "physical_function") %>%
  ggplot(aes(x = days_from_procedure, y = tscore, group = respondent_id)) +
  geom_line(alpha = 0.2) +
  theme_bw()

# What is the occurence of complications?
table(complications$complication)
table(complications$complication[complications$days_from_procedure > 0])

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
    other_or_unspecified = "Other or Unspecified",
    mechanical = "Mechanical",
    break_or_fracture = "Break or Fracture",
    embolism = "Embolism",
    infection = "Infection",
    ssi = "Surgical Site Infection"
  )

# Change any remaining NA values to "no" as these represent patients that do not
# have a record in the complications dataframe
patient_data <- patient_data %>%
  mutate(across(
    c(other_or_unspecified, mechanical, break_or_fracture, embolism, infection, ssi), 
    ~ factor(replace_na(.x, "no"), levels = c("no", "yes"))
  ))

# How are patient characteristics distributed?
patient_data %>% 
  tbl_summary(
    include = c(
      smoking_status, age_cat, bmi_cat, other_or_unspecified, mechanical,
      break_or_fracture, embolism, infection, ssi
      ),
    label = list(
      smoking_status ~ "Smoking behaviour",
      age_cat ~ "Age category",
      bmi_cat ~ "BMI category",
      other_or_unspecified = "complication: other",
      mechanical = "complication: mechanical",
      break_or_fracture = "complication: break/fracture",
      embolism = "complication: embolism",
      infection = "complication: infection",
      ssi = "complication: SSI"
    ),
    statistic = list(
      smoking_status ~ "{n} ({p}%)",
      age_cat ~ "{n} ({p}%)",
      bmi_cat ~ "{n} ({p}%)",
      other_or_unspecified ~ "{n} ({p}%)",
      mechanical ~ "{n} ({p}%)",
      break_or_fracture ~ "{n} ({p}%)",
      embolism ~ "{n} ({p}%)",
      infection ~ "{n} ({p}%)",
      ssi ~ "{n} ({p}%)"
    )
  ) %>%
  bold_labels()


# Comorbidities
comorbidity_counts <- comorbidities %>%
  # Only pivot comorbidity columns (those that start with "cci:" or "elix:")
  pivot_longer(
    cols = matches("^cci:|^elix:"),
    names_to = "full_name",
    values_to = "status"
  ) %>%
  filter(status == "Present") %>%
  mutate(
    index = if_else(str_starts(full_name, "cci:"), "cci", "elix"),
    comorbidity = full_name %>%
      str_remove("^cci: |^elix: ") %>%
      str_trim() %>%
      str_to_lower() %>%
      str_replace_all("[^a-z0-9]+", "_")
  ) %>%
  count(index, comorbidity, name = "count")

# Generate labels
comorbidity_labels <- c(
  "chf_congestive_heart_" = "Congestive heart failure",
  "copd_chronic_obstructive_pulmonary_" = "COPD",
  "diabetes" = "Diabetes",
  "diabetes_complications" = "Diabetes +",
  "mild_ld_liver_" = "Mild liver disease",
  "pvd_peripheral_vascular_" = "Peripheral vascular disease",
  "rd_renal_" = "Renal disease",
  "rheumatoid_disease" = "Rheumatoid disease",
  "alcohol_abuse" = "Alcohol abuse",
  "chronic_pulmonary_disease" = "CPD",
  "congestive_heart_failure" = "Congestive heart failure",
  "depression" = "Depression",
  "diabetes_complicated" = "Diabetes +",
  "diabetes_uncomplicated" = "Diabetes",
  "hypertension_uncomplicated" = "Hypertension",
  "liver disease" = "Liver disease",
  "obesity" = "Obesity",
  "peripheral_vascular disorders" = "Peripheral vascular",
  "renal_failure" = "Renal failure",
  "rheumatoid_arthritis_collagen_vascular" = "Rheumatoid arthritis"
)

comorbidity_counts <- comorbidity_counts %>%
  mutate(label = comorbidity_labels)

# Pie chart
comorbidity_counts %>%
  ggplot(aes(x = "", y= count, fill = label)) +
  geom_bar(stat = "identity", width = 1, colour = "white") +
  facet_wrap(~ index) +
  scale_fill_viridis(discrete = T) +
  coord_polar("y", start = 0) +
  theme_void()

# Lollipop chart
comorbidity_counts %>%
  ggplot(aes(x = label, y = count)) +
  geom_segment(aes(x = label, xend = label, y = 0, yend = count), 
               colour = "grey") +
  facet_wrap(~ index, ncol = 1) +
  geom_point(colour = "orange", size = 4) +
  theme_bw() +
  xlab("Comorbidities") +
  ylab("Count") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Chart PROMIS domains over time
colours <- c(
  "depression" = "blue",
  "pain_interference" = "red",
  "physical_function" = "darkgreen",
  "anxiety" = "purple",
  "self-efficacy" = "darkorange3"
)

promis %>%
  ggplot(aes(x = days_from_procedure, y = tscore, group = instrument, colour = instrument)) +
  geom_smooth(method = "gam", se = TRUE, show.legend = TRUE) +
  scale_colour_manual(values = colours) +
  theme_bw() + 
  coord_cartesian(ylim = c(37, 63)) +
  guides(colour = guide_legend(override.aes = list(fill = NA)))
