# SET-UP =======================================================================
# Load libraries
library(dplyr)        # for data wrangling
library(tidyr)        # for data wrangling
library(here)         # for relative pathways
library(forcats)      # for wrangling factors
library(caret)        # for splitting the data
library(readr)
library(lcmm)         # for LCGA
library(future)
library(furrr)
library(splines)

# Load data
df <- read_rds(here("1.data", "2.processed", "1.cleaned", "merged.rds"))
load(here("1.data", "3.final", "lcga_pf_spl.RData"))

# Prep data
# filter for promis instrument
df_pf <- df %>% filter(instrument == "physical_function")

# filter for respondents with preoperative PF scores
df_pf2 <- df_pf %>%
  group_by(respondent_id) %>%
  filter(sum(promis_time >= -365 & promis_time <= 21) >= 1) %>%
  ungroup()

df_pf3 <- df_pf2 %>%
  group_by(respondent_id) %>%
  filter(sum(promis_time >= 0 & promis_time <= 365) >= 2) %>%
  ungroup()

# Change respondent_ids to numeric, which hlme can handle
df_pf3$respondent_id <- as.character(df_pf3$respondent_id) %>% as.numeric(df_pf3$respondent_id)

# Specify time horizon
df_pf3 <- df_pf3 %>% filter(promis_time >=0 & promis_time <=365)


df_pf2 <- df_pf %>% filter(promis_time >=0 & promis_time <=365)


# LCGA =========================================================================
# Set seed
set.seed(2026)

# Step 1: Create a 1-class unconditional model ---------------------------------
lcga1_pf <- hlme(tscore ~ promis_time, subject = "respondent_id", ng = 1, data = df_pf)
lcga1_pf_spl <- hlme(tscore ~ ns(promis_time, df = 3), subject = "respondent_id", ng = 1, data = df_pf3)

# Step 2: Create k-class iterations based on the 1-class model -----------------
# Gridsearch iterates the algorithm over different starting values and subsequently
# selects the best model. The starting values are also saved, and will be used in 
# later bootstrapping.
# PF
lcga2_pf <- gridsearch(rep = 3, maxit = 70, minit= lcga1_pf,
                    m = hlme(tscore ~ promis_time, 
                             subject = "respondent_id", 
                             ng = 2, 
                             data = df_pf, 
                             mixture = ~ promis_time, 
                             B=lcga1_pf))

# PF - splines
lcga2_pf_spl <- gridsearch(rep = 3, maxit = 70, minit= lcga1_pf_spl,
                       m = hlme(tscore ~ ns(promis_time, df = 3),
                                subject = "respondent_id", 
                                ng = 2, 
                                data = df_pf, 
                                mixture = ~ ns(promis_time, df = 3), 
                                B=lcga1_pf_spl))


# PF
lcga3_pf <- gridsearch(rep = 3, maxit = 70, minit= lcga1_pf,
                    m = hlme(tscore ~ promis_time, 
                             subject = "respondent_id", 
                             ng = 3, 
                             data = df_pf, 
                             mixture = ~ promis_time,
                             random = ~1,
                             B=lcga3_pf))

# PF - splines
lcga3_pf_spl <- gridsearch(rep = 3, maxit = 70, minit= lcga1_pf_spl,
                           m = hlme(tscore ~ ns(promis_time, df = 3),
                                    subject = "respondent_id", 
                                    ng = 3, 
                                    data = df_pf3, 
                                    mixture = ~ ns(promis_time, df = 3), 
                                    B=lcga1_pf_spl))

# PF
lcga4_pf <- gridsearch(rep = 3, maxit = 70, minit= lcga1_pf,
                    m = hlme(tscore ~ promis_time, 
                             subject = "respondent_id", 
                             ng = 4, 
                             data = df_pf, 
                             mixture = ~ promis_time, 
                             B=lcga1_pf))

# PF - splines
lcga4_pf_spl <- gridsearch(rep = 3, maxit = 70, minit= lcga1_pf_spl,
                           m = hlme(tscore ~ ns(promis_time, df = 3),
                                    subject = "respondent_id", 
                                    ng = 4, 
                                    data = df_pf, 
                                    mixture = ~ ns(promis_time, df = 3), 
                                    B=lcga1_pf_spl))

# PF
lcga5_pf <- gridsearch(rep = 3, maxit = 70, minit= lcga1_pf,
                    m = hlme(tscore ~ promis_time, 
                             subject = "respondent_id", 
                             ng = 5, 
                             data = df_pf, 
                             mixture = ~ promis_time, 
                             B=lcga1_pf))

# PF - splines
lcga5_pf_spl <- gridsearch(rep = 3, maxit = 70, minit= lcga1_pf_spl,
                           m = hlme(tscore ~ ns(promis_time, df = 3),
                                    subject = "respondent_id", 
                                    ng = 5, 
                                    data = df_pf, 
                                    mixture = ~ ns(promis_time, df = 3), 
                                    B=lcga1_pf_spl))


lcga6_pf_spl <- gridsearch(rep = 3, maxit = 70, minit= lcga1_pf_spl,
                           m = hlme(tscore ~ ns(promis_time, df = 3),
                                    subject = "respondent_id", 
                                    ng = 6, 
                                    data = df_pf, 
                                    mixture = ~ ns(promis_time, df = 3), 
                                    B=lcga1_pf_spl))


# diagnostics
summaryplot(lcga1_pf_spl, lcga2_pf_spl, lcga3_pf_spl, lcga4_pf_spl, lcga5_pf_spl, 
            which = c("loglik", "entropy", "BIC", "AIC"))
summarytable(lcga1_pf, lcga2_pf, lcga3_pf, lcga4_pf, lcga5_pf)



summarytable(lcga3_pf_spl, lcga3_pf_spl2)
summaryplot(lcga3_pf_spl, lcga3_pf_spl2, which = c("entropy", "BIC"))


# Save objects
save(lcga1_pf_spl, lcga2_pf_spl, lcga3_pf_spl, lcga4_pf_spl, lcga5_pf_spl,
     file = here("1.data", "3.finalized", "lcga_pf_spl.RData"))

save(lcga6_pf_spl,
     file = here("1.data", "3.finalized", "lcga6_pf_spl.RData"))

# Compare models using Bootstrapped LRT
# In order to bootstrap, we must write a function which resamples the data with
# replacement B times, refit the cluster models on the bootstrapped data, and calculate the
# LRT statistic for each k versus k-1 model. Finally, the p-value is the proportion of B
# in which the test statistic exceeds the original LR.
model_km1 <- lcga4_pf_spl
model_k <- lcga5_pf_spl

bootstrap_LRT_parallel <- function(data, k, B, model_k, model_km1, workers = 8) {
  
  plan(multisession, workers = workers)
  
  results <- future_map(
    1:B,
    function(i, data, model_k, model_km1) {
      # Print progress
      cat("Bootstrap iteration:", i, "of", B, "\n")
      
      # Sample unique respondents to preserve longitudinal structure
      sampled_ids <- sample(unique(data$respondent_id), replace = TRUE)
      
      boot_data <- bind_rows(lapply(sampled_ids, function(id) {
        data[data$respondent_id == id, ]
      })
    )
      # Re-fit models on bootstrap sample
      model1_boot <- tryCatch({ # tryCatch makes sure a model not converging doesn;t break the rest of the loop
        hlme(tscore ~ ns(promis_time, df = 3),
             subject = "respondent_id", 
             B=model_km1$best, # use the starting values from the original models
             ng = k-1, 
             data = boot_data,
             mixture = ~ ns(promis_time, df = 3)
             )
      }, error = function(e) NULL)   # part of tryCatch, returns NULL if a model doesn't converge/work
      
      model2_boot <- tryCatch({
        hlme(tscore ~ ns(promis_time, df = 3), 
             subject = "respondent_id", 
             B=model_k$best,
             ng = k, 
             data = boot_data,
             mixture = ~ ns(promis_time, df = 3)
             )
      }, error = function(e) NULL)
      
      # Check if models converged
      if (!is.null(model1_boot) && !is.null(model2_boot) && 
          model1_boot$conv == 1 && model2_boot$conv == 1) {
        
        return(-2 * (model1_boot$loglik - model2_boot$loglik))
        
      } else {
        return(NA_real_)
      }
    },
    data = data,
    model_k = model_k,
    model_km1 = model_km1,
    .options = furrr_options(
      seed = TRUE,
      packages = c("dplyr", "lcmm", "splines")),
    .progress = TRUE
  )
  
  return(na.omit(unlist(results)))
}
  
BLRT_4v5_lcga <- bootstrap_LRT_parallel(
  data=df_pf, 
  k=5, 
  B=250,
  model_k = model_k,
  model_km1 = model_km1,
  workers = 8)



# Save data objects
save(lcga1_pf, lcga2_pf, lcga3_pf, lcga4_pf, lcga5_pf,
     file = here("1.data", "2.processed", "2.analyzed", "lcga_pf_v3.RData"))

save(BLRT_4v5_lcga, file = here("1.data", "3.finalized", "BLRT_4v5.RData"))








# VISUALIZATION ================================================================
# Extract class membership 
# (respondent ids are identical across lcga objects and also identical in their order)
classes <- data.frame(
  respondent_id = lcga1_pf_spl$pprob$respondent_id,
  lcga1_pf = factor(lcga1_pf_spl$pprob$class),
  lcga2_pf = factor(lcga2_pf_spl$pprob$class),
  lcga3_pf = factor(lcga3_pf_spl$pprob$class),
  lcga4_pf = factor(lcga4_pf_spl$pprob$class),
  lcga5_pf = factor(lcga5_pf_spl$pprob$class),
  lcga6_pf = factor(lcga6_pf_spl$pprob$class),
  lcga7_pf = factor(lcga7_pf_spl$pprob$class),
  maxprob2_pf = pmax(lcga2_pf_spl$pprob$prob1, lcga2_pf_spl$pprob$prob2),
  maxprob3_pf = pmax(lcga3_pf_spl$pprob$prob1, lcga3_pf_spl$pprob$prob2),
  maxprob4_pf = pmax(lcga4_pf_spl$pprob$prob1, lcga4_pf_spl$pprob$prob2),
  maxprob5_pf = pmax(lcga5_pf_spl$pprob$prob1, lcga5_pf_spl$pprob$prob2),
  maxprob6_pf = pmax(lcga6_pf_spl$pprob$prob1, lcga6_pf_spl$pprob$prob2),
  maxprob7_pf = pmax(lcga7_pf_spl$pprob$prob1, lcga7_pf_spl$pprob$prob2))
 
  
# Add new class variables to existing dataset
plot_df <- left_join(df_pf, 
                     classes, 
                     by = "respondent_id")

# Standardize LCGA labels
# Wrangle lcga factor variables
plot_df$lcga2_pf <- fct_recode(
  plot_df$lcga2_pf,
  "high recovery" = "1",
  "low recovery" = "2") %>%
  fct_relevel("low recovery", "high recovery")

plot_df$lcga3_pf <- fct_recode(
  plot_df$lcga3_pf,
  "low recovery" = "2",
  "moderate recovery" = "1",
  "high recovery" = "3") %>%
  fct_relevel("low recovery", "moderate recovery", "high recovery")

plot_df$lcga4_pf <- fct_recode(
  plot_df$lcga4_pf,
  "low recovery" = "3",
  "low-moderate recovery" = "2",
  "high-moderate recovery" = "1",
  "high recovery" = "4") %>%
  fct_relevel("low recovery", "low-moderate recovery", "high-moderate recovery", "high recovery")

plot_df$lcga5_pf <- fct_recode(
  plot_df$lcga5_pf,
  "low recovery" = "5",
  "low-moderate recovery" = "4",
  "steep growth recovery" = "1",
  "high-moderate recovery" = "2",
  "high recovery" = "3") %>%
  fct_relevel("low recovery", "low-moderate recovery", "steep growth recovery", "high-moderate recovery", "high recovery")


  
# Median lines plus loess smooth curves
colours <- c("high recovery" = "#FC4E07",
             "moderate recovery" = "#00AFBB",
             "low recovery" = "#E7B800",
             "low-moderate recovery" = "#52854C",
             "high-moderate recovery" = "#00AFBB",
             "steep growth recovery" = "#293352")

summary_df_2 <- plot_df %>%
  filter(!is.na(lcga2_pf)) %>%
  group_by(lcga2_pf, promis_time) %>%
  summarise(
    median_t = median(tscore, na.rm = TRUE),
    mean_t = mean(tscore, na.rm = TRUE),
    n = n(),
    sd_t = sd(tscore, na.rm = TRUE),
    se_t = sd_t / sqrt(n),
    .groups = "drop"
  )

summary_df_3 <- plot_df %>%
  filter(!is.na(lcga3_pf)) %>%
  group_by(lcga3_pf, promis_time) %>%
  summarise(
    median_t = median(tscore, na.rm = TRUE),
    mean_t = mean(tscore, na.rm = TRUE),
    n = n(),
    sd_t = sd(tscore, na.rm = TRUE),
    se_t = sd_t / sqrt(n),
    .groups = "drop"
  )

summary_df_4 <- plot_df %>%
  filter(!is.na(lcga4_pf)) %>%
  group_by(lcga4_pf, promis_time) %>%
  summarise(
    median_t = median(tscore, na.rm = TRUE),
    mean_t = mean(tscore, na.rm = TRUE),
    n = n(),
    sd_t = sd(tscore, na.rm = TRUE),
    se_t = sd_t / sqrt(n),
    .groups = "drop"
  )

summary_df_5 <- plot_df %>%
  filter(!is.na(lcga5_pf)) %>%
  group_by(lcga5_pf, promis_time) %>%
  summarise(
    median_t = median(tscore, na.rm = TRUE),
    mean_t = mean(tscore, na.rm = TRUE),
    n = n(),
    sd_t = sd(tscore, na.rm = TRUE),
    se_t = sd_t / sqrt(n),
    .groups = "drop"
  )

ggplot(summary_df_5, aes(x = promis_time, y = median_t, colour = lcga5_pf)) +
  geom_errorbar(aes(ymin = median_t - sd_t, ymax = median_t + sd_t), width = 0, alpha = 0.7) +
  geom_smooth(se = FALSE) +
  scale_color_manual(values = colours) +
  labs(
    x = "days from surgery",
    y = "T-score",
    colour = "Cluster"
  ) +
  coord_cartesian(ylim = c(18, 67)) +
  theme_bw() +
  theme(legend.position = "top")




# ==============================================================================
# CO-OCCURENCE MATRIX
# ==============================================================================
coocc_boot <- function(ids, index) {
  
  # -------------------------
  # 1 — bootstrap patients
  # -------------------------
  boot_ids <- ids[index]
  
  boot_ids_df <- tibble(
    respondent_id = boot_ids
  )
  
  boot_data <- boot_ids_df %>%
    left_join(model_df, by = "respondent_id", relationship = "many-to-many")
  
  # -------------------------
  # 2 — LCGA
  # -------------------------
  lcga <- tryCatch({
    hlme(
      tscore ~ splines::ns(promis_time, df = 3),
      subject = "respondent_id",
      ng = 3,
      data = boot_data,
      mixture = ~ splines::ns(promis_time, df = 3),
      B = init_B
    )
  }, error = function(e) return(NULL))
  
  if (is.null(lcga)) return(NULL)
  
  # -------------------------
  # 3 — class extraction 
  # -------------------------
  class_df <- lcga$pprob %>%
    dplyr::select(respondent_id, class)
  
  boot_data <- boot_data %>%
    left_join(class_df, by = "respondent_id")
  
  # -------------------------
  # 4 — one row per patient
  # -------------------------
  boot_unique <- boot_data %>%
    distinct(respondent_id, class)
  
  if (nrow(boot_unique) < 2) return(NULL)
  
  # -------------------------
  # 5 — co-occurrence matrix
  # -------------------------
  ids <- boot_unique$respondent_id
  classes <- boot_unique$class
  
  mat <- outer(classes, classes, FUN = "==") * 1
  rownames(mat) <- ids
  colnames(mat) <- ids
  
  return(mat)
}

# Employ
B <- 3
init_B <- lcga3_pf_spl$best
ids <- unique(model_df$respondent_id)

results <- vector("list", B)

for(b in seq_len(B)) {
  index <- sample(seq_along(ids), replace = TRUE)
  results[[b]] <- coocc_boot(ids, index)
}

# initialize matrix
co_mat <- matrix(0, length(ids), length(ids))
counts <- matrix(0, length(ids), length(ids))

rownames(co_mat) <- ids
colnames(co_mat) <- ids
rownames(counts) <- ids
colnames(counts) <- ids

for (res in results) {
  if (is.null(res)) next
  
  common <- intersect(rownames(res), ids)
  
  co_mat[common, common] <- co_mat[common, common] + res[common, common]
  counts[common, common] <- counts[common, common] + 1
}

# avoid division by zero
co_mat[counts > 0] <- co_mat[counts > 0] / counts[counts > 0]
co_mat[counts == 0] <- NA

# Heatmap
library(pheatmap)

pheatmap(
  co_mat,
  color = colorRampPalette(c("white", "blue"))(100),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = FALSE,
  show_colnames = FALSE
)
