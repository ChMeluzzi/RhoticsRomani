# =============================================================================
# RHOTIC VARIATION IN ROMANI — BAYESIAN INFERENTIAL ANALYSIS
# Series 2: DurStand (standardized duration, continuous)
# REVISED: RhoticType included as covariate in all models
#          Variety labels corrected: Burgudži, Kalajdži, Kalderaš
# =============================================================================
# Dependent variable : DurStand (z-score standardized duration)
# Covariate          : RhoticType_f (TAP vs. TRILL+) — in all models
# Random intercept   : (1 | Word)
# Dataset            : TAP + TRILL+ only (APPROX excluded for consistency
#                      with Series 1; PreEU only for Models 2a and 2b)
# Framework          : Bayesian, rstanarm::stan_lmer()
# Prior              : Normal(0, 1) on z-score scale
#
# MODEL 1 (stepping stone):
#   DurStand ~ RhoticType_f + PhonCont_f + Variety_f + (1|Word)
#   Full dataset | All varieties pooled
#   Reference: TAP | Intervocalic | Burgudži
#
# MODEL 2a:
#   DurStand ~ RhoticType_f + PhonCont_f + OriginUltraAggregate_f + (1|Word)
#   Full dataset | Fit separately per variety
#   Reference: TAP | Intervocalic | PreEU
#
# MODEL 2b:
#   DurStand ~ RhoticType_f + PhonCont_f + FormerRhotic_f + (1|Word)
#   PreEU tokens only | Fit separately per variety
#   Reference: TAP | Intervocalic | Vibrant
#
# Model comparison: LOO-CV (Vehtari et al., 2017)
# =============================================================================


# -----------------------------------------------------------------------------
# 0. PACKAGES
# -----------------------------------------------------------------------------
# Uncomment at first run:
# install.packages(c("rstanarm", "bayesplot", "tidybayes",
#                    "tidyverse", "readxl", "patchwork", "scales"))

library(tidyverse)
library(readxl)
library(rstanarm)
library(bayesplot)
library(tidybayes)
library(patchwork)
library(scales)


# -----------------------------------------------------------------------------
# 1. DATASET
# -----------------------------------------------------------------------------

setwd("C:/Users/Fin/Documents/LAVORO/Pubblicazioni/Romani_MeliMeluzzi_AISV26/New Scripts R")

df_raw <- read_xlsx("romani_unified_new_Rstand_FINAL.xlsx") %>%
  mutate(DurStand = as.numeric(DurStand))

# Base dataset: TAP + TRILL+ only, complete cases
df_base <- df_raw %>%
  filter(
    !is.na(DurStand),
    RhoticType %in% c("TAP", "TRILL+"),
    !is.na(PhonCont),
    !is.na(OriginUltraAggregate),
    !is.na(FormerRhotic),
    !is.na(Word)
  ) %>%
  mutate(
    # Correct variety labels
    Variety = recode(Variety,
                     "Burgudzi" = "Burgudži",
                     "Kalajdzi" = "Kalajdži",
                     "Kaldaras" = "Kalderaš"),
    RhoticType_f           = factor(RhoticType,
                                    levels = c("TAP", "TRILL+")),
    PhonCont_f             = factor(PhonCont,
                                    levels = c("Initial", "Intervocalic", "Final")),
    Variety_f              = factor(Variety,
                                    levels = c("Burgudži", "Kalajdži", "Kalderaš")),
    OriginUltraAggregate_f = factor(OriginUltraAggregate,
                                    levels = c("PreEU", "EU")),
    FormerRhotic_f         = factor(FormerRhotic,
                                    levels = c("Vibrant", "Retroflex")),
    Word                   = factor(Word)
  )

# PreEU subset for Model 2b
df_preeu <- df_base %>% filter(OriginUltraAggregate == "PreEU")

varieties <- c("Burgudži", "Kalajdži", "Kalderaš")

cat("=== DATASET ===\n")
cat("Full (TAP + TRILL+, excl. NA):", nrow(df_base), "tokens\n")
cat("PreEU only                    :", nrow(df_preeu), "tokens\n\n")
for (v in varieties) {
  cat(v, "— full:", nrow(df_base[df_base$Variety == v, ]),
      "| PreEU:", nrow(df_preeu[df_preeu$Variety == v, ]), "\n")
}


# -----------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# -----------------------------------------------------------------------------

# Posterior summary table (fixed effects)
posterior_summary_lm <- function(model, variety_name) {
  post <- as.data.frame(model) %>%
    select(starts_with("(Intercept)") | starts_with("PhonCont") |
             starts_with("Variety") | starts_with("FormerRhotic") |
             starts_with("RhoticType") | starts_with("OriginUltra"))

  post %>%
    pivot_longer(everything(), names_to = "Term", values_to = "value") %>%
    group_by(Term) %>%
    summarise(
      Median  = round(median(value), 3),
      SD      = round(sd(value), 3),
      CI_low  = round(quantile(value, 0.025), 3),
      CI_high = round(quantile(value, 0.975), 3),
      p_dir   = round(mean(value > 0), 3),
      .groups = "drop"
    ) %>%
    mutate(
      Variety = variety_name,
      sig = case_when(
        pmin(p_dir, 1 - p_dir) < 0.025 ~ "**",
        pmin(p_dir, 1 - p_dir) < 0.05  ~ "*",
        pmin(p_dir, 1 - p_dir) < 0.1   ~ ".",
        TRUE                            ~ ""
      )
    ) %>%
    select(Variety, Term, Median, SD, CI_low, CI_high, p_dir, sig)
}

# Predicted means (marginal over Word)
get_pred_means <- function(model, grid, df_v, variety_name) {
  grid_full <- grid %>%
    mutate(Word = levels(df_v$Word)[1])
  post_pred <- posterior_epred(model,
                               newdata = grid_full,
                               re.form = NA)
  grid %>%
    mutate(
      Variety = variety_name,
      EMM     = apply(post_pred, 2, median),
      CI_low  = apply(post_pred, 2, quantile, 0.025),
      CI_high = apply(post_pred, 2, quantile, 0.975)
    )
}

# Standard plot theme
theme_romani <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      strip.text         = element_text(face = "bold", size = 11),
      plot.title         = element_text(face = "bold"),
      plot.subtitle      = element_text(colour = "grey40", size = 9),
      legend.position    = "bottom",
      axis.text.x        = element_text(angle = 20, hjust = 1)
    )
}

palette_phonCont <- c("Intervocalic" = "#4E9AF1",
                      "Initial"      = "#E86F51",
                      "Final"        = "#6ABF69")
palette_variety  <- c("Burgudži" = "#4E9AF1",
                      "Kalajdži" = "#E86F51",
                      "Kalderaš" = "#6ABF69")
palette_rhotic   <- c("TAP"    = "#4E9AF1",
                      "TRILL+" = "#E86F51")


##DESCRIPTIVE ANALYSIS##
# Tabella mediane e SD
tab_desc <- df_viol %>%
  group_by(Variety_f, PhonCont_f, RhoticType_f) %>%
  summarise(
    N      = n(),
    Median = round(median(DurStand), 3),
    SD     = round(sd(DurStand), 3),
    .groups = "drop"
  )

print(tab_desc, n = Inf)
write.csv(tab_desc, "Tab_Descriptive_DurStand.csv", row.names = FALSE)
cat("Salvato: Tab_Descriptive_DurStand.csv\n")

# Grafico: varietà affiancate, violini TAP vs TRILL+ per contesto
stat_med <- df_viol %>%
  group_by(Variety_f, PhonCont_f, RhoticType_f) %>%
  summarise(med = median(DurStand), .groups = "drop")

p_finale <- df_viol %>%
  ggplot(aes(x = RhoticType_f, y = DurStand, fill = RhoticType_f)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = "grey40") +
  geom_text(data = stat_med,
            aes(x = RhoticType_f, y = med, label = round(med, 2)),
            vjust = -0.8, size = 2.5,
            fontface = "bold", colour = "grey10") +
  facet_grid(Variety_f ~ PhonCont_f) +
  scale_fill_manual(values = c("TAP" = "#4E9AF1", "TRILL+" = "#E86F51"),
                    name = "Rhotic type") +
  labs(
    title    = "Standardized duration by rhotic type, variety and phonological context",
    subtitle = "Violin + annotated median | TAP vs. TRILL+",
    x        = "Rhotic type",
    y        = "Standardized duration (z-score)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 10, colour = "grey40"),
    strip.text         = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_blank(),
    legend.position    = "bottom"
  )

ggsave("Fig2_Descriptive_DurStand_Violini.jpg", p_finale,
       width = 12, height = 10, dpi = 300)
cat("Salvato: Fig2_Descriptive_DurStand_Violini.jpg\n")

# =============================================================================
# 3. MODEL 1 — Stepping stone
#    DurStand ~ RhoticType_f + PhonCont_f + Variety_f + (1|Word)
#    Full dataset, all varieties pooled
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("MODEL 1 — Stepping stone (all varieties)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

m1 <- stan_lmer(
  DurStand ~ RhoticType_f + PhonCont_f + Variety_f + (1 | Word),
  data            = df_base,
  prior           = normal(0, 1),
  prior_intercept = normal(0, 1),
  seed            = 42,
  iter            = 4000,
  chains          = 4,
  cores           = 4,
  refresh         = 500
)

cat("\n--- Posterior summary: Model 1 ---\n")
print(summary(m1,
              pars = c("(Intercept)",
                       "RhoticType_fTRILL+",
                       "PhonCont_fInitial", "PhonCont_fFinal",
                       "Variety_fKalajdži", "Variety_fKalderaš"),
              probs = c(0.025, 0.975),
              digits = 3))

# Coefficient table
tab_m1 <- posterior_summary_lm(m1, "All varieties")
cat("\n--- Coefficient table: Model 1 ---\n")
print(tab_m1)
write.csv(tab_m1, "Tab_Model1_Bayes_DurStand_coefficients.csv",
          row.names = FALSE)

# LOO-CV
loo_m1 <- loo(m1)
cat("\n--- LOO-CV: Model 1 ---\n")
print(loo_m1)

# PP check
p_pp_m1 <- pp_check(m1, nreps = 100) +
  labs(title = "Model 1 DurStand — Posterior predictive check") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
ggsave("Fig_Model1_Bayes_DurStand_PPcheck.pdf", p_pp_m1,
       width = 7, height = 5, dpi = 300)
ggsave("Fig_Model1_Bayes_DurStand_PPcheck.jpg", p_pp_m1,
       width = 7, height = 5, dpi = 300)

# Predicted means: RhoticType x PhonCont x Variety
grid_m1 <- expand.grid(
  RhoticType_f = levels(df_base$RhoticType_f),
  PhonCont_f   = levels(df_base$PhonCont_f),
  Variety_f    = levels(df_base$Variety_f)
)
pred_m1 <- get_pred_means(m1, grid_m1, df_base, "All") %>%
  mutate(
    Variety_f    = factor(Variety_f, levels = varieties),
    PhonCont_f   = factor(PhonCont_f,
                          levels = c("Initial", "Intervocalic", "Final")),
    RhoticType_f = factor(RhoticType_f, levels = c("TAP", "TRILL+"))
  )

p_pred_m1 <- pred_m1 %>%
  ggplot(aes(x = PhonCont_f, y = EMM,
             colour = RhoticType_f, group = RhoticType_f)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.15, linewidth = 0.8,
                position = position_dodge(width = 0.3)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.3)) +
  facet_wrap(~ Variety_f, nrow = 1) +
  scale_colour_manual(values = palette_rhotic, name = "Rhotic type") +
  labs(
    title    = "Model 1 (Bayesian) — Posterior predicted means of standardized duration",
    subtitle = "stan_lmer | RhoticType + PhonCont + Variety | Prior: Normal(0,1) | Random intercept: Word\nReference: TAP | Intervocalic | Burgudži | Dashed line = grand mean (z = 0)",
    x        = "Phonological context",
    y        = "Posterior predicted mean (z-score)"
  ) +
  theme_romani()


p_pred_m1 <- pred_m1 %>%
  mutate(PhonCont_f = factor(PhonCont_f,
                             levels = c("Initial", "Intervocalic", "Final"))) %>%
  ggplot(aes(x = PhonCont_f, y = EMM,
             colour = RhoticType_f, group = RhoticType_f)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.15, linewidth = 0.8,
                position = position_dodge(width = 0.3)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.3)) +
  facet_wrap(~ Variety_f, nrow = 1) +
  scale_colour_manual(values = palette_rhotic, name = "Rhotic type") +
  labs(
    title    = "Model 1 (Bayesian) — Posterior predicted means of standardized duration",
    subtitle = "stan_lmer | RhoticType + PhonCont + Variety | Prior: Normal(0,1) | Random intercept: Word\nReference: TAP | Intervocalic | Burgudži | Dashed line = grand mean (z = 0)",
    x        = "Phonological context",
    y        = "Posterior predicted mean (z-score)"
  ) +
  theme_romani()

print(p_pred_m1)
ggsave("Fig_Model1_Bayes_DurStand_PredMeans.pdf", p_pred_m1,
       width = 12, height = 6, dpi = 300)
ggsave("Fig_Model1_Bayes_DurStand_PredMeans.jpg", p_pred_m1,
       width = 12, height = 6, dpi = 300)
cat("Model 1 complete.\n")


# =============================================================================
# 4. MODEL 2a — OriginUltraAggregate
#    DurStand ~ RhoticType_f + PhonCont_f + OriginUltraAggregate_f + (1|Word)
#    Full dataset | Fit separately per variety
# =============================================================================

results_2a <- list()

for (v in varieties) {

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("MODEL 2a —", v, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  df_v <- df_base %>% filter(Variety == v)

  cat("N =", nrow(df_v),
      "| TAP:", sum(df_v$RhoticType == "TAP"),
      "| TRILL+:", sum(df_v$RhoticType == "TRILL+"), "\n")
  cat("EU:", sum(df_v$OriginUltraAggregate == "EU"),
      "| PreEU:", sum(df_v$OriginUltraAggregate == "PreEU"), "\n\n")

  m_2a <- stan_lmer(
    DurStand ~ RhoticType_f + PhonCont_f + OriginUltraAggregate_f + (1 | Word),
    data            = df_v,
    prior           = normal(0, 1),
    prior_intercept = normal(0, 1),
    seed            = 42,
    iter            = 4000,
    chains          = 4,
    cores           = 4,
    refresh         = 500
  )

  # LOO
  loo_2a <- loo(m_2a)
  cat("\n--- LOO-CV:", v, "---\n")
  print(loo_2a)

  # Coefficient table
  tab <- posterior_summary_lm(m_2a, v)
  cat("\n--- Coefficient table:", v, "---\n")
  print(tab)
  write.csv(tab,
            paste0("Tab_Model2a_Bayes_DurStand_", v, "_coefficients.csv"),
            row.names = FALSE)

  # Predicted means: RhoticType x PhonCont (marginal over OriginUltraAggregate)
  grid_2a <- expand.grid(
    RhoticType_f           = levels(df_v$RhoticType_f),
    PhonCont_f             = levels(df_v$PhonCont_f),
    OriginUltraAggregate_f = levels(df_v$OriginUltraAggregate_f)
  )
  pred <- get_pred_means(m_2a, grid_2a, df_v, v)

  results_2a[[v]] <- list(model = m_2a, tab = tab,
                          loo = loo_2a, pred = pred)
  cat("\nModel 2a", v, "complete.\n")
}

# LOO comparison: Model 1 vs Model 2a per variety
cat("\n=== LOO-CV COMPARISON: skipped — dataset size mismatch ===\n")
cat("NOTE: LOO comparison between Model 1 and Model 2a not performed\n")
cat("due to unequal number of observations across models.\n")
cat("Model fit assessed via individual LOO-CV for each model.\n")

# PP checks
pp_plots_2a <- list()
for (v in varieties) {
  pp_plots_2a[[v]] <- pp_check(results_2a[[v]]$model, nreps = 100) +
    labs(title = v) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
}
p_pp_2a <- wrap_plots(pp_plots_2a, nrow = 1) +
  plot_annotation(
    title    = "Posterior predictive checks — Model 2a DurStand (Bayesian)",
    subtitle = "Observed (dark) vs. 100 posterior predictive draws (light)",
    theme    = theme(plot.title = element_text(face = "bold"))
  )
ggsave("Fig_Model2a_Bayes_DurStand_PPcheck.pdf", p_pp_2a,
       width = 12, height = 5, dpi = 300)
ggsave("Fig_Model2a_Bayes_DurStand_PPcheck.jpg", p_pp_2a,
       width = 12, height = 5, dpi = 300)

# Predicted means plot: RhoticType x PhonCont, facet Variety x Origin
pred_all_2a <- bind_rows(lapply(results_2a, `[[`, "pred")) %>%
  mutate(
    Variety                = factor(Variety, levels = varieties),
    PhonCont_f             = factor(PhonCont_f,
                                    levels = c("Initial", "Intervocalic", "Final")),
    OriginUltraAggregate_f = factor(OriginUltraAggregate_f,
                                    levels = c("PreEU", "EU")),
    RhoticType_f           = factor(RhoticType_f, levels = c("TAP", "TRILL+"))
  )

p_pred_2a <- pred_all_2a %>%
  ggplot(aes(x = PhonCont_f, y = EMM,
             colour = RhoticType_f, group = RhoticType_f)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.15, linewidth = 0.8,
                position = position_dodge(width = 0.3)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.3)) +
  facet_grid(OriginUltraAggregate_f ~ Variety) +
  scale_colour_manual(values = palette_rhotic, name = "Rhotic type") +
  labs(
    title    = "Model 2a (Bayesian) — Posterior predicted means of standardized duration",
    subtitle = "stan_lmer | RhoticType + PhonCont + OriginUltraAggregate | Prior: Normal(0,1) | Random intercept: Word\nRows: PreEU / EU | Columns: Variety | Dashed line = grand mean (z = 0)",
    x        = "Phonological context",
    y        = "Posterior predicted mean (z-score)"
  ) +
  theme_romani()

print(p_pred_2a)
ggsave("Fig_Model2a_Bayes_DurStand_PredMeans.pdf", p_pred_2a,
       width = 12, height = 8, dpi = 300)
ggsave("Fig_Model2a_Bayes_DurStand_PredMeans.jpg", p_pred_2a,
       width = 12, height = 8, dpi = 300)
cat("\nModel 2a complete.\n")


# =============================================================================
# 5. MODEL 2b — FormerRhotic
#    DurStand ~ RhoticType_f + PhonCont_f + FormerRhotic_f + (1|Word)
#    PreEU tokens only | Fit separately per variety
# =============================================================================

results_2b <- list()

for (v in varieties) {

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("MODEL 2b —", v, "(PreEU only)\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  df_v <- df_preeu %>% filter(Variety == v)

  cat("N =", nrow(df_v),
      "| TAP:", sum(df_v$RhoticType == "TAP"),
      "| TRILL+:", sum(df_v$RhoticType == "TRILL+"), "\n")
  cat("Vibrant:", sum(df_v$FormerRhotic == "Vibrant"),
      "| Retroflex:", sum(df_v$FormerRhotic == "Retroflex"), "\n\n")

  if (v == "Burgudži") {
    cat("NOTE: Retroflex x Final = 1 token only — interpret with caution\n\n")
  }

  m_2b <- stan_lmer(
    DurStand ~ RhoticType_f + PhonCont_f + FormerRhotic_f + (1 | Word),
    data            = df_v,
    prior           = normal(0, 1),
    prior_intercept = normal(0, 1),
    seed            = 42,
    iter            = 4000,
    chains          = 4,
    cores           = 4,
    refresh         = 500
  )

  # LOO
  loo_2b <- loo(m_2b)
  cat("\n--- LOO-CV:", v, "---\n")
  print(loo_2b)

  # Coefficient table
  tab <- posterior_summary_lm(m_2b, v)
  cat("\n--- Coefficient table:", v, "---\n")
  print(tab)
  write.csv(tab,
            paste0("Tab_Model2b_Bayes_DurStand_", v, "_coefficients.csv"),
            row.names = FALSE)

  # Predicted means: RhoticType x PhonCont x FormerRhotic
  grid_2b <- expand.grid(
    RhoticType_f   = levels(df_v$RhoticType_f),
    PhonCont_f     = levels(df_v$PhonCont_f),
    FormerRhotic_f = levels(df_v$FormerRhotic_f)
  )
  pred <- get_pred_means(m_2b, grid_2b, df_v, v)

  results_2b[[v]] <- list(model = m_2b, tab = tab,
                          loo = loo_2b, pred = pred)
  cat("\nModel 2b", v, "complete.\n")
}

# PP checks
pp_plots_2b <- list()
for (v in varieties) {
  pp_plots_2b[[v]] <- pp_check(results_2b[[v]]$model, nreps = 100) +
    labs(title = v) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
}
p_pp_2b <- wrap_plots(pp_plots_2b, nrow = 1) +
  plot_annotation(
    title    = "Posterior predictive checks — Model 2b DurStand (Bayesian)",
    subtitle = "Observed (dark) vs. 100 posterior predictive draws (light)",
    theme    = theme(plot.title = element_text(face = "bold"))
  )
ggsave("Fig_Model2b_Bayes_DurStand_PPcheck.pdf", p_pp_2b,
       width = 12, height = 5, dpi = 300)
ggsave("Fig_Model2b_Bayes_DurStand_PPcheck.jpg", p_pp_2b,
       width = 12, height = 5, dpi = 300)

# Predicted means plot: RhoticType x PhonCont, facet Variety x FormerRhotic
pred_all_2b <- bind_rows(lapply(results_2b, `[[`, "pred")) %>%
  mutate(
    Variety        = factor(Variety, levels = varieties),
    PhonCont_f     = factor(PhonCont_f,
                            levels = c("Initial", "Intervocalic", "Final")),
    FormerRhotic_f = factor(FormerRhotic_f,
                            levels = c("Vibrant", "Retroflex")),
    RhoticType_f   = factor(RhoticType_f, levels = c("TAP", "TRILL+"))
  )

p_pred_2b <- pred_all_2b %>%
  ggplot(aes(x = PhonCont_f, y = EMM,
             colour = RhoticType_f, group = RhoticType_f)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.15, linewidth = 0.8,
                position = position_dodge(width = 0.3)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.3)) +
  facet_grid(FormerRhotic_f ~ Variety) +
  scale_colour_manual(values = palette_rhotic, name = "Rhotic type") +
  labs(
    title    = "Model 2b (Bayesian) — Posterior predicted means of standardized duration (PreEU only)",
    subtitle = "stan_lmer | RhoticType + PhonCont + FormerRhotic | Prior: Normal(0,1) | Random intercept: Word\nRows: Vibrant / Retroflex | Columns: Variety | Dashed line = grand mean (z = 0)\n* Burgudži Retroflex x Final = 1 token",
    x        = "Phonological context",
    y        = "Posterior predicted mean (z-score)"
  ) +
  theme_romani()

print(p_pred_2b)
ggsave("Fig_Model2b_Bayes_DurStand_PredMeans.pdf", p_pred_2b,
       width = 12, height = 8, dpi = 300)
ggsave("Fig_Model2b_Bayes_DurStand_PredMeans.jpg", p_pred_2b,
       width = 12, height = 8, dpi = 300)

cat("\n=== SERIES 2 — BAYESIAN ANALYSIS COMPLETE ===\n")
cat("Output files:\n")
cat("  Tables : Tab_Model1/2a/2b_Bayes_DurStand_*_coefficients.csv\n")
cat("  Figures: Fig_Model1/2a/2b_Bayes_DurStand_*.pdf/.jpg\n")


##newgraphs##
# Ricolora: x = RhoticType, colore = FormerRhotic, facet = PhonCont x Variety

palette_former <- c("Vibrant"   = "#6ABF69",
                    "Retroflex" = "#B07CC6")

pred_all_2b <- bind_rows(lapply(results_2b, `[[`, "pred")) %>%
  mutate(
    Variety        = factor(Variety, levels = varieties),
    PhonCont_f     = factor(PhonCont_f,
                            levels = c("Initial", "Intervocalic", "Final")),
    FormerRhotic_f = factor(FormerRhotic_f,
                            levels = c("Vibrant", "Retroflex")),
    RhoticType_f   = factor(RhoticType_f, levels = c("TAP", "TRILL+"))
  )

p_pred_2b <- pred_all_2b %>%
  ggplot(aes(x = RhoticType_f, y = EMM,
             colour = FormerRhotic_f, group = FormerRhotic_f)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.15, linewidth = 0.8,
                position = position_dodge(width = 0.3)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.3)) +
  facet_grid(PhonCont_f ~ Variety) +
  scale_colour_manual(values = palette_former, name = "Former rhotic") +
  labs(
    title    = "Model 2b (Bayesian) — Posterior predicted means of standardized duration (PreEU only)",
    subtitle = "stan_lmer | RhoticType + PhonCont + FormerRhotic | Prior: Normal(0,1) | Random intercept: Word\nRows: Phonological context | Columns: Variety | Dashed line = grand mean (z = 0)\n* Burgudži Retroflex x Final = 1 token",
    x        = "Rhotic type",
    y        = "Posterior predicted mean (z-score)"
  ) +
  theme_romani()

ggsave("Fig_Model2b_Bayes_DurStand_PredMeans2.jpg", p_pred_2b,
       width = 12, height = 10, dpi = 300)
cat("Salvato!\n")