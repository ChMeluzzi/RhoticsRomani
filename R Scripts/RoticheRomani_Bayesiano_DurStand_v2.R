# =============================================================================
# RHOTIC VARIATION IN ROMANI — BAYESIAN INFERENTIAL ANALYSIS
# Series 2: DurStand (standardized duration, continuous)
# =============================================================================
# Dependent variable: DurStand (z-score standardized duration)
# Random intercept : (1 | Word)
# Framework        : Bayesian, rstanarm::stan_lmer()
# Prior            : Normal(0, 1) on z-score scale
#                    (weakly informative — effects > 2 SD implausible a priori)
#
# MODEL 1 (stepping stone):
#   DurStand ~ PhonCont + Variety + (1|Word)
#   Full dataset | All varieties pooled
#   Reference: Intervocalic | Burgudzi
#
# MODEL 2a:
#   DurStand ~ PhonCont * OriginUltraAggregate + (1|Word)
#   Full dataset | Fit separately per variety
#   Reference: Intervocalic | PreEU
#
# MODEL 2b:
#   DurStand ~ PhonCont * FormerRhotic + (1|Word)
#   PreEU tokens only | Fit separately per variety
#   Reference: Intervocalic | Vibrant
#   NOTE: FormerRhotic has no variation in EU stratum (all Vibrant)
#   NOTE: Burgudzi Final x Retroflex = 1 token — interpret with caution
#
# Model comparison: LOO-CV (Vehtari et al., 2017)
#
# RATIONALE FOR BAYESIAN APPROACH:
# Methodological consistency with Series 1 (RhoticType); regularization
# for cells with sparse data; LOO-CV enables direct model comparison
# across all models within a unified framework.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. PACKAGES
# -----------------------------------------------------------------------------
# Uncomment at first run:
# install.packages(c("rstanarm", "bayesplot", "tidybayes",
#                    "tidyverse", "patchwork", "scales"))

library(readxl)
library(tidyverse)
library(rstanarm)
library(bayesplot)
library(tidybayes)
library(patchwork)
library(scales)


# -----------------------------------------------------------------------------
# 1. DATASET
# -----------------------------------------------------------------------------

setwd("C:/Users/Fin/Documents/LAVORO/Pubblicazioni/AISV2026_Milano_Studi AISV15/Romani_AISV26/")

df_raw <- read_xlsx("romani_unified_new_Rstand290426_GM_230626.xlsx")
                   
                    %>%
  mutate(DurStand = as.numeric(DurStand))

# Base dataset: complete cases
df_base <- df_raw %>%
  filter(
    !is.na(DurStand),
    !is.na(PhonCont),
    !is.na(OriginUltraAggregate),
    !is.na(FormerRhotic),
    !is.na(Word)
  ) %>%
  mutate(
    PhonCont_f             = factor(PhonCont,
                                    levels = c("Intervocalic", "Initial", "Final")),
    Variety_f              = factor(Variety,
                                    levels = c("Burgudzi", "Kalajdzi", "Kaldaras")),
    OriginUltraAggregate_f = factor(OriginUltraAggregate,
                                    levels = c("PreEU", "EU")),
    FormerRhotic_f         = factor(FormerRhotic,
                                    levels = c("Vibrant", "Retroflex")),
    Word                   = factor(Word)
  )

# PreEU subset for Model 2b
df_preeu <- df_base %>% filter(OriginUltraAggregate == "PreEU")

varieties <- c("Burgudzi", "Kalajdzi", "Kaldaras")

cat("=== DATASET ===\n")
cat("Full (excl. NA):", nrow(df_base), "tokens\n")
cat("PreEU only     :", nrow(df_preeu), "tokens\n\n")
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
             starts_with("OriginUltra"))

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
palette_variety  <- c("Burgudzi" = "#4E9AF1",
                      "Kalajdzi" = "#E86F51",
                      "Kaldaras" = "#6ABF69")


# =============================================================================
# 3. MODEL 1 — Stepping stone
#    DurStand ~ PhonCont + Variety + (1|Word)
#    Full dataset, all varieties pooled
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("MODEL 1 — Stepping stone (all varieties)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

m1 <- stan_lmer(
  DurStand ~ PhonCont_f + Variety_f + (1 | Word),
  data            = df_base,
  prior           = normal(0, 1),
  prior_intercept = normal(0, 1),
  seed            = 42,
  iter            = 4000,
  chains          = 4,
  cores           = 4,
  refresh         = 100
)

cat("\n--- Posterior summary: Model 1 ---\n")
print(summary(m1,
              pars = c("(Intercept)",
                       "PhonCont_fInitial", "PhonCont_fFinal",
                       "Variety_fKalajdzi", "Variety_fKaldaras"),
              probs = c(0.025, 0.975),
              digits = 3))

# NOTE: Word random effect variance expected near zero (singular fit
# observed in frequentist version). Retained for theoretical consistency.
cat("\n--- Random effects ---\n")
print(VarCorr(m1))

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

# Predicted means
grid_m1 <- expand.grid(
  PhonCont_f = levels(df_base$PhonCont_f),
  Variety_f  = levels(df_base$Variety_f)
)
pred_m1 <- get_pred_means(m1, grid_m1, df_base, "All") %>%
  mutate(
    Variety_f  = factor(Variety_f, levels = varieties),
    PhonCont_f = factor(PhonCont_f,
                        levels = c("Intervocalic", "Initial", "Final"))
  )

p_pred_m1 <- pred_m1 %>%
  ggplot(aes(x = PhonCont_f, y = EMM,
             colour = PhonCont_f)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.15, linewidth = 0.8) +
  geom_point(size = 4) +
  facet_wrap(~ Variety_f, nrow = 1) +
  scale_colour_manual(values = palette_phonCont, guide = "none") +
  labs(
    title    = "Model 1 (Bayesian) — Posterior predicted means of standardized duration",
    subtitle = "stan_lmer | PhonCont + Variety | Prior: Normal(0,1) | Random intercept: Word\nDashed line = grand mean (z = 0)",
    x        = "Phonological context",
    y        = "Posterior predicted mean (z-score)"
  ) +
  theme_romani()

print(p_pred_m1)
ggsave("Fig_Model1_Bayes_DurStand_PredMeans.pdf", p_pred_m1,
       width = 10, height = 6, dpi = 300)
ggsave("Fig_Model1_Bayes_DurStand_PredMeans.jpg", p_pred_m1,
       width = 10, height = 6, dpi = 300)
cat("Model 1 complete.\n")


# =============================================================================
# 4. MODEL 2a — OriginUltraAggregate
#    DurStand ~ PhonCont * OriginUltraAggregate + (1|Word)
#    Full dataset | Fit separately per variety
# =============================================================================

results_2a <- list()

for (v in varieties) {

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("MODEL 2a —", v, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  df_v <- df_base %>% filter(Variety == v)

  cat("N =", nrow(df_v),
      "| EU:", sum(df_v$OriginUltraAggregate == "EU"),
      "| PreEU:", sum(df_v$OriginUltraAggregate == "PreEU"), "\n\n")

  m_2a <- stan_lmer(
    DurStand ~ PhonCont_f * OriginUltraAggregate_f + (1 | Word),
    data            = df_v,
    prior           = normal(0, 1),
    prior_intercept = normal(0, 1),
    seed            = 42,
    iter            = 4000,
    chains          = 4,
    cores           = 4,
    refresh         = 100
  )

  cat("\n--- Posterior summary:", v, "---\n")
  print(summary(m_2a, probs = c(0.025, 0.975), digits = 3))

  # LOO-CV
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

  # Predicted means
  grid_2a <- expand.grid(
    PhonCont_f             = levels(df_v$PhonCont_f),
    OriginUltraAggregate_f = levels(df_v$OriginUltraAggregate_f)
  )
  pred <- get_pred_means(m_2a, grid_2a, df_v, v)

  results_2a[[v]] <- list(model = m_2a, tab = tab,
                          loo = loo_2a, pred = pred)
  cat("\nModel 2a", v, "complete.\n")
}

# LOO comparison: Model 1 vs Model 2a per variety
cat("\n=== LOO-CV COMPARISON: Model 1 vs. Model 2a (per variety) ===\n")
for (v in varieties) {
  df_v <- df_base %>% filter(Variety == v)
  m1_v <- stan_lmer(
    DurStand ~ PhonCont_f + (1 | Word),
    data            = df_v,
    prior           = normal(0, 1),
    prior_intercept = normal(0, 1),
    seed            = 42,
    iter            = 4000,
    chains          = 4,
    cores           = 4,
    refresh         = 0
  )
  loo_m1_v  <- loo(m1_v)
  loo_m2a_v <- results_2a[[v]]$loo
  cat("\n---", v, "---\n")
  print(loo_compare(loo_m1_v, loo_m2a_v))
}

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

# Predicted means plot
pred_all_2a <- bind_rows(lapply(results_2a, `[[`, "pred")) %>%
  mutate(
    Variety                = factor(Variety, levels = varieties),
    PhonCont_f             = factor(PhonCont_f,
                                    levels = c("Intervocalic", "Initial", "Final")),
    OriginUltraAggregate_f = factor(OriginUltraAggregate_f,
                                    levels = c("PreEU", "EU"))
  )

p_pred_2a <- pred_all_2a %>%
  ggplot(aes(x = PhonCont_f, y = EMM, colour = PhonCont_f)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.15, linewidth = 0.8) +
  geom_point(size = 4) +
  facet_grid(OriginUltraAggregate_f ~ Variety) +
  scale_colour_manual(values = palette_phonCont, guide = "none") +
  labs(
    title    = "Model 2a (Bayesian) — Posterior predicted means of standardized duration",
    subtitle = "stan_lmer | PhonCont × OriginUltraAggregate | Prior: Normal(0,1) | Random intercept: Word\nRows: PreEU / EU | Columns: Variety | Dashed line = grand mean (z = 0)",
    x        = "Phonological context",
    y        = "Posterior predicted mean (z-score)"
  ) +
  theme_romani()

print(p_pred_2a)
ggsave("Fig_Model2a_Bayes_DurStand_PredMeans.pdf", p_pred_2a,
       width = 12, height = 8, dpi = 300)
ggsave("Fig_Model2a_Bayes_DurStand_PredMeans.jpg", p_pred_2a,
       width = 12, height = 8, dpi = 300)

# Coefficient plot
tab_all_2a <- bind_rows(lapply(results_2a, `[[`, "tab")) %>%
  filter(!grepl("Intercept", Term)) %>%
  mutate(
    Term = case_when(
      Term == "PhonCont_fInitial"                          ~ "Initial\n(vs. Intervocalic)",
      Term == "PhonCont_fFinal"                            ~ "Final\n(vs. Intervocalic)",
      Term == "OriginUltraAggregate_fEU"                   ~ "EU\n(vs. PreEU)",
      Term == "PhonCont_fInitial:OriginUltraAggregate_fEU" ~ "Initial ×\nEU",
      Term == "PhonCont_fFinal:OriginUltraAggregate_fEU"   ~ "Final ×\nEU",
      TRUE ~ Term
    ),
    Term    = factor(Term, levels = rev(c(
      "Initial\n(vs. Intervocalic)",
      "Final\n(vs. Intervocalic)",
      "EU\n(vs. PreEU)",
      "Initial ×\nEU",
      "Final ×\nEU"
    ))),
    Variety = factor(Variety, levels = varieties)
  )

p_coef_2a <- tab_all_2a %>%
  ggplot(aes(x = Median, y = Term, colour = Variety)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high),
                width = 0.25, linewidth = 0.7,
                position = position_dodge(width = 0.6)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = palette_variety, name = "Variety") +
  labs(
    title    = "Model 2a (Bayesian) — Posterior estimates with 95% credible intervals",
    subtitle = "z-score scale | Reference: Intervocalic | PreEU",
    x        = "Posterior median (z-score)",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(colour = "grey40", size = 9),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom")

print(p_coef_2a)
ggsave("Fig_Model2a_Bayes_DurStand_Coef.pdf", p_coef_2a,
       width = 10, height = 6, dpi = 300)
ggsave("Fig_Model2a_Bayes_DurStand_Coef.jpg", p_coef_2a,
       width = 10, height = 6, dpi = 300)
cat("\nModel 2a complete.\n")


# =============================================================================
# 5. MODEL 2b — FormerRhotic
#    DurStand ~ PhonCont * FormerRhotic + (1|Word)
#    PreEU tokens only | Fit separately per variety
# =============================================================================

results_2b <- list()

for (v in varieties) {

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("MODEL 2b —", v, "(PreEU only)\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  df_v <- df_preeu %>% filter(Variety == v)

  cat("N =", nrow(df_v),
      "| Vibrant:", sum(df_v$FormerRhotic == "Vibrant"),
      "| Retroflex:", sum(df_v$FormerRhotic == "Retroflex"), "\n\n")

  if (v == "Burgudzi") {
    cat("NOTE: Final x Retroflex = 1 token.",
        "Interpret Final x Retroflex posterior with caution.\n\n")
  }

  m_2b <- stan_lmer(
    DurStand ~ PhonCont_f * FormerRhotic_f + (1 | Word),
    data            = df_v,
    prior           = normal(0, 1),
    prior_intercept = normal(0, 1),
    seed            = 42,
    iter            = 4000,
    chains          = 4,
    cores           = 4,
    refresh         = 100
  )

  cat("\n--- Posterior summary:", v, "---\n")
  print(summary(m_2b, probs = c(0.025, 0.975), digits = 3))

  # LOO-CV
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

  # Predicted means
  grid_2b <- expand.grid(
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

# Predicted means plot
pred_all_2b <- bind_rows(lapply(results_2b, `[[`, "pred")) %>%
  mutate(
    Variety        = factor(Variety, levels = varieties),
    PhonCont_f     = factor(PhonCont_f,
                            levels = c("Intervocalic", "Initial", "Final")),
    FormerRhotic_f = factor(FormerRhotic_f,
                            levels = c("Vibrant", "Retroflex"))
  )

p_pred_2b <- pred_all_2b %>%
  ggplot(aes(x = PhonCont_f, y = EMM, colour = PhonCont_f)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.15, linewidth = 0.8) +
  geom_point(size = 4) +
  facet_grid(FormerRhotic_f ~ Variety) +
  scale_colour_manual(values = palette_phonCont, guide = "none") +
  labs(
    title    = "Model 2b (Bayesian) — Posterior predicted means of standardized duration (PreEU only)",
    subtitle = "stan_lmer | PhonCont × FormerRhotic | Prior: Normal(0,1) | Random intercept: Word\nRows: Vibrant / Retroflex | Columns: Variety | Dashed line = grand mean (z = 0)",
    x        = "Phonological context",
    y        = "Posterior predicted mean (z-score)"
  ) +
  theme_romani()

print(p_pred_2b)
ggsave("Fig_Model2b_Bayes_DurStand_PredMeans.pdf", p_pred_2b,
       width = 12, height = 8, dpi = 300)
ggsave("Fig_Model2b_Bayes_DurStand_PredMeans.jpg", p_pred_2b,
       width = 12, height = 8, dpi = 300)

# Coefficient plot
tab_all_2b <- bind_rows(lapply(results_2b, `[[`, "tab")) %>%
  filter(!grepl("Intercept", Term)) %>%
  mutate(
    Term = case_when(
      Term == "PhonCont_fInitial"                         ~ "Initial\n(vs. Intervocalic)",
      Term == "PhonCont_fFinal"                           ~ "Final\n(vs. Intervocalic)",
      Term == "FormerRhotic_fRetroflex"                   ~ "Retroflex\n(vs. Vibrant)",
      Term == "PhonCont_fInitial:FormerRhotic_fRetroflex" ~ "Initial ×\nRetroflex",
      Term == "PhonCont_fFinal:FormerRhotic_fRetroflex"   ~ "Final ×\nRetroflex",
      TRUE ~ Term
    ),
    Term    = factor(Term, levels = rev(c(
      "Initial\n(vs. Intervocalic)",
      "Final\n(vs. Intervocalic)",
      "Retroflex\n(vs. Vibrant)",
      "Initial ×\nRetroflex",
      "Final ×\nRetroflex"
    ))),
    Variety = factor(Variety, levels = varieties)
  )

p_coef_2b <- tab_all_2b %>%
  ggplot(aes(x = Median, y = Term, colour = Variety)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high),
                width = 0.25, linewidth = 0.7,
                position = position_dodge(width = 0.6)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = palette_variety, name = "Variety") +
  labs(
    title    = "Model 2b (Bayesian) — Posterior estimates with 95% credible intervals",
    subtitle = "z-score scale | Reference: Intervocalic | Vibrant | PreEU only",
    x        = "Posterior median (z-score)",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(colour = "grey40", size = 9),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom")

print(p_coef_2b)
ggsave("Fig_Model2b_Bayes_DurStand_Coef.pdf", p_coef_2b,
       width = 10, height = 6, dpi = 300)
ggsave("Fig_Model2b_Bayes_DurStand_Coef.jpg", p_coef_2b,
       width = 10, height = 6, dpi = 300)

cat("\n=== SERIES 2 — BAYESIAN ANALYSIS COMPLETE ===\n")
cat("Output files:\n")
cat("  Tables : Tab_Model1/2a/2b_Bayes_DurStand_*_coefficients.csv\n")
cat("  Figures: Fig_Model1/2a/2b_Bayes_DurStand_*.pdf/.jpg\n")
