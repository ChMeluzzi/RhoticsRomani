# =============================================================================
# RHOTIC VARIATION IN ROMANI — EXPLORATORY ANALYSIS OF NEW VARIABLES
# Variables: Origin, OriginAggregate, OriginUltraAggregate, FormerRhotic
# Dataset: romani_unified_new_Rstand290426_GM_230626.xlsx (updated 29/04/2026)
# =============================================================================


# -----------------------------------------------------------------------------
# 0. PACKAGES
# -----------------------------------------------------------------------------
# install.packages(c("tidyverse", "rstatix", "patchwork", "scales"))

library(readxl)
library(tidyverse)
library(rstatix)    # cramer_v()
library(patchwork)
library(scales)


# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------

setwd("C:/Users/Fin/Documents/LAVORO/Pubblicazioni/AISV2026_Milano_Studi AISV15/Romani_AISV26/")

df <- read_xlsx("romani_unified_new_Rstand290426_GM_230626.xlsx")
               
                %>%
  mutate(
    DurStand            = as.numeric(gsub(",", ".", DurStand)),
    Variety_f           = factor(Variety,
                                 levels = c("Burgudzi", "Kalajdzi", "Kaldaras")),
    FormerRhotic_f      = factor(FormerRhotic,
                                 levels = c("Vibrant", "Retroflex")),
    OriginAggregate_f   = factor(OriginAggregate,
                                 levels = c("Indian+", "Greek", "Turkish+",
                                            "Slavic", "Balkan", "Rumenian",
                                            "Uncertain_PreEU", "Uncertain_EU")),
    OriginUltraAggregate_f = factor(OriginUltraAggregate,
                                    levels = c("PreEU", "EU"))
  )

cat("=== DATASET ===\n")
cat("Total tokens:", nrow(df), "\n")
cat("Columns     :", ncol(df), "\n\n")


# -----------------------------------------------------------------------------
# 2. DISTRIBUTION OF ORIGIN VARIABLES
# -----------------------------------------------------------------------------

cat("=== Origin (full, 35 levels) ===\n")
origin_tab <- sort(table(df$Origin), decreasing = TRUE)
print(origin_tab)
cat("\n")

cat("=== OriginAggregate (9 levels, incl. NA) ===\n")
print(table(df$OriginAggregate_f, useNA = "ifany"))
cat("\n")

cat("=== OriginUltraAggregate (2 levels, incl. NA) ===\n")
print(table(df$OriginUltraAggregate_f, useNA = "ifany"))
cat("\n")

cat("=== FormerRhotic ===\n")
print(table(df$FormerRhotic_f, useNA = "ifany"))
cat("\n")


# -----------------------------------------------------------------------------
# 3. CROSSTABS: FormerRhotic x OriginAggregate / OriginUltraAggregate
# -----------------------------------------------------------------------------

cat("=== FormerRhotic x OriginAggregate (counts) ===\n")
tab_fr_oa <- table(df$FormerRhotic_f, df$OriginAggregate_f)
print(addmargins(tab_fr_oa))
cat("\nProportions by row (%):\n")
print(round(prop.table(tab_fr_oa, margin = 1) * 100, 1))

chi_fr_oa <- chisq.test(tab_fr_oa)
cat("\nChi-square:", round(chi_fr_oa$statistic, 3),
    "  df:", chi_fr_oa$parameter,
    "  p =", format.pval(chi_fr_oa$p.value, digits = 4), "\n")
cat("Cramér's V:", round(cramer_v(tab_fr_oa), 3), "\n\n")

cat("=== FormerRhotic x OriginUltraAggregate (counts) ===\n")
tab_fr_oua <- table(df$FormerRhotic_f, df$OriginUltraAggregate_f)
print(addmargins(tab_fr_oua))
cat("\nProportions by row (%):\n")
print(round(prop.table(tab_fr_oua, margin = 1) * 100, 1))

chi_fr_oua <- chisq.test(tab_fr_oua)
cat("\nChi-square:", round(chi_fr_oua$statistic, 3),
    "  df:", chi_fr_oua$parameter,
    "  p =", format.pval(chi_fr_oua$p.value, digits = 4), "\n")
cat("Cramér's V:", round(cramer_v(tab_fr_oua), 3), "\n\n")


# -----------------------------------------------------------------------------
# 4. CROSSTABS: OriginAggregate / OriginUltraAggregate x Variety
# -----------------------------------------------------------------------------

cat("=== OriginAggregate x Variety (counts) ===\n")
tab_oa_var <- table(df$OriginAggregate_f, df$Variety_f)
print(addmargins(tab_oa_var))
cat("\nProportions by row (%):\n")
print(round(prop.table(tab_oa_var, margin = 1) * 100, 1))
cat("\n")

cat("=== OriginUltraAggregate x Variety (counts) ===\n")
tab_oua_var <- table(df$OriginUltraAggregate_f, df$Variety_f)
print(addmargins(tab_oua_var))
cat("\nProportions by row (%):\n")
print(round(prop.table(tab_oua_var, margin = 1) * 100, 1))
cat("\n")


# -----------------------------------------------------------------------------
# 5. PLOTS
# -----------------------------------------------------------------------------

palette_origin <- c(
  "Indian+"        = "#4E9AF1",
  "Greek"          = "#E86F51",
  "Turkish+"       = "#F4C430",
  "Slavic"         = "#6ABF69",
  "Balkan"         = "#B07CC6",
  "Rumenian"       = "#F08080",
  "Uncertain_PreEU"= "#A0A0A0",
  "Uncertain_EU"   = "#D0D0D0"
)

palette_ultra <- c(
  "PreEU" = "#4E9AF1",
  "EU"    = "#E86F51"
)

# --- 5a. OriginAggregate distribution by Variety ---

p_oa_var <- df %>%
  filter(!is.na(OriginAggregate_f)) %>%
  count(Variety_f, OriginAggregate_f) %>%
  group_by(Variety_f) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = Variety_f, y = prop, fill = OriginAggregate_f)) +
  geom_col(position = "stack", width = 0.7,
           colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(prop >= 0.04,
                               paste0(round(prop * 100), "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 3, colour = "white", fontface = "bold") +
  scale_fill_manual(values = palette_origin, name = "Origin") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "OriginAggregate distribution by variety",
    x     = "Variety",
    y     = "Proportion (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "right",
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold")
  )

# --- 5b. OriginUltraAggregate distribution by Variety ---

p_oua_var <- df %>%
  filter(!is.na(OriginUltraAggregate_f)) %>%
  count(Variety_f, OriginUltraAggregate_f) %>%
  group_by(Variety_f) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = Variety_f, y = prop, fill = OriginUltraAggregate_f)) +
  geom_col(position = "stack", width = 0.7,
           colour = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(prop * 100), "%")),
            position = position_stack(vjust = 0.5),
            size = 3.5, colour = "white", fontface = "bold") +
  scale_fill_manual(values = palette_ultra, name = "Origin") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "OriginUltraAggregate distribution by variety",
    x     = "Variety",
    y     = "Proportion (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "right",
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold")
  )

# --- 5c. FormerRhotic x OriginAggregate ---

p_fr_oa <- df %>%
  filter(!is.na(OriginAggregate_f), !is.na(FormerRhotic_f)) %>%
  count(OriginAggregate_f, FormerRhotic_f) %>%
  group_by(OriginAggregate_f) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = OriginAggregate_f, y = prop, fill = FormerRhotic_f)) +
  geom_col(position = "stack", width = 0.7,
           colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(prop >= 0.05,
                               paste0(round(prop * 100), "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 3, colour = "white", fontface = "bold") +
  scale_fill_manual(values = c("Vibrant" = "#4E9AF1",
                               "Retroflex" = "#E86F51"),
                    name = "Former rhotic") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "FormerRhotic distribution by OriginAggregate",
    x     = "Origin (aggregated)",
    y     = "Proportion (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "right",
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 25, hjust = 1)
  )

# --- 5d. FormerRhotic x OriginUltraAggregate ---

p_fr_oua <- df %>%
  filter(!is.na(OriginUltraAggregate_f), !is.na(FormerRhotic_f)) %>%
  count(OriginUltraAggregate_f, FormerRhotic_f) %>%
  group_by(OriginUltraAggregate_f) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = OriginUltraAggregate_f, y = prop, fill = FormerRhotic_f)) +
  geom_col(position = "stack", width = 0.5,
           colour = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(prop * 100), "%")),
            position = position_stack(vjust = 0.5),
            size = 3.5, colour = "white", fontface = "bold") +
  scale_fill_manual(values = c("Vibrant" = "#4E9AF1",
                               "Retroflex" = "#E86F51"),
                    name = "Former rhotic") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "FormerRhotic distribution by OriginUltraAggregate",
    x     = "Origin (ultra-aggregated)",
    y     = "Proportion (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "right",
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold")
  )

# --- Compose and save ---

p_origin <- (p_oa_var | p_oua_var) /
            (p_fr_oa  | p_fr_oua) +
  plot_annotation(
    title = "Distribution of etymological origin and former rhotic variables",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

print(p_origin)
ggsave("Fig_Exploratory_Origin_FormerRhotic.pdf", p_origin,
       width = 14, height = 10, dpi = 300)
ggsave("Fig_Exploratory_Origin_FormerRhotic.jpg", p_origin,
       width = 14, height = 10, dpi = 300)
cat("Plot saved: Fig_Exploratory_Origin_FormerRhotic\n")

cat("\n=== EXPLORATORY ANALYSIS COMPLETE ===\n")
