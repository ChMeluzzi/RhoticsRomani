# =============================================================================
# VARIAZIONE DELLE ROTICHE IN ROMANI — ANALISI DESCRITTIVA
# Variabili dipendenti: Rhotic (categoriale) | DurStand (continua)
# Predittori principali: PhonCont, Variety
# =============================================================================


# -----------------------------------------------------------------------------
# 0. PACCHETTI
# -----------------------------------------------------------------------------
# Decommentare solo alla prima esecuzione:
# install.packages(c("tidyverse", "rstatix", "ggplot2", "patchwork", "scales"))

library(readxl)
library(tidyverse)   # manipolazione dati + ggplot2
library(rstatix)     # cramer_v(), chisq_test()
library(patchwork)   # composizione multipanel dei grafici


# -----------------------------------------------------------------------------
# 1. CARICAMENTO E PULIZIA DEL DATASET
# -----------------------------------------------------------------------------

setwd("C:/Users/Fin/Documents/LAVORO/Pubblicazioni/AISV2026_Milano_Studi AISV15/Romani_AISV26/")

df_raw <- read_xlsx("romani_unified_new_Rstand290426_GM_230626.xlsx")
                   
                   
                   

cat("Token totali caricati:", nrow(df_raw), "\n")
cat("Colonne:", ncol(df_raw), "\n\n")


# --- 1a. DurStand: assicura tipo numerico (xlsx già in formato float) ---
df_raw <- df_raw %>%
  mutate(
    DurStand = as.numeric(DurStand),
    RhoticDur = as.numeric(RhoticDur)
  )

cat("DurStand — range:", round(range(df_raw$DurStand, na.rm = TRUE), 3), "\n\n")


# --- 1b. Ricodifica Rhotic: livelli ambigui → AMBIGUOUS ---
df_raw <- df_raw %>%
  mutate(
    Rhotic_clean = case_when(
      Rhotic %in% c("TAP + FRIC", "TRILL + FRIC") ~ "AMBIGUOUS",
      TRUE ~ Rhotic
    )
  )

cat("Distribuzione Rhotic_clean (inclusi NA):\n")
print(table(df_raw$Rhotic_clean, useNA = "ifany"))
cat("\n")


# --- 1c. Dataset pulito per analisi Rhotic (escludi NA su variabili chiave) ---
df <- df_raw %>%
  filter(
    !is.na(Rhotic_clean),
    !is.na(PhonCont),
    !is.na(Variety)
  ) %>%
  mutate(
    # Fattori con livelli di riferimento espliciti
    Rhotic_f   = factor(Rhotic_clean,
                        levels = c("TAP", "APPROX", "TRILL", "FRIC", "AMBIGUOUS")),
    PhonCont_f = factor(PhonCont,
                        levels = c("Intervocalic", "Initial", "Final")),
    Variety_f  = factor(Variety,
                        levels = c("Burgudzi", "Kalajdzi", "Kaldaras"))
  )

cat("Token dopo rimozione NA su Rhotic/PhonCont/Variety:", nrow(df), "\n\n")
cat("=== Distribuzione finale ===\n")
cat("Rhotic_f:\n");    print(table(df$Rhotic_f))
cat("\nPhonCont_f:\n"); print(table(df$PhonCont_f))
cat("\nVariety_f:\n");  print(table(df$Variety_f))


# =============================================================================
# 2. ANALISI DESCRITTIVA — RHOTIC (CATEGORIALE)
# =============================================================================

# --- 2a. Funzione tabella di contingenza + χ² + V di Cramér ---

tabella_contingenza <- function(var_riga, var_col, nome_riga, nome_col, data) {
  cat("\n==========================================\n")
  cat(nome_riga, "×", nome_col, "\n")
  cat("------------------------------------------\n")

  tab <- table(data[[var_riga]], data[[var_col]])

  # Tabella con totali
  cat("\nConteggi (con totali di riga e colonna):\n")
  print(addmargins(tab))

  # Proporzioni per riga (%)
  cat("\nProporzioni per riga (%):\n")
  print(round(prop.table(tab, margin = 1) * 100, 1))

  # Chi-square
  chi <- chisq.test(tab)
  cat("\nChi-square:", round(chi$statistic, 3),
      "  df:", chi$parameter,
      "  p =", format.pval(chi$p.value, digits = 4), "\n")

  # Cramér's V (via rstatix)
  v <- cramer_v(tab)
  cat("Cramér's V:", round(v, 3),
      " (0.1 = piccolo, 0.3 = medio, 0.5 = grande)\n")
  cat("==========================================\n\n")

  invisible(list(tab = tab, chi = chi, V = v))
}


# --- Rhotic × PhonCont ---
res_rph <- tabella_contingenza("Rhotic_f", "PhonCont_f",
                               "Rhotic", "PhonCont", df)

# --- Rhotic × Variety ---
res_rv  <- tabella_contingenza("Rhotic_f", "Variety_f",
                               "Rhotic", "Variety", df)


# --- 2b. Grafico: proporzioni di Rhotic per PhonCont × Variety (stacked bar) ---

# Palette colori per i tipi di rotica
palette_rhotic <- c(
  "TAP"       = "#4E9AF1",
  "TRILL"     = "#E86F51",
  "APPROX"    = "#6ABF69",
  "FRIC"      = "#F4C430",
  "AMBIGUOUS" = "#B0B0B0"
)

p_rhotic <- df %>%
  count(Variety_f, PhonCont_f, Rhotic_f) %>%
  group_by(Variety_f, PhonCont_f) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = PhonCont_f, y = prop, fill = Rhotic_f)) +
  geom_col(position = "stack", width = 0.75, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(prop >= 0.05,
                               paste0(round(prop * 100), "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 3, colour = "white", fontface = "bold") +
  facet_wrap(~ Variety_f, nrow = 1) +
  scale_fill_manual(values = palette_rhotic, name = "Tipo di rotica") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(
    title    = "Distribuzione dei tipi di rotica per contesto fonologico e varietà",
    subtitle = "Proporzioni calcolate per cella Varietà × Contesto",
    x        = "Contesto fonologico",
    y        = "Proporzione (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.major.x = element_blank(),
    strip.text       = element_text(face = "bold", size = 12),
    plot.title       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 20, hjust = 1)
  )

print(p_rhotic)
ggsave("Fig1_Rhotic_PhonCont_Variety.pdf", p_rhotic,
       width = 10, height = 6, dpi = 300)
cat("Grafico salvato: Fig1_Rhotic_PhonCont_Variety.pdf\n")


# =============================================================================
# 3. ANALISI DESCRITTIVA — DURSTAND (CONTINUA)
# =============================================================================

# Dataset per DurStand: escludi NA su DurStand, PhonCont, Variety
df_dur <- df_raw %>%
  filter(
    !is.na(DurStand),
    !is.na(PhonCont),
    !is.na(Variety)
  ) %>%
  mutate(
    PhonCont_f = factor(PhonCont,
                        levels = c("Intervocalic", "Initial", "Final")),
    Variety_f  = factor(Variety,
                        levels = c("Burgudzi", "Kalajdzi", "Kaldaras"))
  )

cat("\nToken per analisi DurStand:", nrow(df_dur), "\n\n")

# --- 3a. Statistiche descrittive per Variety × PhonCont ---
cat("=== Statistiche descrittive DurStand (Variety × PhonCont) ===\n")
df_dur %>%
  group_by(Variety_f, PhonCont_f) %>%
  summarise(
    N      = n(),
    Media  = round(mean(DurStand), 3),
    Mediana= round(median(DurStand), 3),
    SD     = round(sd(DurStand), 3),
    Min    = round(min(DurStand), 3),
    Max    = round(max(DurStand), 3),
    .groups = "drop"
  ) %>%
  print(n = Inf)


# --- 3b. Grafici a violino: 3 panel (uno per varietà), x = PhonCont ---

palette_phonCont <- c(
  "Intervocalic" = "#4E9AF1",
  "Initial"      = "#E86F51",
  "Final"        = "#6ABF69"
)

# Funzione per costruire il violino di una singola varietà
violino_variety <- function(variety_name, colore_fill) {

  df_sub <- df_dur %>% filter(Variety_f == variety_name)

  # Statistiche per le etichette mediana
  stat_med <- df_sub %>%
    group_by(PhonCont_f) %>%
    summarise(med = median(DurStand), .groups = "drop")

  ggplot(df_sub, aes(x = PhonCont_f, y = DurStand, fill = PhonCont_f)) +
    geom_violin(trim = FALSE, alpha = 0.6, colour = "grey40", linewidth = 0.4) +
    geom_boxplot(width = 0.12, outlier.shape = NA,
                 colour = "grey20", fill = "white", linewidth = 0.5) +
    geom_jitter(width = 0.08, size = 1, alpha = 0.35, colour = "grey30") +
    geom_text(data = stat_med,
              aes(x = PhonCont_f, y = med, label = round(med, 2)),
              vjust = -0.6, size = 3, fontface = "bold", colour = "grey10") +
    scale_fill_manual(values = palette_phonCont, guide = "none") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    labs(
      title = variety_name,
      x     = "Contesto fonologico",
      y     = "Durata standardizzata (z-score)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
      panel.grid.major.x = element_blank(),
      axis.text.x      = element_text(size = 10)
    )
}

p_burg  <- violino_variety("Burgudzi",  "#4E9AF1")
p_kalaj <- violino_variety("Kalajdzi",  "#E86F51")
p_kald  <- violino_variety("Kaldaras",  "#6ABF69")

# Composizione multipanel con patchwork
p_dur <- p_burg | p_kalaj | p_kald

p_dur_finale <- p_dur +
  plot_annotation(
    title    = "Distribuzione della durata standardizzata per varietà e contesto fonologico",
    subtitle = "Violini + boxplot + singole osservazioni | mediana annotata",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, colour = "grey40")
    )
  )

print(p_dur_finale)
ggsave("Fig2_DurStand_Violini.pdf", p_dur_finale,
       width = 14, height = 6, dpi = 300)
cat("Grafico salvato: Fig2_DurStand_Violini.pdf\n")

cat("\n=== ANALISI DESCRITTIVA COMPLETATA ===\n")
