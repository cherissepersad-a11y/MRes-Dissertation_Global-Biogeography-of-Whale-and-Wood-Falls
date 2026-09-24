# ============================================================
# Supplementary Figure S2
# Publication history of analytical organic-fall evidence base
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# 1. Analytical whale/wood events
# ------------------------------------------------------------

event_key <- events %>%
  filter(
    Analysis_Include == 1,
    Fall_Type %in% c("Whale", "Wood")
  ) %>%
  select(Event_ID, Fall_Type)

stopifnot(!anyDuplicated(event_key$Event_ID))

# ------------------------------------------------------------
# 2. Link analytical occurrences to their sources
# ------------------------------------------------------------

analytical_occ_sources <- occurrences %>%
  inner_join(event_key, by = "Event_ID") %>%
  filter(
    !is.na(Source_ID),
    Source_ID != ""
  )

stopifnot(nrow(analytical_occ_sources) == 6112)

# ------------------------------------------------------------
# 3. Classify each source as Whale, Wood or Both
# ------------------------------------------------------------

source_fall <- analytical_occ_sources %>%
  distinct(Source_ID, Fall_Type) %>%
  group_by(Source_ID) %>%
  summarise(
    Evidence_Type = case_when(
      all(c("Whale", "Wood") %in% Fall_Type) ~ "Both",
      "Whale" %in% Fall_Type ~ "Whale",
      "Wood"  %in% Fall_Type ~ "Wood",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  )

stopifnot(nrow(source_fall) == 226)

# ------------------------------------------------------------
# 4. Join source metadata and recover publication year
# ------------------------------------------------------------

source_history <- source_fall %>%
  left_join(sources, by = "Source_ID") %>%
  mutate(
    Year_explicit =
      suppressWarnings(as.integer(Year)),
    
    Year_from_citation =
      suppressWarnings(
        as.integer(
          str_extract(
            Source_Citation,
            "(?<!\\d)(18|19|20)\\d{2}(?!\\d)"
          )
        )
      ),
    
    Publication_Year =
      coalesce(Year_explicit, Year_from_citation)
  )

# ------------------------------------------------------------
# 5. Restrict publication-history analysis to dated articles
# ------------------------------------------------------------

publication_history <- source_history %>%
  filter(
    Source_Type == "Article",
    !is.na(Publication_Year)
  ) %>%
  distinct(Source_ID, .keep_all = TRUE) %>%
  mutate(
    Evidence_Type = factor(
      Evidence_Type,
      levels = c("Whale", "Both", "Wood")
    )
  )

cat("\nDated peer-reviewed articles used in S2:",
    nrow(publication_history), "\n")

cat("\nBy evidence type:\n")
print(table(publication_history$Evidence_Type))

cat("\nYear range:\n")
print(range(publication_history$Publication_Year))

# ------------------------------------------------------------
# 6. Annual publication counts
# ------------------------------------------------------------

year_range <- seq(
  min(publication_history$Publication_Year),
  max(publication_history$Publication_Year)
)

annual <- publication_history %>%
  count(Publication_Year, Evidence_Type, name = "n") %>%
  complete(
    Publication_Year = year_range,
    Evidence_Type,
    fill = list(n = 0)
  )

# ------------------------------------------------------------
# 7. Cumulative publication count
#
# Each publication counted ONCE, including papers that support
# both whale- and wood-fall records.
# ------------------------------------------------------------

cumulative <- publication_history %>%
  count(Publication_Year, name = "n") %>%
  complete(
    Publication_Year = year_range,
    fill = list(n = 0)
  ) %>%
  arrange(Publication_Year) %>%
  mutate(Cumulative_Publications = cumsum(n))

# ------------------------------------------------------------
# 8. Shared theme
# ------------------------------------------------------------

theme_s2 <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0
    ),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(8, 12, 8, 8)
  )

# ------------------------------------------------------------
# 9. Panel A - annual publications
# ------------------------------------------------------------

pA <- ggplot(
  annual,
  aes(
    x = Publication_Year,
    y = n,
    fill = Evidence_Type
  )
) +
  geom_col(width = 0.9) +
  scale_fill_manual(
    values = c(
      "Whale" = "#4477AA",
      "Both"  = "#8B7D6B",
      "Wood"  = "#CC8844"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = seq(1880, 2020, by = 20),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "A  Annual contributing publications",
    x = NULL,
    y = "Number of publications",
    fill = "Evidence type"
  ) +
  theme_s2

# ------------------------------------------------------------
# 10. Panel B - cumulative evidence base
# ------------------------------------------------------------

final_n <- max(cumulative$Cumulative_Publications)

pB <- ggplot(
  cumulative,
  aes(
    x = Publication_Year,
    y = Cumulative_Publications
  )
) +
  geom_step(
    linewidth = 0.8
  ) +
  geom_point(
    data = cumulative %>%
      filter(Publication_Year == max(Publication_Year)),
    size = 2.5
  ) +
  annotate(
    "text",
    x = max(cumulative$Publication_Year),
    y = final_n,
    label = paste0(final_n, " publications"),
    hjust = 1.05,
    vjust = -0.8,
    size = 3.5
  ) +
  scale_x_continuous(
    breaks = seq(1880, 2020, by = 20),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(),
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    title = "B  Cumulative contributing publications",
    x = "Publication year",
    y = "Cumulative publications"
  ) +
  theme_s2 +
  theme(
    legend.position = "none"
  )

# ------------------------------------------------------------
# 11. Combine panels
# ------------------------------------------------------------

fig_s2 <- pA / pB +
  plot_layout(
    heights = c(1.15, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )

# ------------------------------------------------------------
# 12. Save
# ------------------------------------------------------------

dir.create(
  "Outputs/Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "Outputs/Figures/Supplementary_Figure_S2_data_source_history.png",
  fig_s2,
  width = 10,
  height = 7.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Outputs/Figures/Supplementary_Figure_S2_data_source_history.pdf",
  fig_s2,
  width = 10,
  height = 7.5,
  bg = "white"
)

# ------------------------------------------------------------
# 13. Export underlying publication-history table
# ------------------------------------------------------------

write.csv(
  publication_history %>%
    select(
      Source_ID,
      Publication_Year,
      Evidence_Type,
      Authors,
      Title,
      Journal_or_Repository,
      DOI
    ),
  "Outputs/Tables/supplementary_publication_history.csv",
  row.names = FALSE
)

cat("\nS2 complete.\n")
cat("Publications plotted:", nrow(publication_history), "\n")
cat(
  "Years:",
  min(publication_history$Publication_Year),
  "-",
  max(publication_history$Publication_Year),
  "\n"
)
cat("Final cumulative total:", final_n, "\n")