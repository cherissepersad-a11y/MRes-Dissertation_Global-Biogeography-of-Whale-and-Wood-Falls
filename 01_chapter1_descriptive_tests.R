# 01_chapter1_descriptive_tests.R
# Chapter 1: evidence-base description, depth comparison and supporting tests.
source("R/00_setup_and_validate.R")

# Restrict event-level comparisons to records explicitly included for analysis.
events_analysis <- events %>%
  filter(Analysis_Include == 1, Fall_Type %in% c("Whale", "Wood")) %>%
  mutate(Depth_m = case_when(
    !is.na(Event_Depth_Min_m) & !is.na(Event_Depth_Max_m) ~ (Event_Depth_Min_m + Event_Depth_Max_m)/2,
    !is.na(Event_Depth_Min_m) ~ Event_Depth_Min_m,
    !is.na(Event_Depth_Max_m) ~ Event_Depth_Max_m,
    TRUE ~ NA_real_
  ))

# Summaries are reported before hypothesis tests so sample size and missingness
# remain visible. This prevents a statistically significant p-value from being
# mistaken for a large ecological effect.
depth_summary <- events_analysis %>%
  group_by(Fall_Type) %>%
  summarise(n_events = n(), n_depth = sum(!is.na(Depth_m)),
            median_depth_m = median(Depth_m, na.rm = TRUE),
            IQR_depth_m = IQR(Depth_m, na.rm = TRUE), .groups = "drop")
write_csv(depth_summary, "Outputs/Tables/ch1_depth_summary.csv")

# The event-depth distributions are non-normal and strongly unbalanced in
# sample size, so a Wilcoxon/Mann-Whitney comparison is used as a descriptive
# two-group test. Interpret the effect size and distributions, not p alone.
depth_test <- wilcox.test(Depth_m ~ Fall_Type, data = events_analysis,
                          exact = FALSE, conf.int = TRUE)
writeLines(capture.output(depth_test), "Outputs/Tables/ch1_depth_wilcoxon.txt")

# Research-footprint tables explicitly describe sampling effort, not habitat
# prevalence. A high count therefore means "well represented in the database",
# not "ecologically more common".
node_basin <- site_universe %>% count(Fall_Type, Basin, name = "Analytical_nodes")
origin_counts <- events_analysis %>% count(Fall_Type, Origin, name = "Events")
write_csv(node_basin, "Outputs/Tables/ch1_nodes_by_basin.csv")
write_csv(origin_counts, "Outputs/Tables/ch1_events_by_origin.csv")
