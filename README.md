# MRes Organic-Fall Dissertation Reproducibility Package

## What this package is
This package is designed to reproduce the dissertation analyses from the frozen `MRes_Dissertation_Database_FINAL_AnalysisReady` database. R is the analytical authority. Gephi files are supplied for interactive exploration and figure prototyping, but Gephi is not required to reproduce the statistical results.

## Scientific interpretation
A site-to-site network edge means that two analytical nodes share recorded taxa. A bipartite site-to-taxon edge means that a taxon was recorded at that analytical node. Neither edge demonstrates dispersal, gene flow or a realised migration corridor. Throughout the workflow, use terms such as *compositional association*, *shared taxon*, *cross-module taxon* or *taxon contributing to compositional connectivity*.

## Folder structure
- `Database/` — frozen final database and a nested ZIP of the same files. Do not edit these files.
- `R/` — fully commented R scripts.
- `Gephi/nodes/` and `Gephi/edges/` — ready-to-import CSVs at species/OTU, genus and family resolution for Whale, Wood and Combined scopes.
- `Gephi/guide/` — Gephi instructions.
- `Outputs/` — destination for regenerated figures, tables and diagnostics.
- `Documentation/` — analysis crosswalk and interpretation notes.

## Software
Use R 4.3 or later where possible. Required packages are listed in `R/00_setup_and_validate.R`. The script does not install packages silently. Install missing packages yourself, then rerun. Gephi is optional and should be a current stable desktop release.

## Reproduce the analysis
1. Unzip this package without changing the internal folder names.
2. Open R or RStudio and set the working directory to the package root. In RStudio, `Session > Set Working Directory > Choose Directory` is sufficient.
3. Run `source("R/00_setup_and_validate.R")`. Inspect `Outputs/Diagnostics/validation_summary.csv`. Stop if validation fails.
4. For the complete core workflow, run `source("R/99_run_all.R")`.
5. Run `R/07_network_nulls_modules_sensitivity.R` separately when you are ready for the computationally heavier null-model and threshold-sensitivity analyses. The 999-permutation example is deliberately commented out so it cannot start accidentally.
6. Compare regenerated outputs with the dissertation figures/tables. `Outputs/Diagnostics/sessionInfo.txt` records the software environment.

## Script order
`00` setup/validation → `01` Chapter 1 tests → `02` Chapter 1 figures → `03` overlap/beta diversity → `04` distance decay/multivariate template → `05` site-similarity networks → `06` taxon-resolved bipartite/connectors → `07` nulls/modules/sensitivity → `08` Gephi exports.

## Why there are two network types
The site-similarity network answers **which communities have similar recorded composition?** The bipartite network answers **which named taxa are responsible for the observed sharing?** The second is essential for the biological story requested for this dissertation. The same construction is repeated at species/OTU, genus and family levels.

## Taxon labels in figures
The R bipartite figure does not manually cherry-pick attractive names. The example labels taxa that span at least two basins and at least three analytical nodes. Those rules are explicit in `06_network_bipartite_taxon_connectors.R` and can be sensitivity-tested. If a different rule is used in the final thesis, record it in the caption and Methods.

## Zeros
Zeros in the presence matrices mean **taxon not recorded**, not confirmed ecological absence. This affects all incidence-based similarity and beta-diversity interpretation.
