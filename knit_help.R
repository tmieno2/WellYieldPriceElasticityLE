
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
 knitr,
 here,
 rmarkdown,
 officedown,
 officer,
 data.table,
 tidyverse,
 flextable,
 stringr,
 sf,
 lfe,
 modelsummary,
 patchwork
)

render(here("Manuscript.Rmd"))
