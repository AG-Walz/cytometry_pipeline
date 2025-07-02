# devtools::install_github("Close-your-eyes/fcexpr")
library(tidyverse)
library(fcexpr)

ps_row <- fcexpr::wsx_get_popstats(ws = "/home/ngr24/Desktop/Maerklin/cytometry_pipeline/FlowJo tests/ICK-049-051.wsp" )
View(ps_row)
ps <- ps_row[["counts"]]
View(ps)
openxlsx::write.xlsx(ps, file ="/home/ngr24/Desktop/Maerklin/cytometry_pipeline/FlowJo tests/ps.xlsx")

# Add metadata (xlsx with file names & metadata columns--> cohort file for ICK)
file_path <- "/home/ngr24/Desktop/Maerklin/cytometry_pipeline/FlowJo tests/"
fcs_files <- list.files(path = file_path, pattern = "fcs", recursive = T, full.names= T)
file_name <- basename(fcs_files)
table <- data.frame(fcs_files, file_name)
openxlsx::write.xlsx(table, file = paste0(file_path, "metadata_table.xlsx"))


# Filter opulations
for(unique_pop in unique(Population)){
  filtered_pop[unique_pop] <-dplyr::filter(ps_row[["counts"]], Population==unique_pop)
}
ps_emtregs <- dplyr::filter(ps_row[["counts"]], Population=="em Treg")
ps_emtregs <- dplyr::filter(ps_row[["counts"]], Population=="effector")
ps_emtregs <- dplyr::filter(ps_row[["counts"]], Population=="memory")
ps_emtregs <- dplyr::filter(ps_row[["counts"]], Population=="cm Treg")
ps_emtregs <- dplyr::filter(ps_row[["counts"]], Population=="em Treg")
