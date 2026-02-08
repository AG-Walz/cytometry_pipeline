# Taget FCS filename
target_fcs <- samp_name_list[[1]][1]

# Custom name for the workspace (keep identifier as short as possible)
custom_name <- custom_name

# User provided dictionary for translating labels
# This dictionary should map the original labels to the desired labels
label_dict <- c(
  "root" = "Removed", #"Erys"
  "/no erys" = "Removed", #"Doublets"
  "/no erys/Single Cells" = "Removed", #"Myelos"
  "/no erys/Single Cells/CD45+" = "Removed",

  "/no erys/Single Cells/lymphocytes" = "Lymphos", #"ncLymphos"
  "/no erys/Single Cells/lymphocytes/CD45+" = "Lymphos",
  "/no erys/Single Cells/CD45+/lymphocytes" = "Lymphos",

  "/no erys/Single Cells/lymphocytes/CD3+" = "CD3+",
  "/no erys/Single Cells/lymphocytes/CD45+/CD3+"  = "CD3+",
  "/no erys/Single Cells/CD45+/lymphocytes/CD3+" = "CD3+",

  "/no erys/Single Cells/lymphocytes/CD3+/T cells"  = "T cells",
  "/no erys/Single Cells/lymphocytes/CD45+/CD3+/T cells"  = "T cells",
  "/no erys/Single Cells/CD45+/lymphocytes/CD3+/T cells"  = "T cells",

  "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD4+" = "CD4+ T cells",
  "/no erys/Single Cells/lymphocytes/CD45+/CD3+/T cells/CD4+" = "CD4+ T cells",
  "/no erys/Single Cells/CD45+/lymphocytes/CD3+/T cells/CD4+" = "CD4+ T cells",

  "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD4+/FAS+" = "FAS+CD4+ T cells",
  "/no erys/Single Cells/lymphocytes/CD45+/CD3+/T cells/CD4+/FAS+" = "FAS+CD4+ T cells",
  "/no erys/Single Cells/CD45+/lymphocytes/CD3+/T cells/CD4+/FAS+" = "FAS+CD4+ T cells",

  "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD4+/PD1+" = "PD1+CD4+ T cells",
  "/no erys/Single Cells/lymphocytes/CD45+/CD3+/T cells/CD4+/PD1+" = "PD1+CD4+ T cells",
  "/no erys/Single Cells/CD45+/lymphocytes/CD3+/T cells/CD4+/PD1+" = "PD1+CD4+ T cells",

  "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD4+/FAS+ PD1+" = "FAS+PD1+CD4+ T cells",
  "/no erys/Single Cells/lymphocytes/CD45+/CD3+/T cells/CD4+/FAS+ PD1+" = "FAS+PD1+CD4+ T cells",
  "/no erys/Single Cells/CD45+/lymphocytes/CD3+/T cells/CD4+/FAS+ PD1+" = "FAS+PD1+CD4+ T cells",

  "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD8+" = "CD8+ T cells",
  "/no erys/Single Cells/lymphocytes/CD45+/CD3+/T cells/CD8+" = "CD8+ T cells",
  "/no erys/Single Cells/CD45+/lymphocytes/CD3+/T cells/CD8+" = "CD8+ T cells",

  "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD8+/FAS+" = "FAS+CD8+ T cells",
  "/no erys/Single Cells/lymphocytes/CD45+/CD3+/T cells/CD8+/FAS+" = "FAS+CD8+ T cells",
  "/no erys/Single Cells/CD45+/lymphocytes/CD3+/T cells/CD8+/FAS+" = "FAS+CD8+ T cells",

  "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD8+/PD1+" = "PD1+CD8+ T cells",
  "/no erys/Single Cells/lymphocytes/CD45+/CD3+/T cells/CD8+/PD1+" = "PD1+CD8+ T cells",
  "/no erys/Single Cells/CD45+/lymphocytes/CD3+/T cells/CD8+/PD1+" = "PD1+CD8+ T cells",

  "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD8+/FAS+ PD1+" = "FAS+PD1+CD8+ T cells",
  "/no erys/Single Cells/lymphocytes/CD45+/CD3+/T cells/CD8+/FAS+ PD1+" = "FAS+PD1+CD8+ T cells",
  "/no erys/Single Cells/CD45+/lymphocytes/CD3+/T cells/CD8+/FAS+ PD1+" = "FAS+PD1+CD8+ T cells",


  "/no erys/Single Cells/lymphocytes/CD3-" = "CD3-",
  "/no erys/Single Cells/lymphocytes/CD45+/CD3-" = "CD3-",
  "/no erys/Single Cells/CD45+/lymphocytes/CD3-" = "CD3-"

)
