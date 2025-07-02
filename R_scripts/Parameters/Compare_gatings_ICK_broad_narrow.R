# Taget FCS filename
target_fcs <- "ICK-043-02_Unmixed.fcs"

# Custom name for the workspace (keep identifier as short as possible)
custom_name <- c("Lea", "Jonas", "Samu", "Broad", "Narrow")

# User provided dictionary for translating labels
# This dictionary should map the original labels to the desired labels
label_dict <- c(
  "root" = "Removed", #"Erys"
  
  "/no erys" = "Removed", #"Doublets"
  "/No Erythrocytes" = "Removed", #"Doublets"
  "/No_erys" = "Removed", #"Doublets"
  
  "/no erys/Single Cells" = "Removed", #"Myelos"
  "/No Erythrocytes/Single Cells" = "Removed", #"Myelos"
  "/No_erys/singlets"  = "Removed", #"Myelos"
  
  "/no erys/Single Cells/lymphocytes" = "Removed", #"ncLymphos"
  "/No Erythrocytes/Single Cells/Lymphocytes"  = "Removed", #"ncLymphos"
  "/No_erys/singlets/lymphocytes"  = "Removed", #"ncLymphos"
  
  "/no erys/Single Cells/lymphocytes/CD3+ TCRyd-" = "ncT cells",
  "/no erys/Single Cells/lymphocytes/CD3+ TCRyd-/CD3+ CD56-"  = "ncT cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+gdTCR-" = "ncT cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+gdTCR-/CD3+CD56-" =
    "ncT cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+ TCRgd-" = "ncT cells",
  "/No_erys/singlets/lymphocytes/CD3+/gd-/TCs (CD3+CD56-)" = "ncT cells",
  "/No_erys/singlets/lymphocytes/CD3+/gd-" = "ncT cells",
  "/No_erys/singlets/lymphocytes/CD3+" = "ncT cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+ TCRgd-/CD3+ CD56-" =
    "ncT cells",
  
  "/no erys/Single Cells/lymphocytes/CD3+ TCRyd-/CD3+ CD56-/CD4 T cells" =
    "CD4T cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+gdTCR-/CD3+CD56-/CD4+" =
    "CD4T cells",
  "/No_erys/singlets/lymphocytes/CD3+/gd-/TCs (CD3+CD56-)/CD4TCs" =
    "CD4T cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+ TCRgd-/CD3+ CD56-/CD4+" =
    "CD4T cells",
  
  "/no erys/Single Cells/lymphocytes/CD3+ TCRyd-/CD3+ CD56-/CD8 Tcells" =
    "CD8T cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+gdTCR-/CD3+CD56-/CD8+" =
    "CD8T cells",
  "/No_erys/singlets/lymphocytes/CD3+/gd-/TCs (CD3+CD56-)/CD8TCs" =
    "CD8T cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+ TCRgd-/CD3+ CD56-/CD8+" =
    "CD8T cells",
  
  "/no erys/Single Cells/lymphocytes/CD3+ TCRyd-/NKTcells" = "NKT cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+gdTCR-/CD3+CD56+" =
    "NKT cells",
  "/No_erys/singlets/lymphocytes/CD3+/gd-/NKTCs" = "NKT cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+ TCRgd-/CD3+ CD56+ " =
    "NKT cells",
  
  "/no erys/Single Cells/lymphocytes/yd Tcells" = "gdT cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+ TCRgd+" = "gdT cells",
  "/No_erys/singlets/lymphocytes/CD3+/gdTCs" = "gdT cells",
  "/No_erys/singlets/lymphocytes/gdTCs" = "gdT cells",
  
  "/no erys/Single Cells/lymphocytes/CD3-" = "ncCD3-",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3-" = "ncCD3-",
  "/No_erys/singlets/lymphocytes/CD3-" = "ncCD3-",
  "/No_erys/singlets/lymphocytes/CD3-/CD19-CD56-" = "ncCD3-",
  "/No_erys/singlets/lymphocytes/CD3-/CD19-CD56-/CD3-CD56-CD19-CD11b-" =
    "ncCD3-",
  "/No_erys/singlets/lymphocytes/CD3-/CD19-CD56-/DCs" = "ncCD3-",
  
  "/no erys/Single Cells/lymphocytes/CD3-/B cells" = "B cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3-/CD19+ CD56-" = "B cells",
  "/No_erys/singlets/lymphocytes/CD3-/BCs" = "B cells",
  
  "/no erys/Single Cells/lymphocytes/CD3-/NK cells" = "NK cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3-/CD19- CD56+" = "NK cells",
  "/No_erys/singlets/lymphocytes/CD3-/NKCs" = "NK cells"
)