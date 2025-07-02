### ============================ GENERATE BARPLOTS =========================== ###
# This script allows to generate barplots from a list of flow_frames or a matrix

source('Compare.R', chdir = TRUE)

ensure_and_load_packages(c("readxl", "ggplot2", "dplyr", "stringr", "Polychrome"))

# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# =============================== Data preparation  =============================== #
# Replace 'file_path.xlsx' with the actual path to your Excel file
data <- read_excel("Metadaten ICK-018-xx (1).xlsx")

df_subset <- data %>%
  mutate(
    # Extract the patient number from characters 5 to 7
    PatientNumber = str_sub(`Proben-ID`, 5, 7),
    
    # Extract the timepoint from characters 9 to 10
    Timepoint = str_sub(`Proben-ID`, 10),
    
    `T-cells(%)` = `T-cells(%)`*`Lymphocytes (%)`,
    `B cells(%)` = `B cells(%)`*`Lymphocytes (%)`,
    `NK cells (%)` = `NK cells (%)`*`Lymphocytes (%)`,
    `gd T cells (%)` = `T-cells(%)`*`gd T cells (%)`,
    `CD4+ (%)` = `T-cells(%)`*`CD4+ (%)`,
    `CD8+(%)` = `T-cells(%)`*`CD8+(%)`
  ) %>% 
  select(PatientNumber, Timepoint, `Granulocytes (%)`,`Monocytes (%)`, `Dentritic cells (%)`, 
         `B cells(%)`, `NK cells (%)`, 
         `gd T cells (%)`, `CD4+ (%)`, `CD8+(%)`)%>%
  na.omit()
  
# Reshape data to long format using tidyr::pivot_longer
df_long <- df_subset %>%
  pivot_longer(cols = c(`Granulocytes (%)`,`Monocytes (%)`, `Dentritic cells (%)`,`B cells(%)`, `NK cells (%)`, 
                        `gd T cells (%)`, `CD4+ (%)`, `CD8+(%)`), names_to = "Variable", values_to = "Value")


# =============================== Plotting === =============================== #

# Custom x-axis labels (you can define this based on your data)
custom_labels <- c("18.2","18.3","18.4","19.2","19.3","19.4","20.1","20.2","20.3","20.4",
                   "21.1","21.2","21.3","21.4","22.2","22.3","22.4","23.1","23.2","23.3","23.4","24.3",
                    "24.4","25.1","25.2","25.3","25.4","26.1","27.1","27.2","27.3",
                    "27.4","28.1","28.2","28.3","28.4","29.1","29.2","29.3","30.1","30.2","30.3","30.4","31.1",
                    "31.2","31.3","31.4","32.1","32.2","32.3","32.4"
) # Example labels

# Custom colors for each variable (change the hex codes as needed)
custom_colors <- c(`Monocytes (%)`="red",`Dentritic cells (%)`="orange",`B cells(%)`="purple", 
                   `NK cells (%)`="limegreen",`gd T cells (%)`="gold",  
                   `CD4+ (%)`="navy",  `CD8+(%)`= "#009FFF", `Granulocytes (%)`="darkgreen")

# Create the stacked barplot using ggplot2
theme_set(theme_bw())
ggplot(df_long, aes(x = interaction(Timepoint, PatientNumber), y = Value, fill = Variable)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  scale_x_discrete(labels = custom_labels) +   # Rename x-axis tick labels
  theme_minimal() +
  labs(
    x = "Patient.Timepoint",
    y = "Value",
    fill = "Cell Population",
    title = "Stacked Barplot by Patient Number and Timepoint"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Custom colors for each variable (change the hex codes as needed)
custom_colors <- c(`Monocytes (%)`="red",`Dentritic cells (%)`="orange",`B cells(%)`="purple", 
                   `NK cells (%)`="limegreen",`gd T cells (%)`="gold",  
                   `CD4+ (%)`="navy",  `CD8+(%)`= "#009FFF")


# Create the stacked barplot using ggplot2
theme_set(theme_bw())
ggplot(df_long, aes(x = interaction(Timepoint, PatientNumber), y = Value, fill = Variable)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  scale_x_discrete(labels = custom_labels) +   # Rename x-axis tick labels
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1.0)) + 
  theme_minimal() +
  labs(
    x = "Patient.Timepoint",
    y = "Value",
    fill = "Cell Population",
    title = "Stacked Barplot by Patient Number and Timepoint"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
