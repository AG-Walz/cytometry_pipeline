source('Compare.R', chdir = TRUE)
ensure_and_load_packages(c("readxl", "ggplot2", "dplyr", "stringr", "Polychrome"))

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# Step 1: Read the Excel file
# Replace 'file_path.xlsx' with the actual path to your Excel file
data <- read_excel("Metadaten ICK-018-xx (1).xlsx")

# Step 2: Extract columns of interest (Proben-ID, T, CD4+(counts), CD8+(counts))
# Rename columns if they have spaces or special characters
data <- data %>% 
  select("Proben-ID", "CD4+ (count)", "CD8+ (count)")

# Step 3: Calculate CD4/CD8 ratio
data <- data %>%
  mutate(Patient_ID = stringr::str_sub(`Proben-ID`, 5, 7),
         Timepoint = as.numeric(stringr::str_sub(`Proben-ID`, -2, -1)))%>%  # Convert timepoint to numeric
  mutate(CD8CD4_ratio = `CD8+ (count)`/`CD4+ (count)` ) %>%
  mutate(CD4CD8_ratio = `CD4+ (count)`/`CD8+ (count)` ) %>%
  filter(!is.na(`CD4+ (count)`) & !is.na(`CD8+ (count)`))  # Remove rows with missing CD4+ or CD8+ counts

# Adjust the Glasbey color palette to match the number of unique patients
glasbey_colors <-Polychrome::glasbey.colors(length(unique(data$Patient_ID)))
custom_palette <- c("darkred","red","darkorange","gold", "greenyellow","limegreen","darkgreen",  "#1F9698",
                    "navy","cyan2", "#009FFF", "deeppink", "purple",  "chocolate4", "black")

# Assign colors explicitly to the levels of Patient_ID
color_mapping <- setNames(custom_palette, levels(data$Patient_ID))

data_edit <- data[c(1:16,20:23,26:29,31:38,42:53),]
# Step 4: Plot the CD4/CD8 ratio for each Proben-ID over time
theme_set(theme_bw())
plot<-ggplot(data_edit, aes(x = Timepoint, y = CD8CD4_ratio, color = `Patient_ID`)) +
  geom_line(size=0.8)+
  labs(x = "Timepoint", y = "CD8/CD4 Cell Count Ratio")+
  scale_color_manual(values = color_mapping) + 
  theme(plot.title = element_text(size = 14,   hjust = 0.5),# face = "bold",
        
        # Customize axis labels and text sizes
        axis.title.x = element_text(vjust = 0, size = 12),
        axis.title.y = element_text(vjust = 2, size = 12),
        axis.text = element_text( size = 10),
        
        # Customize the legend: remove title and adjust text size and key size
        legend.title = element_text(size = 10, hjust = 0.5), 
        legend.title.align = 0.5,
        legend.text = element_text(size = 10),  # Adjust legend text size
        legend.key.size = unit(0.5, "cm")
        ) +  # Adjust the size of legend keys (ticks)
  guides(color = guide_legend(ncol = 1))
ggsave("CD8CD4_ratio.pdf", plot = plot, width = 4.8, height = 3.6, dpi = 400)

plot<-ggplot(data_edit, aes(x = Timepoint, y = CD4CD8_ratio, color = `Patient_ID`)) +
  geom_line(size=0.8)+
  labs(x = "Timepoint", y = "CD4/CD8 Cell Count Ratio")+
  scale_color_manual(values = color_mapping) + 
  theme(plot.title = element_text(size = 14,   hjust = 0.5),# face = "bold",
        
        # Customize axis labels and text sizes
        axis.title.x = element_text(vjust = 0, size = 12),
        axis.title.y = element_text(vjust = 2, size = 12),
        axis.text = element_text( size = 10),
        
        # Customize the legend: remove title and adjust text size and key size
        legend.title = element_text(size = 10, hjust = 0.5), 
        legend.title.align = 0.5,
        legend.text = element_text(size = 10),  # Adjust legend text size
        legend.key.size = unit(0.5, "cm")
  ) +  # Adjust the size of legend keys (ticks)
  guides(color = guide_legend(ncol = 1))
ggsave("CD4CD8_ratio.pdf", plot = plot, width = 4.8, height = 3.6, dpi = 400)




