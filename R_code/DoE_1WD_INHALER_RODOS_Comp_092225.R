#Dispersion_EMD_exercise
#updated 09-22-25

getwd()
setwd("/Users/grace/Inhaled Biologics DoE/LD")
install.packages("transport")
install.packages("tidyverse")
install.packages("janitor")
library(transport)
library(stringr)
library(tidyverse)
library(janitor)

#Started by creating nested folder system in my wd (Inhaled Biologics DoE ->
#folders named after each formulation (e.g., 402040_IMT, 311059_IMT, etc) ->
#within each folder, had sub-folders named "Rep_1", "Rep_2", "Rep_3" (e.g., "Run4_baseline_INHALER_09072025_rep1")
#within each of those, RODOS and INHALER cvs files named "RunX_INHALER_repY" and "RunX_RODOS_repY"
# List all CSV files in working directory, searching recursively through subdirectories
file_paths <- list.files(
  path = ".",
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# Read all files, skipping the initial rows and selecting the data columns
combined_data <- read_csv(
  file_paths,
  id = "source_file",
  skip = 2, # <--- Change this number based on your file inspection
  col_types = cols(.default = "c")
) %>%
  # Clean column names
  clean_names() %>%
  mutate(
    xo_mm = as.numeric(xo_mm),
    q3_percent = as.numeric(q3_percent)
  ) %>%
  filter(!is.na(xo_mm)) %>%
  
  mutate(
    # Extract formulation using the `Run` number from the filename
    Formulation = str_extract(source_file, "Run\\d+"),
    # Extract replicate number from filename
    Replicate = str_extract(source_file, "Rep_\\d+"),
    # Extract module from the filename
    Module = str_extract(source_file, "INHALER|RODOS")
  )


# Now, plot the CDFs
ggplot(combined_data, aes(x = xo_mm, y = q3_percent, color = Module, linetype = as.factor(Replicate))) +
  geom_step(linewidth = 0.25) +
  facet_wrap(~ Formulation) + # Separate plots for each formulation
  scale_x_log10() +
  labs(title = NULL,
       x = "Particle Size (log(µm))", y = "Cumulative Distribution (%)", linetype = NULL) +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold")) +
  scale_linetype_discrete(labels = c("Rep_1" = "Rep 1",
                                     "Rep_2" = "Rep 2",
                                     "Rep_3" = "Rep 3"))

# Data for INHALER in Run7
inhaler_run7_q3 <- combined_data %>%
  filter(Formulation == "Run7", Module == "INHALER") %>%
  pull(q3_percent) # Extract just the q3_percent values

# Data for RODOS in Run7
rodos_run7_q3 <- combined_data %>%
  filter(Formulation == "Run7", Module == "RODOS") %>%
  pull(q3_percent) # Extract just the q3_percent values

# Compute 1-Wasserstein distance for Run7
wd1_run7 <- wasserstein1d(inhaler_run7_q3, rodos_run7_q3, p = 1)
print(paste("1-Wasserstein distance for Run7 (INHALER vs RODOS):", wd1_run7))
wd1_run7

# Data for INHALER in Run3
inhaler_run3_q3 <- combined_data %>%
  filter(Formulation == "Run3", Module == "INHALER") %>%
  pull(q3_percent) # Extract just the q3_percent values

# Data for RODOS in Run3
rodos_run3_q3 <- combined_data %>%
  filter(Formulation == "Run3", Module == "RODOS") %>%
  pull(q3_percent) # Extract just the q3_percent values

# Compute 1-Wasserstein distance for Run3
wd1_run3 <- wasserstein1d(inhaler_run3_q3, rodos_run3_q3, p = 1)
print(paste("1-Wasserstein distance for Run3 (INHALER vs RODOS):", wd1_run3))
wd1_run3

# Data for INHALER in Run2
inhaler_run2_q3 <- combined_data %>%
  filter(Formulation == "Run2", Module == "INHALER") %>%
  pull(q3_percent) # Extract just the q3_percent values

# Data for RODOS in Run2
rodos_run2_q3 <- combined_data %>%
  filter(Formulation == "Run2", Module == "RODOS") %>%
  pull(q3_percent) # Extract just the q3_percent values

# Compute 1-Wasserstein distance for Run2
wd1_run2 <- wasserstein1d(inhaler_run2_q3, rodos_run2_q3, p = 1)
print(paste("1-Wasserstein distance for Run2 (INHALER vs RODOS):", wd1_run2))
wd1_run2

# Data for INHALER in Run4
inhaler_run4_q3 <- combined_data %>%
  filter(Formulation == "Run4", Module == "INHALER") %>%
  pull(q3_percent) # Extract just the q3_percent values

# Data for RODOS in Run4
rodos_run4_q3 <- combined_data %>%
  filter(Formulation == "Run4", Module == "RODOS") %>%
  pull(q3_percent) # Extract just the q3_percent values

# Compute 1-Wasserstein distance for Run4
wd1_run4 <- wasserstein1d(inhaler_run4_q3, rodos_run4_q3, p = 1)
print(paste("1-Wasserstein distance for Run4 (INHALER vs RODOS):", wd1_run4))
wd1_run4

# Data for INHALER in Run6
inhaler_run6_q3 <- combined_data %>%
  filter(Formulation == "Run6", Module == "INHALER") %>%
  pull(q3_percent) # Extract just the q3_percent values

# Data for RODOS in Run6
rodos_run6_q3 <- combined_data %>%
  filter(Formulation == "Run6", Module == "RODOS") %>%
  pull(q3_percent) # Extract just the q3_percent values

# Compute 1-Wasserstein distance for Run6
wd1_run6 <- wasserstein1d(inhaler_run6_q3, rodos_run6_q3, p = 1)
print(paste("1-Wasserstein distance for Run6 (INHALER vs RODOS):", wd1_run6))
wd1_run6

# Data for INHALER in Run8
inhaler_run8_q3 <- combined_data %>%
  filter(Formulation == "Run8", Module == "INHALER") %>%
  pull(q3_percent) # Extract just the q3_percent values

# Data for RODOS in Run8
rodos_run8_q3 <- combined_data %>%
  filter(Formulation == "Run8", Module == "RODOS") %>%
  pull(q3_percent) # Extract just the q3_percent values

# Compute 1-Wasserstein distance for Run8
wd1_run8 <- wasserstein1d(inhaler_run8_q3, rodos_run8_q3, p = 1)
print(paste("1-Wasserstein distance for Run8 (INHALER vs RODOS):", wd1_run8))
wd1_run8

# Data for INHALER in Run9
inhaler_run9_q3 <- combined_data %>%
  filter(Formulation == "Run9", Module == "INHALER") %>%
  pull(q3_percent) # Extract just the q3_percent values

# Data for RODOS in Run9
rodos_run9_q3 <- combined_data %>%
  filter(Formulation == "Run9", Module == "RODOS") %>%
  pull(q3_percent) # Extract just the q3_percent values

# Compute 1-Wasserstein distance for Run9
wd1_run9 <- wasserstein1d(inhaler_run9_q3, rodos_run9_q3, p = 1)
print(paste("1-Wasserstein distance for Run9 (INHALER vs RODOS):", wd1_run9))
wd1_run9