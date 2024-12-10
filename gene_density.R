

# plot the distribution of genes/TEs 
library(dplyr)
library(ggplot2)

setwd("C:/Dawei_UCB/1_UCB_Research/Aging/gene_TE_distribution/")

# Read the BED file (replace "TEs_Chr1.bed" with the actual file path)
bed1 <- read.table("TEs_Chr5.bed", header = FALSE, col.names = c("chromosome", "start", "stop", "name", "score", "strand"))
bed2 <- read.table("genes_Chr5.bed", header = FALSE, col.names = c("chromosome", "start", "stop", "name", "score", "strand"))

# Combine data from both files
combined_data <- bind_rows(
  mutate(bed1, file = "TEs_Chr5"),
  mutate(bed2, file = "genes_Chr5")
)

# Find the maximum stop site for each chromosome
max_stops <- combined_data %>%
  group_by(chromosome, file) %>%
  summarize(max_stop = max(stop))

# Create 10kbp bins and calculate the sum of rows in each bin
bin_size <- 250000
bed_data_bins <- combined_data %>%
  mutate(bin = (floor((stop - 1) / bin_size) + 1) * bin_size) %>%
  group_by(chromosome, file, bin) %>%
  summarize(sum_rows = n())

# Keep only the first data point for each bin
bed_data_first_point <- bed_data_bins %>%
  group_by(chromosome, file, bin) %>%
  summarize(sum_rows = sum(sum_rows)) %>%
  ungroup()

# Define custom colors
my_colors <- c("TEs_Chr5" = "#ca539f", "genes_Chr5" = "#3a4da1")

# Plot the data points and connect them with a line
ggplot(bed_data_first_point, aes(x = bin, y = sum_rows, color = file)) + theme_classic()+
  geom_line() +
  labs(x = "Bin", y = "Sum of Rows", color = "File") +
  ggtitle("Sum of Rows in 10kbp Bins (Connected Lines by File)") +
  scale_color_manual(values = my_colors) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M"),expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(
    name = "number of genes/TEs",
    breaks = seq(0, 300, by = 50),
    labels = seq(0, 300, by = 50),
    limits = c(0, 300), # Set the y-axis limits
    expand = expansion(mult = c(0, 0)) ) 















#distribution of number of unchanged cytosines (100%) in each bin 
library(dplyr)
library(ggplot2)

# List of your bed files
bed_files <- c(
  "sperm.bed", "embryo.bed", "D13.bed", "D18.bed",
  "D23.bed", "D30.bed", "D45.bed", "D60.bed"
)

# Function to read and preprocess bed files
read_and_preprocess_bed <- function(file) {
  # Read the BED file
  bed_data <- read.table(file, header = FALSE, col.names = c("chromosome", "start", "stop", "name", "score", "strand"))
  
  # Add an extra column (column 7) and set the value as 1 for all rows
  bed_data <- mutate(bed_data, extra_column = 1)
  
  # Extract chromosome 1 data
  chromosome1_data <- bed_data %>% filter(chromosome == "Chr1")
  
  # Add a file column with the current file name
  chromosome1_data$file <- file
  
  return(chromosome1_data)
}

# Read and preprocess chromosome 1 data from all bed files
all_chrom1_data <- lapply(bed_files, read_and_preprocess_bed) %>%
  bind_rows()

# Combine data from all files
combined_data <- all_chrom1_data

# Set the maximum stop site for each chromosome
max_stop_values <- c(30.43e6, 19.70e6, 23.46e6, 18.59e6, 26.97e6)
combined_data$max_stop <- max_stop_values[match(combined_data$chromosome, c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"))]

# Create 100 kbp bins and calculate the sum of column 7 in each bin
bin_size <- 250000
combined_data_bins <- combined_data %>%
  mutate(bin = (floor((start - 1) / bin_size) + 1) * bin_size) %>%
  group_by(bin, file) %>%
  summarize(sum_column7 = sum(extra_column)) %>%
  ungroup()

# Keep only the first data point for each bin
combined_data_first_point <- combined_data_bins %>%
  group_by(bin, file) %>%
  summarize(sum_column7 = sum(sum_column7)) %>%
  ungroup()

# Define custom colors for each dataset
my_colors <- setNames(
  c("sperm.bed" = "red", "embryo.bed" = "cyan3", "D13.bed" = "blue",
    "D18.bed" = "orange", "D23.bed" = "green", "D30.bed" = "purple",
    "D45.bed" = "pink", "D60.bed" = "brown"),
  bed_files
)

# Plot the data points and connect them with a line using chromosome location as x axis
ggplot(combined_data_first_point, aes(x = bin, y = sum_column7, color = file)) +
  geom_line(size = 1) +
  labs(x = "Chromosome Location", y = "Sum of Column 7", title = "Sum of Column 7 in 100kbp Bins for Chromosome 1") +
  theme_minimal() +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_color_manual(values = my_colors)













#plot the location of all data points in a number of files

library(dplyr)
library(ggplot2)

# List of your bed files
bed_files <- c(
  "sperm.bed", "embryo.bed", "D13.bed", "D18.bed",
  "D23.bed", "D30.bed", "D45.bed", "D60.bed"
)

# Function to read and preprocess bed files
read_and_preprocess_bed <- function(file) {
  # Read the BED file
  bed_data <- read.table(file, header = FALSE, col.names = c("chromosome", "start", "stop", "name", "score", "strand"))
  
  # Add an extra column (column 7) and set the value as 1 for all rows
  bed_data <- mutate(bed_data, extra_column = 1)
  
  # Extract chromosome 1 data
  chromosome1_data <- bed_data %>% filter(chromosome == "Chr1")
  
  # Add a file column with the current file name
  chromosome1_data$file <- file
  
  return(chromosome1_data)
}

# Read and preprocess chromosome 1 data from all bed files
all_chrom1_data <- lapply(bed_files, read_and_preprocess_bed) %>%
  bind_rows()

# Combine data from all files
combined_data <- all_chrom1_data

# Set the maximum stop site for each chromosome
max_stop_values <- c(30.43e6, 19.70e6, 23.46e6, 18.59e6, 26.97e6)
combined_data$max_stop <- max_stop_values[match(combined_data$chromosome, c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"))]

# Define custom colors for each dataset
my_colors <- setNames(
  c("sperm.bed" = "red", "embryo.bed" = "cyan3", "D13.bed" = "blue",
    "D18.bed" = "orange", "D23.bed" = "green", "D30.bed" = "purple",
    "D45.bed" = "pink", "D60.bed" = "brown"),
  bed_files
)

# Plot all data points individually using chromosome location as x-axis
ggplot(combined_data, aes(x = start, y = score, color = file)) +
  geom_point(size = 2) +
  labs(x = "Chromosome Location", y = "Score (row 6)", title = "Individual Data Points for Chromosome 1") +
  theme_minimal() +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_color_manual(values = my_colors)











setwd("C:/Dawei_UCB/1_UCB_Research/Plant_Epigenetic_Clock/gene_TE_distribution/stable_cytosine_distribution/")

#plot the location of all data points in a single file


library(dplyr)
library(ggplot2)

# Read the BED file
bed_data <- read.table("unchanged_cytosine.bed", header = FALSE, 
                       col.names = c("chromosome", "start", "stop", "unmethy", "methy", "value"))

# Filter data for chromosome 1
chromosome1_data <- bed_data %>% filter(chromosome == "Chr3")

# Set the maximum stop site for chromosome 1
max_stop_chr1 <- 23.46e6
chromosome1_data$max_stop <- max_stop_chr1

# Define custom color
my_color <- "blue"

# Plot all data points individually using chromosome location as x-axis
ggplot(chromosome1_data, aes(x = start, y = value, color = my_color)) +
  geom_point(size = 1) +
  labs(x = "Chromosome Location", y = "CG methylation", title = "Individual Data Points for Chromosome 1") +
  theme_classic() +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M"), 
                     limits = c(0, max(chromosome1_data$max_stop, na.rm = TRUE)), 
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(90, 110), expand = expansion(mult = c(0, 0))) +
  scale_color_manual(values = my_color)



#plot distribution of data

setwd("C:/Dawei_UCB/1_UCB_Research/Plant_Epigenetic_Clock/DMCs")

library(ggplot2)
library(dplyr)

# Read the BED file
DMCs_nochange <- as.data.frame(read.table("DMCs_nochange_TEs.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

colnames(DMCs_nochange) <- c("chromosome", "start", "end", "methylation")

# Create a scatter plot for each chromosome

ggplot(DMCs_nochange, aes(x = start, y = methylation)) +
  geom_point() +
  labs(x = "Chromosome Location", y = "Methylation") +
  ggtitle("Scatter Plot of Methylation vs. Chromosome Location") +
  facet_grid(chromosome ~ ., scales = "free_x", space = "free_x") +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  theme_classic()

