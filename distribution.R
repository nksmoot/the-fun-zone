# plot the distribution of genes/TEs/DMCs 
#---------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
#---------------------------------------------------------------------------------
### User set variables
path <- "/Volumes/Extreme\ SSD/Williams_Lab/EM_seq/may_2024_resequenced/dmrs without Ws/split_chromosomes/"
output_dir <- "/Users/nksmoot/Desktop/Grad/r studio things/dmr_dist"
samples <- c('4-1_R0', '4-2_R0', '4-3_R0', '4-4_R0', '4-11_R0')
#my_colors <- c("4-1_R0" = "#6CDDEC", "4-2_R0" = "#5DC4E4", 
#               "4-3_R0" = "#4DABDA", "4-4_R0" = "#4092D3", "4-11_R0" = "#3279CB")

my_colors <- c("genes" = "#C85133", "TEs" = "#D9AD4B")
bin_size <- 350000
title <- "Genes TEs"
context <- "CHH"

#---------------------------------------------------------------------------------
chrom <- paste0("Chr", 1:5)
chrom

#---------------------------------------------------------------------------------
#Import files
#file names need to be in the format "samplename_all_CpG_Chr1.bed"
#be careful about this function, ken wrote it and it's buggy 
#main problem: this function will import ur samples ALPHABETICALLY 
#and then rename them the order of ur sample list, so the wrong name will be associated with the imported data
#to solve this, we will order sample names, apply ordered names to the ordered files, then reorder to the order initially set in samples
samples_sort <- sort(samples)
samples_sort


#from whole map methylation plots (i think this is the closest)
import_data <- function(path, tissue, chromosome) {
  bed_paths <- list.files(path, recursive = T, 
                          pattern = paste0("^(", paste(tissue, collapse = "|"), 
                                           ")", "_v_WT_DMRsBins", context, "_without_Ws_", chromosome, ".bed"), 
                          full.names = T)
  data <- lapply(bed_paths, read.table, header = F, 
                 col.names = c("chromosome", "start", "stop", "length", "star", 
                               "six", "seven","eight", "nine", "ten", "eleven", "twelve", 
                               "thirteen","fourteen", "fifteen", "sixteen"), 
                 colClasses = c("character", "numeric", "numeric", rep("NULL", 13)))
  
  names(data) <- tissue
  return(data)
}

samples_sort
all_data <- lapply(chrom, import_data, path = path, tissue = samples_sort)
names(all_data) <- chrom
#head(all_data)


#add genes and TEs
all_data <- list()
for (chr in chrom) {
  path_genes <- paste0("/Volumes/Extreme SSD/Williams_Lab/EM_seq/may_2024_resequenced/dmc_by_chr/split_chromosomes/genes_", chr, ".bed")
  all_data[[chr]]$genes <- read.table(path_genes, 
                                      col.names = c("chromosome", "start", "stop", "blah", "ugh", "goaway"), 
                                      colClasses = c("character", "numeric", "numeric", rep("NULL", 3)))
}
for (chr in chrom) {
  path_genes <- paste0("/Volumes/Extreme SSD/Williams_Lab/EM_seq/may_2024_resequenced/dmc_by_chr/split_chromosomes/TEs_", chr, ".bed")
  all_data[[chr]]$TEs <- read.table(path_genes, 
                                    col.names = c("chromosome", "start", "stop", "blah", "ugh", "goaway"), 
                                    colClasses = c("character", "numeric", "numeric", rep("NULL", 3)))
}


#ordering back to original order given in samples list
#not sure how to actually do this? maybe combine with the combining things together
#add column to sample that is the sample name 
#rbind each chromosome into one long data frame 

comb_data <- list()

add_sample_name <- function(df, sample_name) {
  name_col <- c(rep(sample_name, nrow(df)))
  df$sample_name <- name_col
  df
}

#not using genes and TEs: 
for (chr in chrom) {
  for (sample in samples_sort) {
    named <- add_sample_name(all_data[[chr]][[sample]], sample)
    comb_data[[chr]] <- rbind(comb_data[[chr]], named)
  }
}

#genes and TEs:
samples_plus_gt <- c("4-1_R0", "4-11_R0", "4-2_R0",  "4-3_R0",  "4-4_R0", "genes", "TEs") 
samples_plus_gt <- c("genes", "TEs")
for (chr in chrom) {
  for (sample in samples_plus_gt) {
    named <- add_sample_name(all_data[[chr]][[sample]], sample)
    comb_data[[chr]] <- rbind(comb_data[[chr]], named)
  }
}

comb_data
#comb_data should be a list of five data frames, each with all samples labeled and concatenated 

#new <- add_sample_name(all_data$Chr1$`4-1_R0`, "bob")
#new

#add sorting step later? if needed idk

#read counts: add summed_reads column which is unmeth + meth
for (chr in chrom) {
  comb_data[[chr]]$summed_reads <- comb_data[[chr]]$name + comb_data[[chr]]$score
}
comb_data

#---------------------------------------------------------------------------------

# Create bins and calculate the sum of rows in each bin
final_bins <- list()
for (chr in chrom) {
  bed_data_bins <- comb_data[[chr]] %>%
    mutate(bin = (floor((stop - 1) / bin_size) + 1) * bin_size) %>% 
    group_by(chromosome, sample_name, bin) %>%
    summarize(sum_rows = n())
  final_bins[[chr]] <- bed_data_bins
}
#final_bins

#set y axis length for all plots to be the highest value in any dataset, rounded to nearest 100
y_axis_len <- c()
for (chr in chrom) {
  y_axis_len <- append(y_axis_len, max(final_bins[[chr]]$sum_rows))
}

final_y <- max(y_axis_len)
final_y <- ceiling(max(final_y)/100)*100
final_y

# Keep only the first data point for each bin
#from what I can tell, this does nothing 
bed_data_first_point <- bed_data_bins %>%
  group_by(chromosome, file, bin) %>%
  summarize(sum_rows = sum(sum_rows)) %>%
  ungroup()
head(bed_data_first_point)



#---------------------------------------------------------------------------------
# Plot the data points and connect them with a line (i think this is the right one to use)

#Set plot width by chr in pixels
plot_width <- list()
plot_width$'Chr1' <- 1500
plot_width$'Chr2' <- 972
plot_width$'Chr3' <- 1159
plot_width$'Chr4' <- 918
plot_width$'Chr5' <- 1332



for (chr in chrom) {
  ggplot(final_bins[[chr]], aes(x = bin, y = sum_rows, color = sample_name)) + theme_classic()+
    geom_line() +
    labs(x = chr, y = "Sum of Rows", color = "Sample") +
    ggtitle(paste("Distribution of DMRs, ", context)) +
    scale_color_manual(values = my_colors) +
    scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M"),expand = expansion(mult = c(0, 0))) + 
    #scale_y_log10() 
    scale_y_continuous(
      name = "number of features",
      #trans = 'log2'
      breaks = seq(0, final_y, by = 50),
      labels = seq(0, final_y, by = 50),
      limits = c(0, final_y), # Set the y-axis limits
      expand = expansion(mult = c(0, 0)) ) 
  ggsave(path = output_dir, filename = paste(title, chr, context, " 350kb genes TEs.svg"), 
         width = plot_width[[chr]], height = 750, unit = "px")
}


#---------------------------------------------------------------------------------
# Plot the data points and connect them with a line
#plotting summed_rows 

#old? don't need this? 


y_axis_len <- c()
for (chr in chrom) {
  y_axis_len <- append(y_axis_len, max(comb_data[[chr]]$summed_reads))
}

final_y <- max(y_axis_len)
final_y <- ceiling(max(final_y)/1000)*1000

for (chr in chrom) {
  ggplot(comb_data[[chr]], aes(x = start, y = summed_reads, color = sample_name)) + 
    geom_line(alpha = .5) + 
    scale_color_manual(values = my_colors) +
    theme_classic() + 
    labs(x = chr) + 
    scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M"),expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(
      name = "total read counts",
      breaks = seq(0, final_y, by = 100),
      labels = seq(0, final_y, by = 100),
      limits = c(0, final_y), # Set the y-axis limits
      expand = expansion(mult = c(0, 0)) )
  ggsave(path = output_dir, filename = paste("R0 dmr distribution", title, chr, context, ".png"), 
         width = plot_width[[chr]], height = 750, unit = "px")
}

p <- ggplot(comb_data$Chr3, aes(x = start, y = summed_reads, color = sample_name)) + 
  geom_line(alpha = .5) + 
  scale_color_manual(values = my_colors) +
  theme_classic() + 
  labs(x = "chr 3") + 
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M"),expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(
    name = "total read counts",
    breaks = seq(0, 2000, by = 1000),
    labels = seq(0, 2000, by = 1000),
    limits = c(0, 2000), # Set the y-axis limits
    expand = expansion(mult = c(0, 0)) )

p

