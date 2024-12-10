#Script to build a complex heat map of a list of genomic features across chromosomes
#input: a .txt file of DMRs from DMRcaller
#documentation:
#https://jokergoo.github.io/ComplexHeatmap-reference/book/genome-level-heatmap.html

#---------------------------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("EnrichedHeatmap")
BiocManager::install("rtracklayer")
install.packages('circlize')
BiocManager::install("GenomicRanges", force = TRUE)

library(circlize)
library(GenomicRanges)
library(EnrichedHeatmap)
library(rtracklayer)
library(readr)
library(dplyr)

#---------------------------------------------------------------------------------
###User set variables
path = "/Users/nksmoot/Desktop/Grad/r studio things"
subgroup_col = c("4-1_R0" = "#F5B085",
                 "4-2_R0" = "#51D0C1")
context = c('CG')
samples = c('4-1_R0', '4-2_R0')

#---------------------------------------------------------------------------------
#create GRanges object of chromosomes
chromes = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrMt", "ChrCp")
start = c(0, 0, 0, 0, 0, 0, 0)
end = c(30427671, 26975502, 23459830, 19698289, 18585056, 366924, 154478)
chr_df = data.frame(chromes, start, end)
chr_gr = GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))

#Divide chromosomes into windows of specified size 
chr_window = makeWindows(chr_gr, w=250000)

#Import files
X4_1_R0_v_WT_DMRsBinsCG <- read_table("4-1_R0_v_WT_DMRsBinsCG.txt")
head(X4_1_R0_v_WT_DMRsBinsCG)
new <- make_clean(X4_1_R0_v_WT_DMRsBinsCG)
head(new)

#---------------------------------------------------------------------------------
#Import files
#file names need to be in the format "samplename_v_WT_DMRsBinsCG.txt"
import_data <- function(path, tissue, context) {
  bed_paths <- list.files(path, recursive = T,
                          pattern = paste0("^(", paste(tissue, collapse = "|"),
                                           ")", "_v_WT_DMRsBins", context, ".bed"),
                          full.names = T)
  data <- lapply(bed_paths, read.table, header = F, 
                 col.names = c('seqnames', 'start', 'end', 'width', 'strand', 'sumReadsM1', 'sumReadsN1', 'proportion1', 'sumReadsM2', 'sumReadsN2', 'proportion2', 'cytosinesCount', 'context', 'direction', 'pValue', 'regionType'))
  names(data) <- tissue
  return(data)
}

all_data <- lapply(context, import_data, path = path, tissue = samples)
names(all_data) <- "samples"
head(all_data)


#remove columns that are not needed, * will make an error
clean_file <- select(X4_1_R0_v_WT_DMRsBinsCG, c('seqnames', 'start', 'end', 'proportion1', 'proportion2'))

#add column which is the percent change from WT
clean_file$delta <- (clean_file$proportion2 - clean_file$proportion1)
gr_obj = GRanges(clean_file)
vector <- clean_file$delta

#Need to iterate over a list of samples, and for each one produce a clean version
#and append the delta column to a matrix
make_clean <- function(df) {
  clean <- select(df, c('seqnames', 'start', 'end', 'proportion1', 'proportion2'))
  clean$delta <- (clean$proportion2 - clean$proportion1)
  gr_obj = GRanges(clean)
  gr_obj
}

make_vectors <- function(df) {
  df$delta
}

#works for a single file
X4_1_R0_v_WT_DMRsBinsCG <- read_table("4-1_R0_v_WT_DMRsBinsCG.txt")
head(X4_1_R0_v_WT_DMRsBinsCG)
new <- make_clean(X4_1_R0_v_WT_DMRsBinsCG)
head(new)

#These values are correct 
clean_CG <- lapply(all_data$samples, make_clean)
clean_CG
vectors_CG <- lapply(clean_CG, make_vectors)
vectors_CG

###Oh... problem with combining things that have different numbers of DMRs. Need to call average in window on each sample, 
#and then try to combine them into one matrix 

#this function calculates the average change in methylation in each window that we previously specified 
average_in_window = function(window, gr, v, method = "w0", empty_v = NA) {
  
  if(missing(v)) v = rep(1, length(gr))
  if(is.null(v)) v = rep(1, length(gr))
  if(is.atomic(v) && is.vector(v)) v = cbind(v)
  
  v = as.matrix(v)
  if(is.character(v) && ncol(v) > 1) {
    stop("`v` can only be a character vector.")
  }
  
  if(length(empty_v) == 1) {
    empty_v = rep(empty_v, ncol(v))
  }
  
  u = matrix(rep(empty_v, each = length(window)), nrow = length(window), ncol = ncol(v))
  
  mtch = as.matrix(findOverlaps(window, gr))
  intersect = pintersect(window[mtch[,1]], gr[mtch[,2]])
  w = width(intersect)
  v = v[mtch[,2], , drop = FALSE]
  n = nrow(v)
  
  ind_list = split(seq_len(n), mtch[, 1])
  window_index = as.numeric(names(ind_list))
  window_w = width(window)
  
  if(is.character(v)) {
    for(i in seq_along(ind_list)) {
      ind = ind_list[[i]]
      if(is.function(method)) {
        u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
      } else {
        tb = tapply(w[ind], v[ind], sum)
        u[window_index[i], ] = names(tb[which.max(tb)])
      }
    }
  } else {
    if(method == "w0") {
      gr2 = reduce(gr, min.gapwidth = 0)
      mtch2 = as.matrix(findOverlaps(window, gr2))
      intersect2 = pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])
      
      width_intersect = tapply(width(intersect2), mtch2[, 1], sum)
      ind = unique(mtch2[, 1])
      width_setdiff = width(window[ind]) - width_intersect
      
      w2 = width(window[ind])
      
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        x = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
        u[window_index[i], ] = (x*width_intersect[i] + empty_v*width_setdiff[i])/w2[i]
      }
      
    } else if(method == "absolute") {
      for(i in seq_along(ind_list)) {
        u[window_index[i], ] = colMeans(v[ind_list[[i]], , drop = FALSE])
      }
      
    } else if(method == "weighted") {
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        u[window_index[i], ] = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
      }
    } else {
      if(is.function(method)) {
        for(i in seq_along(ind_list)) {
          ind = ind_list[[i]]
          u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
        }
      } else {
        stop("wrong method.")
      }
    }
  }
  
  return(u)
}

#call average_in_window
output <- average_in_window(chr_window, gr_obj, vector, method = "weighted", empty_v = NA)
type(output)

chr_window
total_data <- matrix(NA, nrow = length(chr_window), ncol = 1)
total_data
df_data <- data.frame(NA, nrow = length(chr_window), ncol = 1)
df_data

for (i in seq_along(clean_CG)) {
  gr_obj <- clean_CG[[i]]
  delta_vector <- vectors_CG[[i]]
  
  # Calculate average methylation level in each window
  avg_methylation <- average_in_window(chr_window, gr_obj, delta_vector, method = "weighted", empty_v = NA)
  
  # Add the calculated values to total_data
  total_data <- cbind(total_data, avg_methylation)
}

# Remove the first column of NA values from total_data
total_data <- total_data[, -1]

total_data
colnames(total_data) <- samples

df <- as.data.frame(total_data)
names(df) <- samples
head(df)


image(1:ncol(total_data), 1:nrow(total_data), t(total_data), col = terrain.colors(60), axes = FALSE)
axis(1, 1:ncol(total_data), colnames(total_data))
axis(2, 1:nrow(total_data), rownames(total_data))
for (x in 1:ncol(total_data))
  for (y in 1:nrow(total_data))
    text(x, y, total_data[y,x])


library(ggplot2)



new_blue = "#4eaadb"
new_red = "#db4d4c"
neutral = "#f5fab7"

chr = as.vector(seqnames(chr_window))
chr_level = paste0("Chr", 1:5)
chr = factor(chr, levels = chr_level)
subgroup = samples
subgroup
sub_level = factor(subgroup, subgroup)


library(ComplexHeatmap)
ht_opt$TITLE_PADDING = unit(c(4, 4), "points")
ht_list = Heatmap(total_data, name = "% change in meth", col = colorRamp2(c(-1, 0, 1), c(new_blue, neutral, new_red)),
                  row_split = chr, cluster_rows = FALSE, show_column_dend = FALSE,
                  column_split = subgroup, cluster_column_slices = FALSE,
                  column_title = paste(context, "DMRs"),
                  top_annotation = HeatmapAnnotation(subgroup = subgroup, annotation_name_side = "left", 
                                                     col = list(subgroup = subgroup_col)),
                  row_title_rot = 0, row_title_gp = gpar(fontsize = 10), border = TRUE,
                  row_gap = unit(0, "points")) 
 

draw(ht_list, merge_legend = TRUE)









