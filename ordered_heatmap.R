#!/usr/bin/env Rscript

#SBATCH --job-name=chromosome_heatmap

#SBATCH --account=fc_williamslab

#SBATCH --partition=savio2_bigmem

#SBATCH --time=01:00:00

#Script to build a complex heat map of a list of genomic features across chromosomes
#input: a .txt file of DMRs from DMRcaller
#documentation:
#https://jokergoo.github.io/ComplexHeatmap-reference/book/genome-level-heatmap.html

#---------------------------------------------------------------------------------
#Ensure you have the following packages installed: 
#through BiocManager: GenomicRanges, EnrichedHeatmap, rtracklayer, circlize, 
#just install.packages: circlize

library(circlize)
library(GenomicRanges)
library(EnrichedHeatmap)
library(rtracklayer)
library(readr)
library(dplyr)
library(ComplexHeatmap)

#---------------------------------------------------------------------------------
###User set variables
####if you're changing the file type naming schema, go down to the import files chunk and make appropriate changes there
path = "/Users/nksmoot/Desktop/Grad/r studio things"
#path_to_files = "/Volumes/Extreme\ SSD/Williams_Lab/EM_seq/may_2024_resequenced/dmrs/"
path_to_files = "/Volumes/Extreme\ SSD/Williams_Lab/EM_seq/may_2024_resequenced/bed/"

#filter values that are NA for all samples
filter_NA = TRUE

#colors for the heat map
new_blue = "#36AFC1"
new_red = "#CA572B"
neutral = "#FFFAE0"

low ="#FBF3CE"
mid="#C85133"
high="#8D104C"

context = c('_all_CpG')
samples = c('DRDD', 'gen4_drdd', '4-1_R0', '4-2_R0', '4-3_R0', '4-4_R0', '4-11_R0')
#samples = c('gen2_drdd', '4-2_R0', '4-2_R1', '4-14_R0', '4-15_R0')
#samples = c('gen2_drdd', '4-2_R0', '4-2_R1', '4-4_R0', '4-4_R1')

window_size = 150
context
title = paste("R0s_", context, "CG DMRs in preR")

#---------------------------------------------------------------------------------
#Import files
#file names need to be in the format "samplename_v_WT_DMRsBinsContext.bed"
#be careful about this function, ken wrote it and it's buggy 
#main problem: this function will import ur samples ALPHABETICALLY 
#and then rename them the order of ur sample list, so the wrong name will be associated with the imported data
#to solve this, we will order sample names, apply ordered names to the ordered files, then reorder to the order initially set in samples
samples_sort <- sort(samples)
samples_sort

import_data <- function(path, tissue, context) {
  bed_paths <- list.files(path, recursive = F,
                          pattern = paste0("^(", paste(tissue, collapse = "|"),
                                           ")", context, ".bed"),
                          full.names = T)

  #data <- lapply(bed_paths, read.table, header = F, 
  #              col.names = c('seqnames', 'start', 'end', 'width', 'strand', 'sumReadsM1', 'sumReadsN1', 'proportion1', 'sumReadsM2', 'sumReadsN2', 'proportion2', 'cytosinesCount', 'context', 'direction', 'pValue', 'regionType'))
  data <- lapply(bed_paths, read.table, header = F, 
                col.names = c('seqnames', 'start', 'end', 'unmeth', 'meth', 'prop_meth'))
  names(data) <- tissue
  #names(data) <- gsub("(.*)_v_WT.*", "\\1", bed_paths)
  return(data)
}

all_data <- lapply(context, import_data, path = path_to_files, tissue = samples_sort)
names(all_data) <- "samples"

#reordering the imported data to the order provided in samples
reordered_data <- list()
reordered_data$samples <- lapply(samples, function(x) {
  return(all_data$samples[[x]])})
names(reordered_data$samples) <- samples
all_data <- reordered_data


#---------------------------------------------------------------------------------
#create GRanges object of chromosomes
#### This section is setting the windows which will be used to generate the heat map

#Use this chunk if you want to have the entire genome, divided into even windows of a size defined earlier in the user set variable section
chromes = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
start = c(0, 0, 0, 0, 0)
end = c(30427671, 19700000, 23465000, 18585000, 26975800)
chr_df = data.frame(chromes, start, end)
chr_gr = GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))

#Divide chromosomes into windows of specified size 
chr_window = makeWindows(chr_gr, w=window_size)


#### Alternative option: provide a set of windows (ie, your differentially methylated regions) to be used: 

#I've just kept various chunks of code that I've run before; modify or add your own here. 
#The important thing is that the windows you want to use are saved as GRanges object chr_window

#import superset of DMRs to use as windows
superset_dmrs <- read.table("/Volumes/Extreme\ SSD/Williams_Lab/EM_seq/may_2024_resequenced/bed/named_merged_R0_CG_DMRs.bed", header = FALSE, fill = TRUE)
colnames(superset_dmrs) = c('seqnames', 'start', 'end','sample','seqname2','start2','end2', 'width', 'strand', 'sumReadsM1', 'sumReadsN1', 
                            'proportion1', 'sumReadsM2', 'sumReadsN2', 'proportion2', 'cytosinesCount', 
                            'context', 'direction', 'pValue', 'regionType')

#updated superset 10/30/24
superset_dmrs <- read.table(
  "/Volumes/Extreme\ SSD/Williams_Lab/EM_seq/may_2024_resequenced/unique_dmrs/gen4_drdd_v_WT_DMRsBinsCG_without_Ws.bed", 
  header = FALSE, fill = TRUE)
colnames(superset_dmrs) = c('seqnames', 'start', 'end')
chr_window <- GRanges(superset_dmrs, ignore.strand=TRUE)

#superset file loses information about strand, which granges requires, so i'm just adding it back idk? 
#update: there's an option to tell granges to ignore strand information (I think you need to use makeGrangesFromDF function?)
superset_dmrs$test = c(rep("*", nrow(superset_dmrs)))
test_df <- data.frame(superset_dmrs$seqnames, superset_dmrs$start, superset_dmrs$end, superset_dmrs$sample, superset_dmrs$test)
colnames(test_df) <- c('seqnames', 'start', 'end', 'sample', 'strand')
chr_window <- GRanges(test_df)
nrow(test_df)

#---------------------------------------------------------------------------------
#Iterate over a list of samples, and for each one produce a clean version without columns with asterixes 
#and append the delta column, which is the difference between proportion 2 (methylation in mutant) and 1 (methylation in WT)
#Can edit this delta column to be whatever you want to be heat mapped 

#delta is a column that measures whatever you want to be heat mapped (ie, methylation level, as is written here)

make_clean <- function(df) {
  clean <- select(df, c('seqnames', 'start', 'end', 'prop_meth'))
  clean$delta <- (clean$prop_meth)
  gr_obj = GRanges(clean)
  gr_obj
}

make_vectors <- function(df) {
  df$delta
}

clean_CG <- lapply(all_data$samples, make_clean)
vectors_CG <- lapply(clean_CG, make_vectors)

clean_CG

#---------------------------------------------------------------------------------
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

#---------------------------------------------------------------------------------
total_data <- matrix(NA, nrow = length(chr_window), ncol = 1)

#iterate through each sample, calculate average methylation change in each window, and add that to the total_data matrix
for (i in seq_along(clean_CG)) {
  gr_obj <- clean_CG[[i]]
  delta_vector <- vectors_CG[[i]]
  
  # Calculate average methylation change in each window
  avg_methylation <- average_in_window(chr_window, gr_obj, delta_vector, method = "weighted", empty_v = NA)
  
  # Add the calculated values to total_data matrix
  total_data <- cbind(total_data, avg_methylation)
}

# Remove the first column of NA values from total_data (created when total_data was initialized)
total_data <- total_data[, -1]
colnames(total_data) <- samples

#---------------------------------------------------------------------------------
#setting levels to control ordering of things in the heatmap
#I don't totally remember but I think you need to run the subgroups lines
#this will order by chromosome if clustering is off in heatmap function 
chr = as.vector(seqnames(chr_window))
chr_level = paste0("Chr", 1:5)
chr = factor(chr, levels = chr_level)
subgroup <- samples
subgroup = factor(subgroup, levels = subgroup)
subgroup
head(total_data)

#removing rows which are NA for all samples
if (filter_NA == TRUE) {
  print("removing NAs")
  new_data <- total_data[rowSums(is.na(total_data)) != ncol(total_data), ]
  head(new_data)
  total_data <- new_data
}


#set NA values to zero (if not, they're not ordered correctly)
total_data[is.na(total_data)] <- 0


length(total_data)
df <- data.frame(total_data)
head(df)
df$number <- c(1:nrow(df))
df

#order the matrix by level of methylation change in drdd gen2 
df_ordered <- df[order(-df$gen2_drdd),]
head(df_ordered)
nrow(df_ordered)
matrix_ordered <- data.matrix(df_ordered, rownames.force = NA)
length(matrix_ordered)
matrix_ordered <- matrix_ordered[,!colnames(matrix_ordered) %in% 'number']
colnames(matrix_ordered) <- samples



head(matrix_ordered)
rownames(matrix_ordered) <- NULL
length(matrix_ordered)


#chr = factor(chr, levels = ordered_gen2)

#---------------------------------------------------------------------------------
#constructing the heatmap
svglite(file = paste0("/Users/nksmoot/Desktop/Grad/r studio things/Manuscript/Fig\ 2\ methylome/", title," preR_CG_dmrs.svg"))



ht_opt$TITLE_PADDING = unit(c(4, 4), "points")
ht_list = Heatmap(total_data, name = "% meth", col = colorRamp2(c(0, 50, 100), c(low, mid, high)), na_col = "#edeae9",
                  #row_split = chr, 
                  cluster_rows = TRUE, show_column_dend = TRUE,
                  column_split = subgroup, 
                  cluster_column_slices = FALSE,
                  column_title = title,
                  #top_annotation = HeatmapAnnotation(subgroup = subgroup, annotation_name_side = "left", 
                  #                                   col = list(subgroup = subgroup_col)),
                  #row_title_rot = 0, 
                  #row_title_gp = gpar(fontsize = 10), 
                  border = TRUE,
                  row_gap = unit(0, "points")) 


draw(ht_list, merge_legend = TRUE)

dev.off()






