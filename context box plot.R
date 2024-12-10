#trying to split DMRs into "CG only" and "multi context" 
#for each CG DMR, does it also show >5% gain in other contexts? 
#output: fig 2D from plant cell. should have WT CG, CHG, CHH and then sample CG, CHG, CHH
#might have been a good idea to work with smaller versions of these files lol
#---------------------------------------------------------------------------------
library(GenomicRanges)
library(reshape2)
library(ggplot2)


samples = c("DRDD", "gen4_drdd", "4-4_R0", "4-4_R1")
sorted_samples <- sort(samples)

title <- paste("all contexts CHG at 4-4 CHG uniaue dmrs")

path_to_files = "/Volumes/Extreme\ SSD/Williams_Lab/EM_seq/may_2024_resequenced/bed/"

#---------------------------------------------------------------------------------
#Import CG, CHG, and CHH info for a sample

all_data <- list()

contexts = c("_all_CpG", "_all_CHG", "_all_CHH")

import_data <- function(path, tissue, context) {
  bed_paths <- list.files(path, recursive = F,
                          pattern = paste0("^(", paste(tissue, collapse = "|"),
                                           ")", context, ".bed"),
                          full.names = T)
  data <- lapply(bed_paths, read.table, header = F, 
                 col.names = c("Chr", "start", "end", "unmeth", "meth", "prop_meth"))
  names(data) <- tissue
  return(data)
}

all_data <- lapply(contexts, import_data, path = path_to_files, tissue = sorted_samples)
head(all_data)
names(all_data) <- contexts
#data get organized back to the samples list order later, in average window calculations

#---------------------------------------------------------------------------------
######superset windows
superset_dmrs <- read.table("/Volumes/Extreme\ SSD/Williams_Lab/EM_seq/may_2024_resequenced/bed/named_merged_R0_CG_DMRs.bed", header = FALSE, fill = TRUE)
colnames(superset_dmrs) = c('seqnames', 'start', 'end','sample','seqname2','start2','end2', 'width', 'strand', 'sumReadsM1', 'sumReadsN1', 
                            'proportion1', 'sumReadsM2', 'sumReadsN2', 'proportion2', 'cytosinesCount', 
                            'context', 'direction', 'pValue', 'regionType')
#superset file loses information about strand, which granges requires, so i'm just adding it back idk? 
superset_dmrs$test = c(rep("*", nrow(superset_dmrs)))
test_df <- data.frame(superset_dmrs$seqnames, superset_dmrs$start, superset_dmrs$end, superset_dmrs$sample, superset_dmrs$test)
colnames(test_df) <- c('seqnames', 'start', 'end', 'sample', 'strand')
gr_windows <- GRanges(test_df)


##### single sample DMRs
dmrs_single <- read.table("/Volumes/Extreme SSD/Williams_Lab/EM_seq/may_2024_resequenced/unique_dmrs/CHG/changes_from_pre_4-4_R0_CHG.bed")
colnames(dmrs_single) <- c('seqnames', 'start', 'end', 'width', 'strand', 'sumReadsM1', 'sumReadsN1', 
                        'proportion1', 'sumReadsM2', 'sumReadsN2', 'proportion2', 'cytosinesCount', 
                        'context', 'direction', 'pValue', 'regionType')
gr_windows <- GRanges(dmrs_single)

##### single sample DMRs
dmrs_single <- read.table("/Volumes/Extreme SSD/Williams_Lab/EM_seq/may_2024_resequenced/unique_dmrs/CHG/unique_4-4_R0_v_R0s_CHG.bed")
colnames(dmrs_single) <- c('seqnames', 'start', 'end', 'width', 'strand', 'sumReadsM1', 'sumReadsN1', 
                           'proportion1', 'sumReadsM2', 'sumReadsN2', 'proportion2', 'cytosinesCount', 
                           'context', 'direction', 'pValue', 'regionType')
gr_windows <- GRanges(dmrs_single)


#---------------------------------------------------------------------------------
#this function calculates the average change in methylation in each window that we previously specified 
average_in_window = function(window, gr, v, method = "absolute", empty_v = NA) {
  
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
#calculate average methylation across windows: 

#initialize data frame
df = data.frame()
df = data.frame(nrow = length(gr_windows))

#for each sample, convert to a GRanges object, extract the level of methylation at each cytosine, 
#and calculate the average level of cytosine methylation across your specified windows
#append to the dataframe 
#possibly switch this to lapply, apparently for loops are slow in r?
for (sample in samples) {
  for (context in contexts) {
    print("now calculating")
    print(context)
    print(sample)
    gr <- sort(GRanges(all_data[[context]][[sample]]))
    level_meth <- gr$prop_meth
    avg_in_win <- average_in_window(gr_windows, gr, level_meth, method = "absolute", empty_v = NA)
    df <- data.frame(df, avg_in_win)
  } }
  
#I don't really know why this works.. have to initialize the df with some dimension?? but it's not the correct dimension it's just some random column. oh well. 
df$nrow <- NULL

#name them CpG_DRDD CpG_4-4
names <- c()
for (sample in samples) {
  for (context in contexts) {
    names <- c(names, paste0(context, "_", sample) )
  }
  
}
colnames(df) <- names
head(df)

#save this file as a thing I don't use the name of
saved_unique_4_4_CG <- df
saved_gains_CHG_4_4 <- df
saved_unique_CHG_4_4 <- df
#rddm_targets: 135 out of 475 unique regions seem to be drdd targets. >5% gain in chg and 3% in chh


#---------------------------------------------------------------------------------
#reshape data
melted_df <- melt(df)

#---------------------------------------------------------------------------------
#make plot
p <- ggplot(melted_df, aes(x=variable, y=value)) + 
  geom_boxplot(outlier.alpha = 0, aes(fill = variable)) +
  #geom_violin(aes(fill = variable)) + geom_boxplot(outlier.alpha=0, width = .1) + 
   
 # geom_jitter(size = .001, alpha = .1, aes(fill = variable)) + 
  xlab("samples") + ylab("methylation in each region") + 
  ggtitle("Methylation level at 4-4 unique CHG DMRs") + 
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +  
  scale_fill_manual(values = c("#F5B085", "#EC8469", "#C73347","#BDDD4B", "#76BB32", "#3F9B16",  
                               "#6BDCEB",  "#4DAADA", "#3178CB", "#6BDCEB",  "#4DAADA", "#3178CB"))
p + theme_classic() + theme(legend.position = "none")

title <- "all contexts CHG at unique 4-4 CG dmrs"

#dim = length(samples)*2 + 1.5

ggsave(paste0("/Users/nksmoot/Desktop/Grad/r studio things/Manuscript/Fig\ 2\ methylome/", title, "multi context box_with R1_CORRECT_box.svg"))




