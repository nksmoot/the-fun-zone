#parse the Xu time lapse data and spit out two nice graphs 
gene_name <- "PLT3"
gene_ID <- "AT5G10510"

library(reshape2)



df <- read_excel("~/Desktop/Grad/RNA seq data/xu_timelapse/xu_timelapse.xlsx")


head(df)

gene_of_interest <- df[df$gene == gene_ID, ]
gene_of_interest


early_timepoints <- gene_of_interest[, c(3:18)]
late_timepoints <- gene_of_interest[, c(21:36)]
library(data.table)

colnames(early_timepoints) <- c(0, 0, .17, .17, .5, .5, 1, 1, 2, 2, 4, 4, 8, 8, 12, 12)
colnames(late_timepoints) <- c(0, 0, .25, .25, .5, .5, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5)

melted_early <- melt(as.data.table(early_timepoints))
head(melted_early)
melted_early$variable <- as.numeric(levels(melted_early$variable))[melted_early$variable]
head(melted_early)


melted_late <- melt(as.data.table(late_timepoints))
melted_late$variable <- as.numeric(levels(melted_late$variable))[melted_late$variable]

ggplot(melted_early, aes(x = variable, y = value)) + 
  geom_point() + 
  theme_classic() + 
  stat_summary(fun = "mean", geom = "line", color="#AC1E38") + 
  ggtitle(paste("Early wound response,", gene_name)) + 
  xlab("time in hours") + ylab("TPM") + 
  scale_y_continuous(expand=c(0,0))


ggplot(melted_late, aes(x = variable, y = value)) + 
  geom_point() + 
  theme_classic() + 
  stat_summary(fun = "mean", geom = "line", color="#8D104C") + 
  ggtitle(paste("Late wound response,", gene_name)) + 
  xlab("time in days") + ylab("TPM") + 
  scale_y_continuous(expand=c(0,0))

