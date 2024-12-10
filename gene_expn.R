##somehow do the gene expression analysis? 

library(readr)
library(ggplot2)

WT_vs_dR_21dpg_2kb_expn <- read_delim("WT vs dR 21dpg_2kb_expn.txt", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)
head(WT_vs_dR_21dpg_2kb_expn)

df <- WT_vs_dR_21dpg_2kb_expn

p <- ggplot(df, aes(x = delta_meth, y = log2foldchange, fill = delta_meth)) + 
  geom_violin() +
  geom_jitter(size = .1, alpha = .5) + 
  geom_boxplot(width = .15, fill = "white", outlier.size = 0) + 
  xlab("change in methylation \nof nearby regions") + ylab("expression change") + 
  scale_fill_manual(values = c(new_red, new_blue)) + 
  ggtitle("Change in expression \nof genes associated \nwith DMRs")
p + theme_classic() + theme(legend.position = "none")
ggsave("violin_meth_expn.png", width = 2.5, height = 4, dpi = 400, units = "in")


p <- ggplot(df, aes(x = delta_prop_meth, y = log2foldchange)) + 
  geom_point(alpha = .5, color = "#A5A59B") + 
  xlab("change in methylation from WT") + 
  ylab("change in expression from WT") + 
  ggtitle("Correlation between DMRs within 2kb of DEGs")
p + theme_classic()


perc_change <- c(41, 30)
delta <- c("down", "up")
df_change <- data.frame(delta, perc_change)
df_change

p <- ggplot(df_change, aes(x = delta, y = perc_change, fill = delta)) + 
  geom_bar(stat = "identity", color = "black") + 
  scale_fill_manual(values = c(new_blue, new_red)) + 
  ggtitle("Percent of DEGs \nwhich are within \n2kb of a DMR") + 
  xlab("change in expression") + ylab("percent of DEGs within 2kb") + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 50))

p + theme_classic() + theme(legend.position = "none")
ggsave("dR_correlation.svg", width = 2.5, height = 4, units = "in", dpi = 400)  
  
  