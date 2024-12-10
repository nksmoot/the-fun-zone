#distribution, type, and size of 4g4.vcf
library(readxl)
df <- read_excel("/Volumes/Extreme SSD/Williams_Lab/pac_bio/vcf_analysis/4g4_graph.xlsx")


head(df)
test <- df[df$`#CHROM`=="Chr2",]


chrom <- paste0("Chr", 1:5)
chrom

all_data <- list()

for (chr in chrom) {
  all_data[[chr]] <- df[df$`#CHROM`==chr,]
}
all_data


ggplot(all_data$Chr1, aes(x=POS, y=type, color=type, size=log)) + 
  geom_point(alpha=.1) + 
  theme_classic() + 
  ylim(c(-5, 5)) + 
  scale_color_gradient(low="#798234", high="#D4677F") + 
  ggtitle("chr 1")
ggsave(paste0("structural_variation/", "Chr1", " variant calls.svg"))

#Set plot width by chr in pixels
plot_width <- list()
plot_width$'Chr1' <- 1500
plot_width$'Chr2' <- 972
plot_width$'Chr3' <- 1159
plot_width$'Chr4' <- 918
plot_width$'Chr5' <- 1332


for (chr in chrom) {
  ggplot(all_data[[chr]], aes(x=POS, y=type, color=type, size=log)) + 
    geom_point(alpha=.1) + 
    theme_classic() + 
    ylim(c(-5, 5)) + 
    scale_color_gradient(low="#798234", high="#D4677F") +
    ggtitle(chr)
  ggsave(paste0("structural_variation/", chr, " variant calls.svg"), 
         width = plot_width[[chr]], height = 750, unit = "px")
}
