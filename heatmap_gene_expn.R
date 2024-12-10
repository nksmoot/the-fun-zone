coord = read.table("/Users/nksmoot/Desktop/Grad/RNA seq data/combined_d4_dR/sp_updated_combined_allDEG_counts.txt", 
                       header=F, sep="\t")
coord <- coord[!duplicated(coord),]
rownames(coord) <- coord$V1
coord$V1 <- NULL
colnames(coord) <- c("DW_1", "DW_2", "DW_3", "DW_4", "DW_5", "DW_6", "DW_7", 
                     "d4_1", "d4_2", "d4_3", "d4_4", "d4_5", "d4_6", "d4_7", 
                     "dR_1", "dR_2", "dR_3", "dR_4", "dR_5", "dR_6", "dR_7")

b1 <- coord[c("DW_1", "DW_2", "DW_3", "d4_1", "d4_2", "d4_3", "dR_1", "dR_2", "dR_3")]
b2 <- coord[c(4:7,11:14,18:21)]
coord_scale <- t(scale(t(coord), center=T, scale=T))
coord_scale <- as.data.frame(coord_scale)
coord_scale <- na.omit(coord_scale)
cor.exp <- as.data.frame(cor(coord_scale))
cor.exp$samples <- c("DW_1", "DW_2", "DW_3", "DW_4", "DW_5", "DW_6", "DW_7", 
                     "d4_1", "d4_2", "d4_3", "d4_4", "d4_5", "d4_6", "d4_7", 
                     "dR_1", "dR_2", "dR_3", "dR_4", "dR_5", "dR_6", "dR_7")

cor.exp$samples <- factor(cor.exp$samples, levels = rev(c("DW_1", "DW_2", "DW_3", "DW_4", "DW_5", "DW_6", "DW_7", 
                                                      "d4_1", "d4_2", "d4_3", "d4_4", "d4_5", "d4_6", "d4_7", 
                                                      "dR_1", "dR_2", "dR_3", "dR_4", "dR_5", "dR_6", "dR_7")))


svg(filename="superset_allDEGs_batch1.svg")
pheatmap(cor.exp, color=colorRampPalette(c("#4DAADA", "white", "#DB4D4B"))(100),
         breaks=seq(-1,1, length.out=101), cluster_cols = F, cluster_rows = F)

dev.off()

melted <- melt(cor.exp, id=c("samples"))

ggplot(melted, aes(x=samples, y=variable, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#4DAADA", 
                       mid = "white", 
                       high = "#DB4D4B") + 
  theme_classic()
ggsave("whole set superset all degs.svg")



pca <- read.table("/Users/nksmoot/Desktop/Grad/RNA seq data/combined_d4_dR/test_WT_v_d4_pca_coord.txt")
pca$batch <- c(rep("WT_b1", 3), rep("WT_b2", 4), rep("d4_b1", 3), rep("d4_b2", 4), rep("dR_b1", 3), rep("dR_b2", 4))


ggplot(pca, aes(x=PC1, y=PC2, color=batch)) + 
  geom_point(size=3, alpha=.3) + 
  theme_classic() + 
  scale_color_manual(values=c("#8CC83F", "#8CC83F", "#4DABDA", "#4DABDA", "#F5B085", "#F5B085"))

ggplot(pca, aes(x=genotype, y=PC1, color=batch)) + 
  geom_boxplot() +
  geom_point(position = position_dodge(width=.75)) + 
  theme_classic() + 
  scale_color_manual(values=c("#D6F568", "#8CC83F", "#6CDDEC", "#4DABDA", "#FAE0BF", "#F5B085"))


ggplot(pca, aes(x=genotype, y=PC1, color=genotype)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_point() + 
  theme_classic() + 
  scale_color_manual(values=c("#8CC83F",  "#4DABDA", "#F5B085"))
ggsave("PC2 7 samples.svg")



all_pcs$sdev^2/sum(all_pcs$sdev^2)
[1] 3.298940e-01 1.767068e-01 1.389459e-01 6.532650e-02 4.845396e-02 3.794012e-02 3.019000e-02 2.599116e-02 2.132233e-02 1.806796e-02 1.610951e-02 1.364746e-02
[13] 1.185792e-02 1.097569e-02 1.056056e-02 1.051745e-02 9.480193e-03 8.772247e-03 8.082145e-03 7.158139e-03 1.767645e-30