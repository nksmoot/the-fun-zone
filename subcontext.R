#make box plot of CHH subcontexts
library(ggplot2)
path="/Volumes/Extreme SSD/Williams_Lab/EM_seq/may_2024_resequenced/bed/centromere_chh_sub/"
input="4-2_"
gen="R1"
sample=paste0(input, gen, "_uhCHHdmrs_")
location="not_centromere"
title=paste0(sample,location)
title
title <- "4-2_R0_uhCHHdmrs_not_centromere.bed"

#need to load: 4-2_R1_uhCHHdmrs_not_centromere.bed
regen_test <- read.table(paste0(path,title,".bed"))
regen <- regen_test

regen <- read.table(paste0(path, "4-2_R0_uhCHHdmrs_not_centromere.bed"))


#need to load: gen2_in_4-2_uhCHHdmrs_not_centromere.bed
gen2_CHH <- read.table(paste0(path,"gen2_in_", input, "uhCHHdmrs_", location,".bed"))
gen2 <- gen2_CHH

gen2 <- read.table(paste0(path, "gen2_in_4-2_uhCHHdmrs_not_centromere.bed"))

names <- c("Chr", "Start", "End", "unmeth", "meth", "prop", 
           "chr_again", "start_again", "end_again", "context", "subcontext")

colnames(regen) <- names
colnames(gen2) <- names

small_names <- c("prop_meth", "subcontext", "sample")
regen_small <- data.frame(regen$prop, regen$subcontext)
regen_small$sample <- "regen"
colnames(regen_small) <- small_names

gen2_small <- data.frame(gen2$prop, gen2$subcontext)
gen2_small$sample <- "gen2"
colnames(gen2_small) <- small_names

both_samples <- rbind(regen_small, gen2_small)

#bias is towards CAA and CTA, order those towards the front
both_samples$subcontext <- factor(both_samples$subcontext, 
                                  levels=c("CAA", "CTA", 
                                           "CAC", "CAT", 
                                           "CCA", "CCC", "CCT", 
                                           "CGA", "CGT", 
                                           "CTC", "CTT"))

table(both_samples$subcontext)
both <- na.omit(both_samples)

#both_samples$subcontext <- factor(both_samples$subcontext, 
#                                  levels=c("CAG", CTG", "CCG", 
#                                           "CAA", "CAC", "CAT", 
#                                           "CCA", "CCC", "CCT", 
#                                           "CGA", "CGT", 
#                                           "CTA", "CTC", "CTT"))


both_2 <- both[both$subcontext != "CGA", ]
both_2 <- both_2[both_2$subcontext != "CGT", ]

p <- ggplot(both_2, aes(x=subcontext, y=prop_meth, fill=sample)) + 
  #geom_violin() + 
 # geom_bar(stat="identity") + 
 # geom_boxplot(outlier.size = 0, outlier.alpha = 0) + 
  stat_summary(
    fun.data = mean_cl_normal, position=position_dodge(.9), geom="errorbar", width = .25) + 
  stat_summary(fun = mean,  geom="bar", position = "dodge2", color="black") + 
  theme_classic() + 
  scale_y_continuous(expand=c(0,0)) + 
  ylab("percent methylation") + 
  ggtitle(title)

p

ggsave(paste0("subcontexts/",title, "4.png"))
ggsave(paste0("subcontexts/",title, "4.svg"))
