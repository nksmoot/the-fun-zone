#volcano plot drdd DEG
library(ggplot2)
install.packages("ggsci")
library(readxl)
library(grid)
WT_v_drdd_21dpg <- read_excel("C:/Users/nsmoo/Desktop/Grad/RNA seq data/WT vs drdd 21dpg.xlsx")
WT_v_drdd_13dpg <- read_excel("C:/Users/nsmoo/Desktop/Grad/RNA seq data/WT v drdd 13dpg.xlsx")
WT_v_drdd_30dpg <- read_excel("C:/Users/nsmoo/Desktop/Grad/RNA seq data/WT v drdd 30dpg.xlsx")
WT_v_dR_21dpg <- read_excel("C:/Users/nsmoo/Desktop/Grad/RNA seq data/WT vs dR 21dpg.xlsx")
WT_v_drdd_21dpg_2 <- read_excel("C:/Users/nsmoo/Desktop/Grad/RNA seq data/WT vs drdd 21dpg 202311.xlsx")


df <- WT_v_drdd_21dpg_2


downregulated <- df[1, 14]
upreglated <- df[2, 14]

print(downregulated)
print(upreglated)

grob <- grobTree(textGrob(downregulated, x=.8, y=.95, hjust=0,
                          gp=gpar(col="#0377bd", fontsize=13)))
grob2 <- grobTree(textGrob(upreglated, x=.8, y=.85, hjust=0,
                          gp=gpar(col="#c62828", fontsize=13)))

#hello <- c("hi", "hello", "gutn tag")
#goodbye <- c("toodles", "see ya", NA)
#df <- data.frame(hello, goodbye)

title <- "21dpg_2"

p <- ggplot(data=df, aes(x=log2FoldChange, y=lg10, col=DEG))+
  geom_point() +
  scale_color_manual(values=c("#0377bd", "grey", "#c62828")) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("-log10(p value)") +
  ggtitle(title)
  

p + theme_classic() + annotation_custom(grob) + annotation_custom(grob2)

ggsave(paste("DEGs WT vs drdd", title, ".png"), width=3.77, height=3.9)


#proximity to hyperMR
#45 out of 113 downregulated genes withing 2kb of a hyperMR
#34 out of 205 upregulated genes within 2kb of a hyperMR
#random chance looks like.. 12% 

genotype <- c("adown", "bup", "cdown", "dup")
value <- c(39.8, 16.5, 12.2, 12.3)
df <- data.frame(genotype, value)
df

p <- ggplot(data=df, aes(x=genotype, y=value, fill=genotype)) + 
  geom_bar(stat="identity", color="#000000", size=1) + 
  scale_fill_manual(values=c("#0377bd",
                             "#c62828",
                             "grey",
                             "grey")) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 50)) 

p + theme_classic() + theme(axis.line=element_line(size=1))
ggsave("hyperMR.png")


teal0 = "#e0f2f1"
teal1 = "#b2dfdb"
teal2 = "#7fcbc4"
teal3 = "#7fcbc4"
teal4 = "#25a69a"
teal5 = "#009688"

sage1 = "#DFEBE3"
yellow0 = "#FDF3E1"
yellow1 = "#FBE7C3"
yellow2 = "#F9DBA5"
yellow3 = "#F5C369"

pink0 = "#fce4ec"
pink1 = "#f8bbd0"
pink2 = "#f48fb1"

purple0 = "#f3e5f5"
purple1 = "#e1bee7"
purple2 = "#ce93d8"



down_blue = "#0377bd"
up_red = "#c62828"


#number of primary roots assay
num_roots <- c("3+", "2", "1", "0")
values <- c(4, 8, 14, 1, 22, 0, 0, 0)
genotype <- c(rep("aWT", 4), rep("drdd", 4))
df <- data.frame(genotype, num_roots, values)
df

p <- ggplot(df, aes(fill=factor(num_roots, levels=c("3+", "2", "1", "0")), y=values, x=genotype)) +
  geom_bar(position="fill", stat="identity", color = "black", width = .5) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=c(yellow3, 
                            yellow2, 
                             yellow1,
                            yellow0)) + 
  labs(color = "Number of primary roots")

p + theme_classic() + theme(axis.line=element_line(size=1)) 
ggsave("primaryroots.png")



#percent of explants with roots that regenerate shoots
percent_shoots <- c(0, 4.3, 12.7)
genotype <- c("aWT", "drdd", "R-drdd")
df_fromroot <- data.frame(genotype, percent_shoots)
df_fromroot

p <- ggplot(df_fromroot, aes(y=percent_shoots, x=genotype, fill=genotype))+
  geom_bar(stat="identity", color = "black", width = .75) + 
  scale_y_continuous(expand = c(0, 0), limits=c(0, 15)) + 
  scale_fill_manual(values=c(yellow0, 
                             yellow2, 
                             yellow3)) + 
  ylab("%explants that had roots that make shoots")

p + theme_classic() + theme(axis.line=element_line(size=1))
ggsave("percentshoots_fromroots.png", width=4.23, height=4.18)


#percent of explants total that regenerate shoots
percent_shoots <- c(0, 1.09, 4.6)
genotype <- c("aWT", "drdd", "R-drdd")
df_total <- data.frame(genotype, percent_shoots)
df_total

p <- ggplot(df_total, aes(y=percent_shoots, x=genotype, fill=genotype))+
  geom_bar(stat="identity", color = "black", width = .75) + 
  scale_y_continuous(expand = c(0, 0), limits=c(0, 8)) + 
  scale_fill_manual(values=c(yellow0, 
                             yellow2, 
                             yellow3)) + 
  ylab("%explants that make shoots")

p + theme_classic() + theme(axis.line=element_line(size=1))

ggsave("percentshoots.png", width=4.23, height=4.18)

#by plate, total, with error bars, just the mean graphed in bars
genotype <- c("aWT", "drdd", "R-drdd")
mean_percent_total <- c(0, 1.43, 4.49)
sd_total <- c(0, 2.02, 3.24)
df_total_plate <- data.frame(genotype, mean_percent_total, sd_total)
df_total_plate

p <- ggplot(df_total_plate, aes(x=genotype, 
                                y=mean_percent_total, 
                                fill=genotype)) + 
  geom_bar(stat="identity", color = "black", width = .75) +
  scale_fill_manual(values=c(yellow0, 
                             yellow2, 
                             yellow3)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 8)) +
  geom_errorbar(aes(ymin=mean_percent_total, 
                    ymax=mean_percent_total+sd_total),
                width=.2)

p + theme_classic()
ggsave("percent with shoots from total by plate with error bars.png")

install.packages("survival")
install.packages("lattice")
install.packages("ggplot2")
install.packages("Hmisc")


#by plate, total, dot plot
genotype <- c(rep("aWT", 2), rep("drdd", 2), rep("R-drdd", 5))
percent_total <- c(0, 0, 0, 2.86, 0, 5.56, 2.78, 8.57, 5.56)
df <- data.frame(genotype, percent_total)
df


library(dplyr)
df.summary <- df %>%
  group_by(genotype) %>%
  summarise(
    sd = sd(percent_total, na.rm = TRUE),
    percent_total = mean(percent_total)
  )
df.summary
geom_errorbar(aes(ymin=mean_percent_total, 
                  ymax=mean_percent_total+sd_total),
              width=.2)


p <- ggplot(df, aes(x=genotype, y=percent_total)) + 
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize=.75, color=grey7, fill=grey7) + 
  
  geom_errorbar(aes(ymin=percent_total-sd, 
                      ymax=percent_total+sd), 
                    width=.2, color=up_red,
                  data=df.summary) + 
  stat_summary(fun.y=mean, geom="point", color=up_red) +
  ggtitle("percent of explants that make shoots, \nfrom the total number, grouped by plate") 
  
p + theme_classic() + theme(plot.title = element_text(size=12))
ggsave("dotplot percent total by plate.png", width=3.86, height=4.81)


--------
p <- ggplot(df, aes(x=genotype, y=percent_total)) + 
  geom_dotplot(binaxis='y', stackdir='center')


p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="red", width=0.2) +
  stat_summary(fun.y=mean, geom="point", color="red") + 
  theme_classic()


ggsave("dots_no_error_bars_lol.png")

grey7 = "#757575"
