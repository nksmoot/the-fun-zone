#Hypocotyl regeneration experiments

library(ggplot)
#Number of hypocotyls with 1+ adventitious roots on induction media (hormone free) at 21dpt

geno <- c("a DRDD", "b dme", "c rdd", "d drdd")
values <- c(2.8, 10, 7.1, 63)

df <- data.frame(geno, values)
p <- ggplot(df, aes(x=geno, y=values, fill=geno)) + 
  geom_bar(stat="identity", color="black") + 
  theme_classic() + 
  scale_y_continuous(expand=c(0,0)) + 
  ggtitle("Hypocotyls with 1+ adv roots at 21dpt on IM") +
  ylab("Percent") + xlab("Genotype") + 
  scale_fill_manual(values=c("#FBE1C0", "#F8EFA2", "#E9F271","#BDDD4B"))
p
ggsave("/Users/nksmoot/Desktop/Grad/r studio things/Manuscript/Fig 1/adv roots hypes IM at 21dpg.svg", 
       width = 2, height = 3, units = "in")



library(readxl)
counts_roots <- read_excel("/Users/nksmoot/Desktop/Grad/r studio things/counts of roots c8 i13.xlsx")
short <- counts_roots[,1:4]
head(short)

library(reshape2)

melted <- melt(short)
head(melted)
melted <- na.omit(melted)

p <- ggplot(melted, aes(x=variable, y=value, fill=variable)) + 
  stat_summary(
    fun.data = mean_cl_normal, position = position_dodge(), geom="errorbar", width = .25) + 
  stat_summary(fun = mean,  geom="bar", color = "black") + 
  theme_classic() + 
  #geom_jitter(size = .1, alpha = .3, width = .15) + 
  scale_y_continuous(expand=c(0,0)) + 
  scale_fill_manual(values=c("#FBE1C0", "#F8EFA2", "#E9F271","#BDDD4B")) + 
  ggtitle("Number of adventitious roots per hypocotyl\n8 days CIM, 13 days IM")

p
ggsave("/Users/nksmoot/Desktop/Grad/r studio things/Manuscript/Fig 1/adv roots c8 i13.svg", 
       width = 2, height = 3, units = "in")



#obviously green hyps
number <- c(9, 28.5, 23.8, 30.7)
df <- data.frame(geno, number)
df
p <- ggplot(df, aes(x=geno, y=number, fill=geno)) + 
  geom_bar(stat="identity", color="black") + 
  theme_classic() + 
  scale_y_continuous(expand=c(0,0), limits=c(0, 40)) + 
  ggtitle("Hypocotyls obviously green at 21dpt on SIM") +
  ylab("Percent") + xlab("Genotype") + 
  scale_fill_manual(values=c("#FBE1C0", "#F8EFA2", "#E9F271","#BDDD4B"))
p
ggsave("green hypes SIM at 21dpg.svg")



##fake data to hold space: 
genotypes <- c(rep("DRDD", 4), rep("dme", 4), rep("rdd", 4), rep("drdd", 4))
percents <- c(80, 10, 10, 0, 50, 30, 20, 0, 40, 20, 20, 20, 0, 15, 45, 40)
types <- c(rep(c("T1", "T2", "T3", "T4"), 4))
df <- data.frame(genotypes, percents, types)
df

ggplot(df, aes(x=genotypes, y=percents, fill=types)) + 
  geom_bar(stat="identity", color="black") + 
  theme_classic() + 
  scale_y_continuous(expand = c(0,0))

ggsave("placeholder_cim_sim.svg")



#SIM direct
genotypes <- c(rep("a DRDD", 4), rep("b dme", 5), rep("c rdd", 6), rep("d drdd", 5))
shoots <- c(0, 
0, 
0, 
0, 
0.2,
0,
0.022727273,
0.085714286,
0.066666667,
0.333333333,
0.078125,
0.068965517,
0.291666667,
0.057142857,
0,
0.25,
0.054054054,
0.195121951,
0,
0.058823529)
df <- data.frame(genotypes, shoots)

p <- ggplot(df, aes(x=genotypes, y=shoots, fill=genotypes)) + 
  stat_summary(
    fun.data = mean_cl_normal, position = position_dodge(), geom="errorbar", width = .25) + 
  stat_summary(fun = mean,  geom="bar", color = "black") + 
  theme_classic() + 
  #geom_point() + 
  #geom_jitter(size = .1, alpha = .3, width = .15) + 
  scale_y_continuous(expand=c(0,0)) + 
  scale_fill_manual(values=c("#FBE1C0", "#F8EFA2", "#E9F271","#BDDD4B")) + 
  ggtitle("Percent of hyps with +1 shoot\nSIM only")

p
ggsave("/Users/nksmoot/Desktop/Grad/r studio things/Manuscript/Fig 1/SIM direct shoots.svg", 
       width = 2, height = 3, units = "in")




#number DEGs
delta <- c("down", "up")
number <- c(360, 280) 
df <- data.frame(delta, number)

ggplot(df, aes(x=delta, y=number, fill=delta)) + 
  geom_bar(stat = "identity", color="black") + 
  theme_classic()

ggsave("/Users/nksmoot/Desktop/Grad/r studio things/Manuscript/Fig 1/DEGS.svg", 
       width = 2, height = 3, units = "in")




#DNSR placeholder
values <- c(0,0,0,1,3)
gt <- c("WT", "dme", "rdd", "drdd", "R1")
df_regen <- data.frame(gt, values)

ggplot(df_regen, aes(x=gt, y=values)) + 
  geom_bar(stat="identity") + 
  theme_classic()


ggsave("/Users/nksmoot/Desktop/Grad/r studio things/Manuscript/Fig 1/dnsr.svg", 
       width=3, height=3, units="in")
