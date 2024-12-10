---
  title: "analysis"
output: html_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## SPCH expression

## notes

```{r}
library(ggplot2)
library(reshape2)
library(ggpubr)
install.packages("RColorBrewer")                   # Install RColorBrewer package
library("RColorBrewer")                            # Load RColorBrewer

my_theme <-  theme_classic(6) + theme(
  legend.position = "none",
  axis.text.x = element_text(angle = 60, hjust = 1),  #rotate & align axis
  axis.text = element_text(size = 6),
  axis.title.x = element_blank(),
  axis.line = element_line(size = 0.2),
  axis.ticks = element_line(size = 0.2),
  strip.background = element_blank(),
  strip.placement = "outside"
)


```
mydata <- read.csv("C:\\Users\\nsmoo\\Desktop\\Grad\\PMB10 midterm f22.csv")

df_whole <- data.frame(mydata)
print(df_whole)
print(colnames(df_whole))
df_q1  <- df_whole$q1
df_totals <- df_whole$ï..total
print(df_totals)

(28, 11, 11, 12, 4, 8, 10)
(1, 2, 3, 4, 5, 6, 7)
2.608465049	44.70464989

df_spch <- data.frame(genotype = c("WT", "drdd"), level = c(2.6, 44.7))
df_skip <- data.frame(questions = c(1, 2, 3, 4, 5, 6, 7), 
                   skip = c(28, 11, 11, 12, 4, 8, 10))
print (df_skip)
print (df_q1)

hist(df_q1, 
     breaks = 20,
     main="Distribution of answers for question 1",
     freq=FALSE
)
#teals: color = "#00796b", fill = "#7fcbc4"
#cyans: color = "#0097a7", fill = "#b2dfdb"
#greys: color = "#636363", fill = "#C9C9C9"

library(ggplot2)
number <- "7"
#holder <- paste("Question ", number)
holder <- "Total Score"
saveas <- gsub(" ", "", paste(holder, ".png"))

p<-ggplot(df_whole, aes(x=ï..total)) + 
  geom_histogram(binwidth=5, color = "#00796b", fill = "#7fcbc4") +
  xlab("Distribution of Scores") + 
  #ylim(0,20) + 
  #xlim(0, 100) + 
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(holder)
p + theme_classic() 
ggsave(saveas)

p<-ggplot(data = df_spch, aes(x=reorder(genotype, level), y = level)) + 
  geom_bar(stat = "identity", color = "#003350", fill = "#005485") +
  xlab("Genotype") + 
  #ylim(0, 30) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("SPCH expression")
p + theme_classic() 
ggsave("spch expn.png", width = 1.5, height = 2)



specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)
print(data)

geno <- c(rep("DRDD 10dpex", 3), rep("drdd 10dpex", 3), rep("DRDD 17dpex", 3), rep("drdd 17dpex", 3))
#geno2 <- c(rep("DRDD 17dpex", 3), rep("drdd 17dpex", 3))
growth <- c("nothing", "callus", "root", "nothing", "callus", "root")
ten_dpex <- c(13, 4, 0, 7, 5, 7, 12, 5, 0, 3, 6, 10)
#seventeen_dpex <- c(12, 5, 0, 3, 6, 10)
#root <- c(7, 1)
df_regen <- data.frame(geno, growth, ten_dpex)
print(df_regen)

ggplot(df_regen, aes(fill=growth, y=ten_dpex, x=geno)) + 
  geom_bar(position="fill", stat="identity")


library(ggplot2)
p <- ggplot(df_regen, aes(fill=factor(growth, levels=c("root", "callus", "nothing")), 
                          y=ten_dpex, x=geno)) + 
  geom_bar(position="fill", stat="identity", width = .5) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("10 and 17 days post explant")+ 
  labs(fill="Growth") +
  scale_fill_manual(values=c("#7fcbc4", "#b2dfdb", "#C9C9C9" ))
p + theme_classic() + scale_x_discrete(limits=c("DRDD 10dpex", "DRDD 17dpex", "drdd 10dpex", "drdd 17dpex"))
ggsave("seventeen days post explant.png")


geno <- c(rep("WT young", 2), rep("WT old", 2), rep("drdd young", 2), rep ("drdd old", 2))
ct <- c(1.04, 1.04, 1.12, 0.89, 3.21, 0.68, 1.33, 0.82)
df_wox5 <- data.frame(geno, ct)
print(df_wox5)
p <- ggplot(df_wox5, aes(x=geno, y=ct))+
  geom_bar(stat="identity")
p +theme_classic()




p <- ggplot(df_wox5, aes(geno, ct)) +
  geom_jitter(position = position_jitter(0.2), color = "#0097a7") +
  ggtitle("WOX5 expression")+
  ylab("Relative expression") +
  xlab("Genotype and Age") +
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="#b2dfdb")
  
p + theme_classic()
ggsave("wox5 2022 10 20.png")

#cyans: color = "#0097a7", fill = "#b2dfdb"
#green: color = "#388e3c"
#teals: color = "#00796b", fill = "#7fcbc4"
install.packages("gghighlight")
library(gghighlight)
library(ggplot2)

age <- "13"
title <- paste(paste("WT vs drdd at ", age), " dpg")
df_13 <- data.frame(WT_v_drdd_13dpg)
head(df_13)
p <- ggplot(df_13, aes(x=log2FoldChange)) + 
  geom_histogram(binwidth=.01, color = "#0097a7")+
  xlim(-10, 10)+
  ggtitle(title)+ 
  xlab("Log 2 Fold Change") + ylab("Number DEGs")+
  scale_y_continuous(expand = c(0, 0))+
  gghighlight(count > 400, label_key = log2FoldChange)
p + theme_classic()
ggsave(paste(title, ".png"))

num <- c(1:5)
genotype <- c("DRDD", "rdd", "drdd 4g4", "drdd 6g3", "drdd R1")
rootpercm <- c(0.909, 1.499, 0.501, 1.506, 0.08)
avglen <- c(1.02, 0.627, 0.798, 0.956, 0.47)

df <- data.frame(num, genotype, rootpercm, avglen)

print (df)

df1 <- df
df1$x <- factor(df1$x, levels = c("DRDD", "rdd", "drdd 4g4", "drdd 6g3", "drdd R1"))
title <- "Avg length of hypocotyl"

p <- ggplot(data=df, aes(x=genotype, y=avglen)) +
  geom_bar(aes(x=factor(genotype, genotype), y=avglen), 
               color = "#0097a7", fill = "#b2dfdb", 
           stat="identity", show.legend = FALSE) + 
  xlab("Genotype") + ylab("Avg length (cm)")+
  ggtitle(title)+
  scale_y_continuous(expand = c(0, 0))
  

p + theme_classic() 
ggsave(paste(title, ".png"))



#graph WOX5 rna seq reads
library(ggplot2)

# create a dataset
genotype <- c(rep("a WT a", 3), rep("drdd", 3))
age <- c("13 days", "21 days", "30 days","13 days", "21 days", "30 days")
reads <- c(1.9, 0, .3, 16.9, 4.9, 1.1)
data <- data.frame(genotype, age, reads)
print(data)

title <- "WOX5 expression"
    
outline="#000000"
color0 = "#b2ebf2"
color1 = "#80deea" #light teal
color2 = "#4dd0e1" #med teal
color3 = "#00bcd4" #dark teal

 
color1 = "#ffcdd2" #light red
color2 = "#f8bbd0" #med teal 
color3 = "#e1bee7" #dark teal


# Grouped
p <- ggplot(data, aes(x=genotype, y=reads, fill=age)) + 
  geom_bar(stat="identity", color=outline, position="dodge", size=1) + 
  scale_y_continuous(expand = c(0, 0)) + 
  ggtitle(title) + 
  scale_fill_manual(values = c(color1, color2, color3)) +
  xlab(" ") + ylab("Read Count")
  

p + theme_classic() + 
ggsave(paste(title, ".png"))
