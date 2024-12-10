# installing/loading the package:
if(!require(installr)) {
  install.packages("installr"); require(installr)} #load / install+load installr
# using the package:
updateR() # this will start the updating process of your R installation.  It will check for newer versions, and if one is available, will guide you through the decisions you'd need to make.



#percent root regeneration graphs
library(ggplot2)

df <- data.frame(dose=c("D0.5", "D1", "D2"),
                 len=c(4.2, 10, 29.5))
df

root_regen <- read_excel("C:/Users/nsmoo/Desktop/Grad/root regen.xlsx")
root_regen


Root_regeneration <- read_excel("C:/Users/nsmoo/Desktop/Grad/Regeneration/Summer 23/Root regeneration.xlsx")
grey7 = "#757575"
teal0 = "#e0f2f1"
teal1 = "#b2dfdb"
teal2 = "#7fcbc4"
teal3 = "#7fcbc4"
teal4 = "#25a69a"
teal5 = "#009688"



#this is the working one
p <- ggplot(data=Root_regeneration, aes(x=dac, y=percent)) +
  geom_line(aes(color = genotype)) +
  scale_color_manual(values=c(grey7, down_blue, up_red)) + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0,0)) + 
  ggtitle("Percent of explants that regenerate roots")


p + theme_classic()
ggsave("root regeneration.png", width=5.41, height=4.81)
