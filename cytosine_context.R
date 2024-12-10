library(ggplot2)



subcontexts <- c("CAG", "CTG", "CCG", "CAG", "CTG", "CCG")
genotypes <- c(rep("all Cs", 3), rep("dmcs", 3))

values <- c(42, 42, 16, 59, 36, 4)
df <- data.frame(subcontexts, genotypes, values)
df

p <- ggplot(df, aes(x = genotypes, y = values, fill = subcontexts)) + 
  geom_bar(position="fill", stat="identity", color = "black") + 
  theme_classic() + 
  xlab("CHG in 4-1 R0") + ylab("proportion") + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = c("#4DAADA", "#BDDD4B", "#DB4D97"))

p
ggsave("4-1 CHG subcontexts.svg")
