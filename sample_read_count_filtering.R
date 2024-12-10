cytes <- c(1, 2, 3)
read_counts <- c(10, 1, 30)
meth <- c(10, 10, 10)

df_counts <- data.frame(cytes, read_counts, meth)


df_counts
df_final <- data.frame(matrix(NA, nrow=1, ncol=3))
colnames(df_final) <- c("cytes", "read_counts", "meth")

for (row in seq_along(df_counts)) {
  row_working <- df_counts[row,]
  print(row_working)
  if (row_working$read_counts < 5) {
    row_working$meth <- NA
  }
  df_final <- rbind(df_final, row_working)
}

df_counts
df_final
df_final_fr <- df_final[-1,]

df_heatmap <- data.frame(df_final_fr$meth)
