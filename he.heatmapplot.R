library(dplyr)
library(tidyverse)
library(ggplot2)
# Adjust path to where your VIMP files are
path <- "D:/personal/work/test326/ACTGout/" 


# Read and combine
n_files <- 20
# Load all VIMP files into a list
vimp_list <- lapply(1:n_files, function(i) {
  as=read.csv(sprintf("D:/personal/work/test326/ACTGout/VIMP.%d.csv", i),)
  colnames(as)=as[1,]
  as=as[-1,]
  as=as[,-1]
  return(as)
  
})
# Convert all to numeric data.frames (they're currently character)
vimp_list <- lapply(vimp_list, function(df) {
  as.data.frame(sapply(df, as.numeric))
})

# Now compute element-wise mean across all data.frames
vimp_mean_df <- Reduce("+", vimp_list) / length(vimp_list)


# Step 4: Normalize each row to [0, 100]
vimp_norm_df <- t(apply(vimp_mean_df, 1, function(row) {
  rng <- range(row, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(row)))  # constant row = 0
  (row - rng[1]) / diff(rng) * 100
}))
vimp_norm_df <- as.data.frame(vimp_norm_df)
colnames(vimp_norm_df)=colnames(vimp_mean_df)
# Step 5: Reorder columns by last row (descending)
last_row <- unlist(vimp_norm_df[nrow(vimp_norm_df), ])
ordered_cols <- order(last_row, decreasing = TRUE)
vimp_final_df <- vimp_norm_df[,ordered_cols ]
vimp_final_df <- vimp_final_df[c(nrow(vimp_norm_df),ordered_cols),]
rownames(vimp_final_df)=c("Full Data",colnames(vimp_final_df ))
#vimp_final_df[2,]=100
# Save result
vimp_norm=vimp_final_df
# === Step 2: Melt to long format (keep original order) ===
# === Step 2: Pivot THEN add Row and dash info ===
vimp_long <- vimp_norm %>%
  tibble::rownames_to_column(var = "Row") %>%
  pivot_longer(cols = -Row, names_to = "Col", values_to = "Score") %>%
  mutate(
    Row = factor(Row, levels = rev(rownames(vimp_norm))),
    Col = factor(Col, levels = colnames(vimp_norm)),
    RowID = as.numeric(Row),
    ColID = as.numeric(Col),
    dash = ifelse( RowID+ColID  == 13, "---", "")
  )

# === Step 3: Plot ===
picheat=ggplot(vimp_long, aes(x = Col, y = Row, fill = Score)) +
  geom_tile(color = "#FFF6E9") +
  geom_text(aes(label = dash), size = 6) +
  scale_fill_gradient(low = "#FFF6E9", high = "#EA5455", name = "VIMP\n[0–100]") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Normalized VIMP for Sensitivity Analysis",
    y = "Variable Removed",
    x = "Normalized VIMP (0–100)"
  )

ggsave(
  filename = paste0("D:/personal/work/test326/ACTGout/pic/", "he.greedy.heat.eps"),
  plot = picheat,
  width = 5, height = 4, device = cairo_ps
)

