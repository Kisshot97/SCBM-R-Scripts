################################################################################
# Module: Greedy Plot
# Description: This script is part of the SCBM framework and implements
# core functionality related to greedy plot.
################################################################################


# ====================
# Load Required Libraries
# ====================
library(reshape2)

# ====================
# Load Required Libraries
# ====================
library(foreach)

# ====================
# Load Required Libraries
# ====================
library(Rmisc)

# ====================
# Load Required Libraries
# ====================
library(ggplot2)

# ====================
# Load Required Libraries
# ====================
library(grid)

# ====================
# Load Required Libraries
# ====================
library(cowplot)

# ====================
# Load Required Libraries
# ====================
library(gridExtra)

quan=function(x) {
  ul=quantile(x, c(0.05,0.95))
  low=which(x < ul[1])
  up=which(x > ul[2])
  x[low]=ul[1]
  x[up]=ul[2]
  #x=x[-c(low,up)]
# --- Return Final Output ---
  return(x)}
######input######
# Define parameters
g_boost <- seq(20, 200, 20)
g_sam <- seq(0.2, 0.8, 0.1)

# Initialize storage for the results
out1_df <- matrix(NA, nrow = length(g_sam), ncol = length(g_boost), 
                  dimnames = list(paste0("sam", g_sam), paste0("boost", g_boost)))
out2_df <- out1_df
out3_df <- out1_df
out4_df <- out1_df

# Iterate over g_boost and g_sam
for (boost in g_boost) {
  for (sam in g_sam) {
    out1=out2=out3=out4=numeric()
    for (i in 1:10) {
      # Construct file path
      output_path <- sprintf(file.path(getwd(), "scenario%s_sam%s_boost%s.csv"), i, sam, boost)
      
      # Check if the file exists
      if (file.exists(output_path)) {
        # Read the data

# --- Load CSV Data ---
        dat0 <- read.csv(output_path)
        
        # Calculate mse for each hte
        out1[i]  <- sum((dat0$hte1 - dat0$true)^2)
        out2[i]  <- sum((dat0$hte2 - dat0$true)^2)
        out3[i]  <- sum((dat0$hte3 - dat0$true)^2)
        out4[i]  <- sum((dat0$hte4 - dat0$true)^2)
      }
    }
      out1=sum(out1)
      out2=sum(out2)
      out3=sum(out3)
      out4=sum(out4)
      print(out1)
      # Store results in respective matrices
      out1_df[paste0("sam", sam), paste0("boost", boost)] <- out1
      out2_df[paste0("sam", sam), paste0("boost", boost)] <- out2
      out3_df[paste0("sam", sam), paste0("boost", boost)] <- out3
      out4_df[paste0("sam", sam), paste0("boost", boost)] <- out4
    }
  }


# --------------------------------------------------
# Function: plot_heatmap1
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
plot_heatmap1 <- function(data, title) {
  data_melt <- melt(data, varnames = c("sam", "boost"), value.name = "CVerror")
  data_melt$sam <- as.numeric(gsub("sam", "", data_melt$sam))
  data_melt$boost <- as.numeric(gsub("boost", "", data_melt$boost))
  
  ggplot(data_melt, aes(x = boost, y = sam, fill = CVerror)) +
    geom_tile() +
    scale_fill_gradient(low = "#FFF6E9", high = "#0B60B0") +  # 红色为低值，蓝色为高值
    labs(title = title,
         x = "Boost",
         y = "Sampling Proportion",
         fill = "CVerror") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}



# --------------------------------------------------
# Function: plot_heatmap
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
plot_heatmap <- function(data, title) {
  data_melt <- melt(data, varnames = c("sam", "boost"), value.name = "CVerror")
  data_melt$sam <- as.numeric(gsub("sam", "", data_melt$sam))
  data_melt$boost <- as.numeric(gsub("boost", "", data_melt$boost))
  
  ggplot(data_melt, aes(x = boost, y = sam, fill = CVerror)) +
    geom_tile() +
    scale_fill_gradient(low = "#FFF6E9", high = "#0B60B0",limits = c(90, 95), oob = scales::squish) +  # 红色为低值，蓝色为高值
    labs(title = title,
         x = "Boost",
         y = "Sampling Proportion",
         fill = "CVerror") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# 使用函数绘制四张热力图
p1 <- plot_heatmap1(out1_df, "Heatmap of MSE for Kmax=1")
p2 <- plot_heatmap(out2_df, "Heatmap of MSE for Kmax=2")
p3 <- plot_heatmap(out3_df, "Heatmap of MSE for Kmax=3")
p4 <- plot_heatmap(out4_df, "Heatmap of MSE for Kmax=random in 1-3")

# 使用 grid.arrange 将四张图排列在一起
ppout=grid.arrange(p1, p2, p3, nrow = 3)
# 保存为 EPS 文件

# --- Save Plot to File ---
ggsave(
  filename = paste0(file.path(getwd(), ""), "greed000.eps"), # 修改扩展名为 .eps
  plot = ppout, # 需要保存的图表对象
  width = 8,    # 图像宽度
  height = 10,  # 图像高度
  device = cairo_ps # 使用 Cairo PostScript 设备生成 EPS 文件
)
