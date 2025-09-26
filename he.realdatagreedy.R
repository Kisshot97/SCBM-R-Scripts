# -------------------
# 仍然使用你的读取方式：getwd() + "/ACTGout/..." + paste0()
# -------------------
library(reshape2)
library(foreach)
library(Rmisc)
library(ggplot2)
library(grid)
library(cowplot)
library(gridExtra)
library(dplyr)
## 1) 只用 boost = 100
g_boost_values <- 100
g_boost <- g_boost_values
g_sam   <- g_sam_values  # 保持不变

## ... 你现有的循环与 df_long 生成代码保持不变 ...
g_boost <- g_boost_values
g_sam   <- g_sam_values

# 收集 out1/out2/out3 的每个 scenario 的 MSE，做成长表
rows <- list()

for (boost in g_boost) {
  for (sam in g_sam) {
    for (i in 1:10) {
      # ----- 保持你的文件读取方式不变 -----
      addr0 <- getwd()
      addr_save  <- sprintf("/ACTGout/scenario%s_sam%s_boost%s.csv", i, sam, boost)
      output_path <- paste0(addr0, addr_save)
      # ----------------------------------
      
      if (!file.exists(output_path)) next
      
      dat0 <- read.csv(output_path)
      
      # 计算每个 HTE 的 MSE（与 true 的平方误差和）
      mse1 <- sum((dat0$hte1 - dat0$true)^2, na.rm = TRUE)
      mse2 <- sum((dat0$hte2 - dat0$true)^2, na.rm = TRUE)
      mse3 <- sum((dat0$hte3 - dat0$true)^2, na.rm = TRUE)
      
      # 存入长表（标记为 K=1/2/3）
      rows[[length(rows) + 1]] <- data.frame(
        boost = boost, sam = sam, scenario = i, Kmax = "K = 1", MSE = mse1
      )
      rows[[length(rows) + 1]] <- data.frame(
        boost = boost, sam = sam, scenario = i, Kmax = "K = 2", MSE = mse2
      )
      rows[[length(rows) + 1]] <- data.frame(
        boost = boost, sam = sam, scenario = i, Kmax = "K = 3", MSE = mse3
      )
    }
  }
}

df_long <- if (length(rows)) do.call(rbind, rows) else stop("未找到任何匹配的 CSV 文件。")

## 2) 画折线图（按 scenario 求每个 K 的 MSE 的 “中央値 + IQR”）
# 仅取 boost=100 的数据（前で factor 化していれば文字でOK）
df100 <- subset(df_long, boost == "100")

## 集計：中央値・Q1・Q3
stat_raw <- aggregate(
  MSE ~ Kmax, data = df100,
  FUN = function(x) c(
    med = median(x, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE)
  )
)



mse_stat <- df100 %>%
  group_by(Kmax) %>%
  summarise(
    med  = median(MSE, na.rm = TRUE),
    q25  = quantile(MSE, 0.25, na.rm = TRUE),
    q75  = quantile(MSE, 0.75, na.rm = TRUE)
  ) %>%
  ungroup()
## 型・順序を整える
mse_stat$Kmax <- factor(mse_stat$Kmax, levels = c("K = 1","K = 2","K = 3"))

## 折れ線（中央値）＋ IQR エラーバー
p_line <- ggplot(mse_stat, aes(x = Kmax, y = med, group = 1)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.15, linewidth = 0.7) +
  labs(
    title = "Median Squared Error with Quartile Range (Boost = 100)",
    x = "Kmax", y = "Median Squared Error"
  ) +
  ylim(0, 15) +
  theme(
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10)
  )

print(p_line)

ggsave(
  filename = paste0("D:/personal/work/test326/ACTGout/pic/", "He.greedy.boost100.eps"),
  plot = p_line,
  width = 6.5, height = 4, device = cairo_ps
)

