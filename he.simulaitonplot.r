library(reshape2)
library(Rmisc)
library(ggplot2)
library(grid)
library(cowplot)
library(gridExtra)

for (sim in 1:6) {
  for (deg in 1:3) {
    for (tp in c("bias","mse")) {
      
      
      pl.type=tp#########bias or mse#########
      
      param_pairs <- data.frame(
        nn = c(1, 2, 3, 3),
        pp = c(1, 2, 3, 4)
      )
      
      grid <- merge(
        param_pairs,
        expand.grid(ddad = 1:12, KEEP.OUT.ATTRS = FALSE),
        by = NULL,
        all = TRUE
      )       
      npnp=t(grid)
      datalist=bartlist=list()
      nop=9 # number of p
      for (i in 1:nop) {
        datalist[[i]]=read.csv(sprintf("D:/personal/work/test326/out/pdat/%s%s_deg%s_var%s.csv",pl.type, sim, deg, i))[-1]######input######
      }
      dat=list()
       
      for (np in 1:100){#####extract first 1000 samples
        aa=list()
        for (k in 1:nop) {
          aa[[k]]=na.omit(na.omit(datalist[[k]][np])[1:ncol(npnp),])
        }
        aa=as.data.frame(as.list(aa))
        colnames(aa)=c("pr","pr_ps","pr_stra","TO","BCM","BCM0","Cbart","CF","VT")[1:nop]
        dat[[np]]=aa
      }
      # Assuming `dat` is a list of 100 data.frames with 48 rows and 8 columns each
      
      # Step 1: Calculate mean matrix
      mean_mat <- Reduce("+", dat) / length(dat)
      
      # Step 2: Calculate standard error matrix
      se_mat <- array(0, dim = dim(mean_mat))
      
      for (i in 1:48) {
        for (j in 1:nop) {
          values <- sapply(dat, function(df) df[i, j])
          se_mat[i, j] <- sd(values) / sqrt(length(values))
        }
      }
      
      # Step 3: Format as "mean (SE)"
      formatted_mat <- matrix("", nrow = 48, ncol = nop)
      for (i in 1:48) {
        for (j in 1:nop) {
          formatted_mat[i, j] <- sprintf("%.5f (%.5f)", mean_mat[i, j], se_mat[i, j])
        }
      }
      
      # Step 4: Convert to data frame and set column names
      formatted_df <- as.data.frame(formatted_mat)
      colnames(formatted_df) <- colnames(dat[[1]])
      
      # Generate row names
      sims0 <- 1:12
      scena0 <- c("n=200,p=25", "n=400,p=50", "n=1000,p=100", "n=1000,p=150")
      
      row_labels <- character()
      for (sims in sims0) {
        for (scena in scena0) {
          row_labels <- c(row_labels, paste0("sim_", sims, "_scenario_", scena))
        }
      }
      
      # Assign to formatted_df
      rownames(formatted_df) <- row_labels
      # Define the file path
      filename0 <- paste0("D:/personal/work/test326/pic/", pl.type, "_sim_", sim, "_deg_", deg, ".csv")
      # Write to CSV
      write.csv(formatted_df, file = filename0, row.names = TRUE)
      
      dat1=dat
      mdat=melt(dat1)
      mdat=na.omit(mdat)
      mdat[,2]=abs(mdat[,2])
      mdat=mdat[-which(mdat$variable=="TO"),]####delete#####nops
      pic1=list()
      
      
      for (i in 1:12) {
        ni=4*(i-1)+1
        mf1=mdat[mdat[,3]==c(ni,ni+1,ni+2,ni+3),]
        colnames(mf1)=c("method","MSE","p")
        lable_name=list("n=200,p=25","n=400,p=50", "n=1000,p=100", "n=1000,p=150")
        lable_name_labeller=function(variable,value){
          return(lable_name[((as.numeric(value)-1)%%4+1)])
        }
        name_mf1=unique(mf1[,1]) # model picking!!!
        mf1=mf1[is.element(mf1[,1],name_mf1),] # model picking!!!
        ymax=quantile(na.omit(mf1[,2]),c(0,0.95))# y scale!!
        
        ######pic ##################
        pic1[[i]]=ggplot(mf1) +
          geom_boxplot(aes( y = MSE, x = method,fill=method)) +
          scale_y_continuous(trans = "sqrt", limits = ymax) +
          ylab(pl.type) +
          #scale_x_discrete(labels = c("n=100,p=50","n=200,p=100", "n=500,p=200", "n=500,p=400")) +
          xlab("scenario") +
          theme(legend.position = "none", text = element_text(size = 12), axis.title.x = element_text(size = 12),
                axis.title.y = element_text(size = 15), axis.text.x  = element_text(size = 8), 
                axis.text.y  = element_text(size = 12), plot.title = element_text(size = 12), 
                plot.margin = unit(c(0, 0, 0, 0), "cm")) +
          ggtitle(sprintf("Simulation%s", i)) +
          scale_fill_manual(values = c("#0B60B0","#B7C9F2", "#D96098", "#1FE5BD", "#EA5455", "#CDF2CA", "#FFDEB4", "#57C5B6","#87A922","#E69F00", "#EE99C2","#FFF6E9","#12372A")) +
          facet_wrap(.~factor(p),labeller=lable_name_labeller,strip.position = "bottom")
      }
      
      
      # Create a DataFrame
      data <- data.frame(
        Xdata = rnorm(nop), Ydata = rnorm(nop),
        Methods =c("PTO","CM","CF","VT","CBART","SCBM","prop_ps","TO","TO.S","LR","LL","Tmars","Tmars.ps")[1:nop])
      data$Methods= factor(data$Methods, levels =c("PTO","CM","CF","VT","CBART","SCBM","prop_ps","TO","TO.S","LR","LL","Tmars","Tmars.ps")[1:nop])
      #ptomse,cmmse,cfmse,bartmse,cbartmse,propmse,prop_psmse,TOmse,TO.Smse,LRmse,LLmse,Tmarsmse,Tmars.psmse
      # Create a Scatter Plot
      gplot <- ggplot(data, aes(Xdata, Ydata, fill = Methods)) +
        geom_boxplot(size = rep(1,nop)) +
        theme(
          legend.position = "bottom",
          text = element_text(size = 20),
          legend.text = element_text(size = 18),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 14)
        ) +guides(fill = guide_legend(nrow = 1))+
        scale_fill_manual(
          values = c("#0B60B0","#B7C9F2", "#D96098", "#1FE5BD", "#EA5455", "#CDF2CA", "#FFDEB4", "#57C5B6","#87A922","#E69F00", "#EE99C2","#FFF6E9","#12372A")
        )
      
      # Draw Only Legend without plot
      # Grab legend from gplot
      legend <- get_legend(gplot)
      
      library(gridExtra)
      
      for (i in c(1,7)) {
        # 将多个图表组合成一个页面
        ppout = grid.arrange(pic1[[i]], pic1[[i+1]], pic1[[i+2]], pic1[[i+3]], pic1[[i+4]], pic1[[i+5]], ncol = 3)
        abnm <- if (i == 1) "a" else "b"
        # 保存为 EPS 文件
        ggsave(
          filename = paste0("D:/personal/work/test326/pic/","He.", pl.type, "_sim_v2_deg_", deg, abnm, ".eps"),
          plot = ppout,
          width = 18,
          height = 10,
          device = cairo_ps  # 使用 Cairo PostScript 设备生成 EPS 文件
        )
      }
    }
  }
}