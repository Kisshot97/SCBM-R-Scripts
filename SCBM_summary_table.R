################################################################################
# Module: Summary Table
# Description: This script is part of the SCBM framework and implements
# core functionality related to summary table.
################################################################################


# ====================
# Load Required Libraries
# ====================
library(reshape2)

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
 
  pl.type="bias"#########bias or mse#########
  
  nn0=c(1,2,3)
  pp0=c(1,2,3,4)
  npnp=numeric()
  nn=rep(nn0,length(pp0))
  pp=rep(pp0,length(nn0))
  npnp0=rbind(nn,pp)[,-c(4:11)]
  ddad0=c(1:12)
  ddad=rep(ddad0,each=ncol(npnp0))
  npnp=matrix(rep(npnp0, length(ddad0)), nrow = nrow(npnp0), byrow = FALSE)
  npnp=rbind(ddad,npnp)
  datalist=bartlist=list()
  nop=8
  for (i in 1:nop) {

# --- Load CSV Data ---
    datalist[[i]]=read.csv(sprintf(paste0(getwd(), "/pdat/%s1_res%s.csv"),pl.type,i))[-1]######input######
  }
  dat=list()
  qp <-as.numeric(c(1:(12*ncol(npnp0)))) 
  for (np in qp){#####extract first 1000 samples
    ddad=as.numeric(npnp[1,np])
    nn=as.numeric(npnp[2,np])
    pp=as.numeric(npnp[3,np])
    nx=c(200,400,1000)
    px=c(25,50,100,150)
    aa=list()
    for (k in 1:nop) {
      aa[[k]]=na.omit(na.omit(datalist[[k]][np])[1:100,])
    }
# --- Construct a Data Frame ---
    aa=as.data.frame(as.list(aa))
    colnames(aa)=c("prop","pr_ps","stra","TO","BCM","BCM0","CF","VT")[1:nop]
    dat[[np]]=aa
  }
  dat1=lapply(dat,function(df) {
  quan=function(x) {
    ul=quantile(x, c(0,0.95))
      low=which(x < ul[1])
      up=which(x > ul[2])
      x[low]=ul[1]
      x[up]=ul[2]
      #x=x[-c(low,up)]
# --- Return Final Output ---
      return(x)}
  df1=numeric()
  for (i in 1:ncol(df)) {
    df1=cbind(df1,quan(df[,i]))
  }
  colnames(df1)=colnames(df)
# --- Construct a Data Frame ---
  df1=as.data.frame(df1)
# --- Return Final Output ---
    return(df1)
  })
  #dat1=dat
  mdat=melt(dat1)
  mdat=na.omit(mdat)
  mdat[,2]=abs(mdat[,2])
  mdat=mdat[-which(mdat$variable=="TO"),]####delete#####nops
  pic1=list()
 #################
  pl.type="mse"
  for (i in 1:nop) {

# --- Load CSV Data ---
    datalist[[i]]=read.csv(sprintf(paste0(getwd(), "/pdat/%s1_res%s.csv"),pl.type,i))[-1]######input######
  }
  dat=list()
  qp <-as.numeric(c(1:(12*ncol(npnp0)))) 
  for (np in qp){#####extract first 1000 samples
    ddad=as.numeric(npnp[1,np])
    nn=as.numeric(npnp[2,np])
    pp=as.numeric(npnp[3,np])
    nx=c(200,400,1000)
    px=c(25,50,100,150)
    aa=list()
    for (k in 1:nop) {
      aa[[k]]=na.omit(na.omit(datalist[[k]][np])[1:100,])
    }
# --- Construct a Data Frame ---
    aa=as.data.frame(as.list(aa))
    colnames(aa)=c("prop","pr_ps","stra","TO","BCM","BCM0","CF","VT")[1:nop]
    dat[[np]]=aa
  }
  dat2=dat
  mdat=melt(dat1)
  mdat=na.omit(mdat)
  mdat[,2]=abs(mdat[,2])
  mdat=mdat[-which(mdat$variable=="TO"),]####delete#####nops
  pic1=list()
  
  
 ###########################
  #############################
  ###########################
# --- Construct a Data Frame ---
  res=res1=res2 <- data.frame(matrix(nrow = 48, ncol = 8))
  
  for (i in seq_along(dat1)) {
   
    mean_dat1 <- colMeans(dat1[[i]], na.rm = TRUE)
    mean_dat2 <- colMeans(dat2[[i]], na.rm = TRUE)
    res1[i, ] <- round(mean_dat1, 5)
    res2[i, ] <- round(mean_dat2, 5)
  }
  colnames(res1) <- colnames(dat1[[1]])
  colnames(res2) <- colnames(dat1[[1]])
for (i in 1:12) {#simulation
  for (i in 1:4) { #ceinaro "n=200,p=25","n=400,p=50", "n=400,p=100", "n=1000,p=150"
  }
}  
  scenarios <- c("n=200,p=25", "n=400,p=50", "n=400,p=100", "n=1000,p=150")
  

  row_names <- c()

  for (i in 1:12) {

    for (j in seq_along(scenarios)) {
      row_names <- c(row_names, paste0("Sim", i, "_", scenarios[j]))
    }
  }
  

  rownames(res1) <- row_names
  rownames(res2) <- row_names
  head(res)
  tail(res) 
  
  tab1=res1
  tab2=res2
  # 查看当前工作目录
  getwd()
  
  # 将 tab1 写出到当前工作目录下的 tab1.csv

# --- Save Results to CSV ---
  write.csv(tab1, file = "tab1.csv", row.names = T)
  
  # 将 tab2 写出到当前工作目录下的 tab2.csv

# --- Save Results to CSV ---
  write.csv(tab2, file = "tab2.csv", row.names = T)
  
  # 读取 CSV 文件

# --- Load CSV Data ---
  data <- read.csv("tab2.csv", row.names = 1)
  
  tb1=matrix()
  for (i in 1:nrow( data)) {
    row=  data[i,]
    tb1[i]=paste(names(latex_rows)[i],
          paste(row, collapse = " & "),
          sep = " & ")
    tb1[i]=paste0(tb1[i]," \\\\")
  } 
    
  # 修正的行拼接逻辑
  writeLines( tb1, "out.tex_saved2.txt")
  # 添加行尾和换行符
  latex_rows <- paste0(latex_rows, " \\\\")
  
  # 准备表头和表尾
  table_header <- c(
    "\\begin{tabular}{l", 
    paste(rep("l", ncol(data)), collapse = ""),
    "}\n\\toprule\n",
    paste("Row Name & ", paste(colnames(data), collapse = " & "), " \\\\\n", sep = ""),
    "\\midrule\n"
  )
  
  table_footer <- "\\bottomrule\n\\end{tabular}\n"
  
  # 合并所有部分
  latex_table <- c(table_header, latex_rows, table_footer)
  
  # 保存为文本文件
  output_file <- "table_output.txt"
  writeLines(latex_table, con = output_file)
  
  # 打印保存路径
  cat("LaTeX table saved to:", output_file, "\n")
  