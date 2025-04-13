################################################################################
# Module: Trend Plot
# Description: This script is part of the SCBM framework and implements
# core functionality related to trend plot.
################################################################################


# ====================
# Load Required Libraries
# ====================
library(BART)

# ====================
# Load Required Libraries
# ====================
library(glmnet)

# ====================
# Load Required Libraries
# ====================
library(grpreg)

# ====================
# Load Required Libraries
# ====================
library(parallel)

# ====================
# Load Required Libraries
# ====================
library(grf)

# ====================
# Load Required Libraries
# ====================
library(foreach)

# ====================
# Load Required Libraries
# ====================
library(ggplot2)

# ====================
# Load Required Libraries
# ====================
library(dplyr)

# ====================
# Load Required Libraries
# ====================
library(gridExtra)

# ====================
# Load Required Libraries
# ====================
library(doParallel)
set.seed(99)
#############data input################
data(ACTG175)
ACTG=ACTG175
new_data <- ACTG[, c("pidnum","arms","age", "wtkg", "hemo", "homo", "drugs",  "karnof" , "preanti", "race", "gender", "symptom", "cd40", "cd80", "cd420")]
new_data <- subset(new_data, arms == 0 | arms == 1)
group=new_data[,2]###
g1 <- sort(unique(group))
gn <- length(group)
g<-numeric(gn)
g[which(group==g1[1])]=0
g[which(group==g1[2])]=1
#############data process############
X=new_data[,c(-1,-2,-ncol(new_data))]
cd40=new_data$cd40
cd420=new_data$cd420
quan=function(x) {
  ul=quantile(x, c(0.05,0.95))
  low=which(x < ul[1])
  up=which(x > ul[2])
  x[low]=ul[1]
  x[up]=ul[2]
# --- Return Final Output ---
  return(x)}
cd40=quan(cd40)
cd420=quan(cd420)
ytt=(cd420-cd40)/cd40
# --- Construct a Data Frame ---
datfortrain=data.frame(ytt,g,X)
ps=0.5
Zpe=g*ytt/ps+(1-g)*(-ytt)/(1-ps)
y_ps=Zpe
psFROEST=regression_forest(X,g)
nboots=80
maxn=ncol(X)
per_resamp=0.4
indat=datfortrain
bcd=list()
deg=1
#######################func#################

# --------------------------------------------------
# Function: earthlasso
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
earthlasso <- function(datain,TO,btype=1,deg,maxn,samp){
  
  ############start####################

# ====================
# Load Required Libraries
# ====================
  library("earth")

# ====================
# Load Required Libraries
# ====================
  library(glmnet)
  DAT=sample(nrow(datain),round(samp*nrow(datain)))
  dat=datain[DAT,]
  oobdat=datain#[-DAT,]
  TTO=TO[DAT]
  oobTO=TO
  X=dat[,c(-1,-2)]
  oobX=oobdat[,c(-1,-2)]
  ytt=dat[,1]
  oobytt=oobdat[,1]
  Z=dat[,2]
  group=oobdat[,2]###??????
# --- Construct a Data Frame ---
  dat0=data.frame(Y=TTO,X)# sampled data
  rawtauar=numeric()
  rss3=numeric()
  for (d in deg) {
    q0=earth(Y~.,data=dat0,degree=d,thresh=0.001,nk=maxn,pmethod="none")
    bx=q0$bx[,-1]
    assign(paste0("bx",d),bx)
    dir=q0$dirs[-1,]
    assign(paste0("dir",d),dir)
    cut=q0$cuts[-1,]
    assign(paste0("cut",d),cut)
    coe=q0$coefficients[-1,]
    assign(paste0("coe",d),coe)
    XX_test=numeric()
    py1=predict(q0,newdata = oobX)
    rss=(py1-oobTO)^2
    #=====================================
    rss3[d]=sum(rss)
    assign(paste0("q0",d),q0)
  }
  df=which(rss3==min(na.omit(rss3)))[1]
  bx=get(paste0("bx",df))
  dir=get(paste0("dir",df))
  cut=get(paste0("cut",df))
  coe=get(paste0("coe",df))
  q0=get(paste0("q0",df))
  bcd=list(bx,cut,dir,coe,q0)
# --- Return Final Output ---
  return(bcd)
}
##################################
####################################
####################################
###################loop####################
vimp=numeric(ncol(X))
cl <- makeCluster(25)  
registerDoParallel(cl)

# Parallel execution with foreach
bcd <- foreach(i = 1:nboots) %dopar% {
  earthlasso(indat, Zpe, deg = deg, maxn = maxn, samp = per_resamp)
}

# Stop the cluster
stopCluster(cl)
cut0=numeric()
dir0=numeric()
coe0=numeric()    
q0=list()
for (i in 1:nboots) {
  cut0=rbind(cut0,bcd[[i]][[2]])
  dir0=rbind(dir0,bcd[[i]][[3]])
# --- Construct a Data Frame ---
  coe0=rbind(coe0,data.frame(bcd[[i]][[4]]))
  q0[[i]]=bcd[[i]][[5]]
}

#============group lasso==============#

XX_test=numeric()
for (i in 1:nrow(X)) {
  xin=matrix(unlist(rep(X[i,],nrow(cut0))),nrow = nrow(cut0),byrow = TRUE)-cut0
  dirin=xin*dir0
  a02=which(dir0==2)##
  dirin[which(dirin==0)]=1
  dirin[dirin<0]=0
  dirin[a02]=xin[a02]
# --- Construct a Data Frame ---
  bs=data.frame(apply(dirin,1,prod))#*coe0#*coe0#
  bs=unlist(bs)
  XX_test=rbind(XX_test,bs)
}
xxg=cbind(XX_test*g,XX_test*(1-g))
groupx=c(c(1:ncol(XX_test)),c(1:ncol(XX_test)))
cv.glmdl1=cv.grpreg(xxg,ytt,group=groupx,nfolds=10)
fittau0_err=predict(cv.glmdl1,X=xxg)
fittau_err=abs(ytt-fittau0_err)

trendat=cbind(ytt,fittau0_err)
# --- Construct a Data Frame ---
sorted_trendat <- data.frame(trendat[order(trendat[, 1]), ])
n <- nrow(sorted_trendat)
quartile_size <- ceiling(n / 4)  # Calculate size for each quartile (round up)

# Split the sorted data into 4 parts
quartiles <- split(sorted_trendat, rep(1:4, each = quartile_size, length.out = n))

quartile_means <- lapply(quartiles, function(df) colMeans(df[, 1:2])) 

# Combine means into a data frame for plotting
quartile_means_df <- do.call(rbind, quartile_means)

# Normalize the matrix by columns
quartile_means_df1 <- apply(quartile_means_df, 2, function(col) (col - min(col)) / (max(col) - min(col)))

# --- Construct a Data Frame ---
quartile_means_df <- data.frame(Quartile = paste0("subgroup", 1:4), quartile_means_df1)

# Reshape data for plotting

# ====================
# Load Required Libraries
# ====================
library(reshape2)
quartile_means_melt <- melt(quartile_means_df, id.vars = "Quartile", variable.name = "Variable", value.name = "Mean")

# Plot line chart for each variable with colors and legend

# ====================
# Load Required Libraries
# ====================
library(ggplot2)

# Correct line plot for the provided data

# ====================
# Load Required Libraries
# ====================
library(ggplot2)

pic=ggplot(quartile_means_melt, aes(x = Quartile, y = Mean, group = Variable, color = Variable)) +
  geom_line(size = 1.8) +  # Add lines
  geom_point(size = 4) +   # Add points
  scale_color_manual(
    values = c("ytt" = "#E69F00", "fittau0_err" = "#0B60B0"),  # Low-saturation yellow and blue
    labels = c("ytt" = "True ATE", "fittau0_err" = "estimated HTE")  # Adjust labels for legend
  ) +
  labs(
    title = "Trend Plot of ATE in subgroups/Mean of HTE in subgroups",
    x = "subgroups",
    y = "ATE/Mean of HTE(normalized)",
    color = "Variable"  # Legend title
  ) +
  theme_minimal()+theme(
    legend.position = "bottom",
    text = element_text(size = 18),
    plot.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  ) +guides(fill = guide_legend(nrow = 1))




# --- Save Plot to File ---
ggsave(
  filename = paste0(paste0(getwd(), "/pic/"), "tendpic", i, ".eps"),
  plot = pic,
  width = 13,
  height = 8,
  device = cairo_ps  # 使用 Cairo PostScript 设备生成 EPS 文件
)