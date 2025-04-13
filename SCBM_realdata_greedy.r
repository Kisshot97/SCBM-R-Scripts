library(parallel)
library(BART)
set.seed(10086)
om=function(sam,i){
  library(doParallel)
  library(foreach)
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
     return(x)}
  ytt=(cd420-cd40)/cd40
  ytt=quan(ytt)
  datalist0=data.frame(0,g,0,X,ytt)
  datalist1=datalist0[-sam,]
  datalist2=datalist0[sam,]
  


  source("E:/SCBM/test326/fun_GL_GREED.R")
  #hte_test1=marslasso_GLPSTO(traindat,testdat,ps,deg=1,boost=10,per_resamp=0.2)

  # write.csv(datalist1,sprintf("E:/SCBM/dat.24.05.39/train/scenario_sam%s_boost%s.csv",1,1),row.names = FALSE)
  # write.csv(datalist2,sprintf("E:/SCBM/dat.24.05.39/test/scenario_sam%s_boost%s.csv",1,1),row.names = FALSE)
  # datalist1=read.csv(sprintf("E:/SCBM/dat.24.05.39/train/scenario_sam%s_boost%s.csv",1,1))
  # datalist2=read.csv(sprintf("E:/SCBM/dat.24.05.39/test/scenario_sam%s_boost%s.csv",1,1))
  traindat=cbind(1,datalist1)
  testdat=cbind(1,datalist2)
  ps=rep(0.5,nrow(datalist1))
 # hte_test1=marslasso_GLPSTO(traindat,testdat,ps,deg=1,boost=10,per_resamp=0.2)
  # 设置并行核心数，根据你的机器调整
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # 创建参数网格
  g_boost_values <- seq(20, 200, 20)
  g_sam_values <- seq(0.2, 0.8, 0.1)
  
  # 使用 foreach 替代嵌套 for 循环
  foreach(g_boost = g_boost_values) %:%
    foreach(g_sam = g_sam_values, .combine = 'c') %dopar% {
      source("E:/SCBM/test326/fun_GL_GREED.R")
      # 创建结果数据框
      res.save <- data.frame(matrix(0, nrow(datalist2), 5))
      colnames(res.save) <- c("true", "hte1", "hte2", "hte3", "hte4")
      res.save$true <- testdat$ytt
      
      # 调用函数
      hte_test1 <- marslasso_GLPSTO(traindat, testdat, ps, deg = 1, boost = g_boost, per_resamp = g_sam)
      hte_test2 <- marslasso_GLPSTO(traindat, testdat, ps, deg = 2, boost = g_boost, per_resamp = g_sam)
      hte_test3 <- marslasso_GLPSTO(traindat, testdat, ps, deg = 3, boost = g_boost, per_resamp = g_sam)
      hte_test4 <- marslasso_GLPSTO(traindat, testdat, ps, deg = c(1:3), boost = g_boost, per_resamp = g_sam)
      
      # 填充结果
      res.save$hte1 <- hte_test1
      res.save$hte2 <- hte_test2
      res.save$hte3 <- hte_test3
      res.save$hte4 <- hte_test4
      
      # 保存结果
      output_path <- sprintf("E:/SCBM/dat.24.05.39/out/scenario%s_sam%s_boost%s.csv", i, g_sam, g_boost)
      write.csv(res.save, output_path)
 
      print(sprintf("Finished boost=%s, sam=%s", g_boost, g_sam))
    }
  
  # 关闭并行计算
  stopCluster(cl)

}




#####################
####################
#####################

  library(rpart)
  library(Rcpp)
  library(grf)
  library(earth)
  library(BART)
  library(bartCause)
  library(randomForest)

#==========================
data(ACTG175)
ACTG=ACTG175
new_data <- ACTG[, c("pidnum","arms","age", "wtkg", "hemo", "homo", "drugs",  "karnof" , "preanti", "race", "gender", "symptom", "cd40", "cd80", "cd420")]
new_data <- subset(new_data, arms == 0 | arms == 1)
ndat=nrow(new_data)
split_size <- ceiling(ndat / 10)

index0= sample(1:ndat)

for (i in 1:10) {
  start_idx <- (i - 1) * split_size + 1
  end_idx <- min(i * split_size, ndat)
  index_list <- index0[start_idx:end_idx]
  om(index_list,i)
  
}

sam=index_list





