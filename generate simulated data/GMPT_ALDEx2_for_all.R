#################### ALDEx2#####################################
rm(list=ls())
library(tidyverse)
library(ALDEx2)
set.seed(123)

father_folder = "../Non_diagStable_connectivity/N30_C0.08_Pheno8_20220429T222328/relative_With_Cdiff_exchange_0.1_0.3_0.5_"
sun_foder = dir(father_folder)
for (ii in 1:length(sun_foder)){
  path<- paste0(father_folder,"/",sun_foder[[ii]])
  filenames<-dir(path,pattern = '.txt')
  filepath<-sapply(filenames, function(x){
    paste(path,x,sep='/')})
  data<-lapply(filepath,function(x){
    read.table(x,header=T)})
  
  ########For loop ###############
  dir.create(file.path(path,"ALDEx2_res"))
  for (i in 1:length(filenames)){
    print(i)
    df<-data[[i]]
    mat<-c(rep("G1",(ncol(df)-1)/2),rep("G2",(ncol(df)-1)/2))
    dat<-as.matrix(df[-c(1),-c(1)])
    class(dat) <- "numeric"
    dat <- ceiling(1e6*dat)
    rownames(dat)=c(1:nrow(dat))
    
    Aldex <- aldex.clr(dat, mat, mc.samples=128, denom="all", verbose=F)
    xtt <- aldex.ttest(Aldex, paired.test=FALSE, verbose=FALSE)
    xkw <- aldex.kw(Aldex)
    xeffect <- aldex.effect(Aldex, CI=T, verbose=FALSE)
    Aldex_res <- data.frame(xtt,xkw,xeffect)
    write.table(Aldex_res,paste0(path,"/ALDEx2_res/","ALDEx2_",strsplit(filenames[i],'_')[[1]][3],"_res.txt"),sep="\t",col.names=NA)
  }
}
