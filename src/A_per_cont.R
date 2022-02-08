setwd("~/HGTnew/")
library(base)
# library(foreach)
library(stringr)
library(seqinr)
library(plotrix)
library(ComplexHeatmap)
library(ggpubr)


args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1 | args[1] == '--help'){
  message("usage: Rscript calculate_A_from_mummer_files.R species1 species2")
  stop()
}

species1=args[1]
species2=args[2]


comm=paste("mkdir -p ~/HGTnew/plots/perContinent/",species1,"_",species2,sep="")
system(comm)
country_per_cont=read.table("~/HGTnew/data/processed/All_countries_file.txt",head=T)
country_per_cont=country_per_cont[order(country_per_cont$Continent),]

Prefactor <- data.frame()
files1 <- Sys.glob(file.path("~/HGTnew","data","processed",species1,paste0("*",species1,"_*.fasta")))
files2 <- Sys.glob(file.path("~/HGTnew","data","processed",species2,paste0("*",species2,"_*.fasta")))


species1_dir=paste("~/HGTnew/data/processed/",species1,sep="")
species2_dir=paste("~/HGTnew/data/processed/",species2,sep="")

system(paste0("mkdir -p ~/HGTnew/plots/",species1,"_",species2))

# files <- list.files(path=paste0("~/HGTnew/data/processed/",species1,"/"),
#                     pattern=paste("*.fasta$"), full.names = TRUE)
# foreach (file = files) %do%
# {
#   system(paste0("seqkit seq -m 40000 ",file," > ",file,".4e4.fasta"))
# }


toRm=c()
for (i in 1:length(files2)) {
  print(i)
  size=file.info(files2[i])$size
  if (size<101){
    toRm=c(toRm,i)
  }
}
if (length(toRm)>0){
  files2=files2[-toRm]
}
toRm=c()
for (i in 1:length(files1)) {
  print(i)
  size=file.info(files1[i])$size
  if (size<101){
    toRm=c(toRm,i)
  }
}
if (length(toRm)>0){
  files1=files1[-toRm]
}
toRm=c()

#temp file to store combination of mummer outputs
tempMum="~/HGTnew/data/tmp/TempMum"
tempHist="~/HGTnew/data/tmp/TempHist"

AllContinents=unique(country_per_cont$Continent)
AllContinents=AllContinents[c(6,1,2,3,4,5,7)]

if (file.exists(paste0("~/HGTnew/data/processed/results/",species1,"_",species2,"_perCont.RData"))==F){
  
  for (ThisContI in AllContinents){
    LiTot=0
    AllCountriesI=subset(country_per_cont,Continent == ThisContI)$Country
    for (country_i in AllCountriesI){
      spfile1 <- Sys.glob(file.path("~/HGTnew","data","processed",species1,paste0(species1,"_",country_i,".fasta")))
      if(length(spfile1)>0 && file.info(spfile1)$size > 5){
        Li <- read.table(paste0(spfile1,".L"))$V1 
        LiTot=LiTot+Li
      }
    }  
    for (ThisContJ in AllContinents){
      LjTot=0
      #reinitialize temp mum file
      comm=paste("rm",tempMum,"; touch",tempMum)
      system(comm)
      comm=paste("rm",tempHist,"; touch",tempHist)
      system(comm)
      AllCountriesJ=subset(country_per_cont,Continent == ThisContJ)$Country
      #Measure total length Species J
      for (country_j in AllCountriesJ){
        spfile2 <- Sys.glob(file.path("~/HGTnew","data","processed",species2,paste0(species2,"_",country_j,".fasta")))
        if(length(spfile2)>0 && file.info(spfile2)$size > 5){
          Lj <- read.table(paste0(spfile2,".L"))$V1 
          LjTot=LjTot+Lj
        }
      }
      
      for (country_i in AllCountriesI){
        for (country_j in AllCountriesJ){
        spfile2 <- Sys.glob(file.path("~/HGTnew","data","processed",species2,paste0(species2,"_",country_j,".fasta")))
        if(length(spfile2)>0 && file.info(spfile2)$size > 5){
          sample_j <- paste0(species2,"_",country_j)
          sample_i <- paste0(species1,"_",country_i)
          print(paste(sample_i,sample_j))
          filename <- Sys.glob(file.path("~/HGTnew/data/processed/mummer*",paste0(species1,"_",species2),
                                     paste0(species1,'-',country_i,"_",species2,'-',country_j,".h")))
          print(filename)
      
          if (length(filename)>0 && file.info(filename)$size > 0){
           comm=paste0("sed '/^>/ d' ",filename," |sed 's/.* //'   >>",tempMum)
           system(comm)
          }
        }
        }
      }
      comm2=paste("sort -n",tempMum,"|uniq -c >",tempHist)
      system(comm2)
      size=file.info(tempHist)$size
      if (size>1){
      L=read.table(tempHist)
      L <- L[!is.na(as.numeric(L$V2)),]
      L$V2 <- as.numeric(L$V2)
      L$V1 <- as.numeric(L$V1) 
      if (nrow(L)>5 & max(L$V2)>100){
        p <- weighted.hist(x=log(L$V2),w=L$V1,breaks=10,plot=FALSE)
        r <- exp(p$mids)
        m <- p$counts/diff(exp(p$breaks))/LiTot/LjTot
        A <- 100*sum(L$V1[L$V2>=100]*L$V2[L$V2>=100])/LiTot/LjTot
        
        pdf(paste0("~/HGTnew/plots/perContinent/",species1,"_",species2,"/",ThisContI,"_",ThisContJ,".pdf"))
        plot(log10(r),log10(m),col="black",pch=3)
        lines(log10(r),log10(A/r^3))
        title(paste0(ThisContJ," ",A))
        dev.off()
        Prefactor <- rbind(Prefactor,data.frame(Cont1=ThisContI,
                                                Cont2=ThisContJ,prefactor=A))
        }
      }
    }
  }
  
  PrefactorMatrix <- matrix(NA,nrow=length(AllContinents),ncol=length(AllContinents))
  PrefactorReciprocal <- data.frame()
  colnames(PrefactorMatrix) <- AllContinents
  rownames(PrefactorMatrix) <- AllContinents
  for (conti1 in AllContinents)
  {
    for (conti2 in AllContinents)
    {
      Ind <- which(Prefactor$Cont1==conti1 & Prefactor$Cont2==conti2)
      if (length(Ind)>0)
      {
        PrefactorMatrix[conti1,conti2] = Prefactor$prefactor[Ind]
      }
    }
  }
  if (!is.na(sd(PrefactorMatrix[!is.na(PrefactorMatrix)]))){
    pdf(paste0("~/HGTnew/plots/perContinent/",species1,"_",species2,"/heatmap.pdf"))
    p <- Heatmap(log10(PrefactorMatrix[rowSums(is.na(PrefactorMatrix))!=ncol(PrefactorMatrix),
                                       colSums(is.na(PrefactorMatrix))!=nrow(PrefactorMatrix)]),    
                 column_names_gp = grid::gpar(fontsize = 5),
                 row_names_gp = grid::gpar(fontsize = 5),
                 show_row_names = T,
                 show_column_names = T,
                 cluster_columns = F, 
                 cluster_rows = F, 
                 row_title=species1,
                 column_title=species2,
                 na_col = "black")
    print(p)
    dev.off()
    # }
    # }
    # }
    # 
    # }
  
  }else{
      print("err won't do heatmap")
    }
  
  Prefactor <- data.frame()
  for (continent1 in rownames(PrefactorMatrix))
  {
    for (continent2 in colnames(PrefactorMatrix))
    {
      if (!is.na(PrefactorMatrix[continent1,continent2]))
      {
        Prefactor <- rbind(Prefactor,data.frame(species1=continent1,species2=continent2,A=log10(PrefactorMatrix[continent1,continent2]),Same=(continent1==continent2)))
      }
    }
  }
  save.image(paste0("~/HGTnew/data/processed/results/",species1,"_",species2,"_perCont.RData"))

}else{
  load(paste0("~/HGTnew/data/processed/results/",species1,"_",species2,"_perCont.RData"))
  }


print("first")
print(dim(Prefactor))

if (dim(Prefactor)[1]>2){
  pdf(paste0("~/HGTnew/plots/perContinent/",species1,"_",species2,"/SameDiffCountries.pdf"))
  p <- ggboxplot(Prefactor, x = "Same", y = "A",
                  color = "Same", palette = "jco")+
     geom_jitter(cex=0.5,aes(col=Same))
  p= p  + stat_compare_means(method = "wilcox.test")
  print(p)
  dev.off()

  print("here")

  Prefactor$species1=as.character(Prefactor$species1)
  Prefactor$species2=as.character(Prefactor$species2)
   
   SameCountriesRatios <- c()
   DiffCountriesRatios <- c()
   for (i in 1:nrow(Prefactor))
   {
     Ind <- which(Prefactor$species1==Prefactor$species2[i] & Prefactor$species2==Prefactor$species1[i])
     DiffCountriesRatios <- c(DiffCountriesRatios,Prefactor$A[i]-Prefactor$A[-c(Ind,i)])
     if (length(Ind)==1 & Prefactor$species1[i]!=Prefactor$species2[i])
     {
       SameCountriesRatios <- c(SameCountriesRatios,Prefactor$A[i]-Prefactor$A[Ind])
     }
   }

  print("there")
  if (length(DiffCountriesRatios) >10000){
    print("supp")
    print(dim(DiffCountriesRatios))
    print(length(DiffCountriesRatios))
    DiffCountriesRatios <- sample(DiffCountriesRatios,10000,replace=FALSE)
  }
  print("right here")

   Ratios <- data.frame(ratios=abs(c(SameCountriesRatios,DiffCountriesRatios)),Same=c(rep("same",length(SameCountriesRatios)),rep("diff",length(DiffCountriesRatios))))
#   
#   
   pdf(paste0("~/HGTnew/plots/perContinent/",species1,"_",species2,"/Ratios.pdf"))
   p <- ggviolin(Ratios, x = "Same", y = "ratios",
                 draw_quantiles = 0.5,
                 color = "Same", palette = "jco")+geom_jitter(cex=0.5,aes(col=Same))
   p = p + stat_compare_means(method = "wilcox.test")
   print(p)
   dev.off()
   
   # Test Stat
   Prefactor2=Prefactor
   for (i in AllContinents){
   Prefactor2[i]=0
   Prefactor2[which(Prefactor2$species1 == i | Prefactor2$species2 == i ),i]=1
   }
   
   
   lm1=lm(data = Prefactor2,log(-A)~NorthAmerica + Africa + Asia + Australia + Europa + Antartica + Same)
   summary(lm1)
   write.table(file =  paste0("~/HGTnew/plots/perContinent/",species1,"_",species2,"/test_stat.txt"),
                x=summary(lm1)$coefficients,quote=F,sep='\t')
   # 
   # p=ggplot(Prefactor2,aes(x=species1,fill=species1,y=log(-A)))+
   #   geom_boxplot()
   # print(p)
   # p=ggplot(Prefactor2,aes(x=species2,fill=species2,y=log(-A)))+
   #   geom_boxplot()
   # print(p)
  #  Prefactor2$res=aa  
  #  lm2=lm(data = Prefactor2,res~Same)
  # summary(lm2)   
 }


# To exclude North America
toRm=which(Prefactor$species1=="NorthAmerica" | Prefactor$species2=="NorthAmerica")
if (length(toRm>0)){
  PrefactornoNA=Prefactor[-toRm,]
  if (dim(PrefactornoNA)[1]>2){
  pdf(paste0("~/HGTnew/plots/perContinent/",species1,"_",species2,"/SameDiffCountriesnoNA.pdf"))
  p <- ggboxplot(PrefactornoNA, x = "Same", y = "A",
                 color = "Same", palette = "jco")+
    geom_jitter(cex=0.5,aes(col=Same))
  p= p  + stat_compare_means(method = "wilcox.test")
  print(p)
  dev.off()
  
  print("here")
  
  PrefactornoNA$species1=as.character(PrefactornoNA$species1)
  PrefactornoNA$species2=as.character(PrefactornoNA$species2)
  
  SameCountriesRatios <- c()
  DiffCountriesRatios <- c()
  for (i in 1:nrow(PrefactornoNA))
  {
    Ind <- which(PrefactornoNA$species1==PrefactornoNA$species2[i] & PrefactornoNA$species2==PrefactornoNA$species1[i])
  DiffCountriesRatios <- c(DiffCountriesRatios,PrefactornoNA$A[i]-PrefactornoNA$A[-c(Ind,i)])
    if (length(Ind)==1 & PrefactornoNA$species1[i]!=PrefactornoNA$species2[i])
    {
      SameCountriesRatios <- c(SameCountriesRatios,PrefactornoNA$A[i]-PrefactornoNA$A[Ind])
    }
  }
  
  print("there")
  if (length(DiffCountriesRatios) >10000){
    print("supp")
    print(dim(DiffCountriesRatios))
    print(length(DiffCountriesRatios))
    DiffCountriesRatios <- sample(DiffCountriesRatios,10000,replace=FALSE)
  }
  print("right here")
  
  Ratios <- data.frame(ratios=abs(c(SameCountriesRatios,DiffCountriesRatios)),Same=c(rep("same",length(SameCountriesRatios)),rep("diff",length(DiffCountriesRatios))))
  #   
  #   
  pdf(paste0("~/HGTnew/plots/perContinent/",species1,"_",species2,"/RatiosnoNA.pdf"))
  p <- ggviolin(Ratios, x = "Same", y = "ratios",
                draw_quantiles = 0.5,
                color = "Same", palette = "jco")+geom_jitter(cex=0.5,aes(col=Same))
  p = p + stat_compare_means(method = "wilcox.test")
  print(p)
  dev.off()
}
}




 
