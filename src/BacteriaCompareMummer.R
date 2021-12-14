setwd("~/HGTnew/")
files <- list.files(path="~/HGTnew/data/processed/EscherichiaColi/",pattern=paste("*.fasta$"), full.names = TRUE)
foreach (file = files) %do%
{
  system(paste0("seqkit seq -m 40000 ",file," > ",file,".4e4.fasta"))
}

sed -i '/^$/d' ~/HGTnew/data/processed/EscherichiaColi/*.fasta
sed -i '/^$/d' ~/HGTnew/data/processed/KlebsiellaPneumoniae/*.fasta


# export OPENBLAS_NUM_THREADS=1
# export GOTO_NUM_THREADS=1
# export OMP_NUM_THREADS=1
setwd("~/HGTnew/")
library(base)
library(foreach)
require(doParallel)
registerDoParallel(30)
library(stringr)
library(seqinr)
setwd("~/HGTnew/")
write("", paste0("~/HGTnew/data/processed/Command.txt"),append=FALSE)
files1 <- Sys.glob(file.path("~/HGTnew","data","processed","EscherichiaColi","*EscherichiaColi_*.fasta"))
files2 <- Sys.glob(file.path("~/HGTnew","data","processed","KlebsiellaPneumoniae","*KlebsiellaPneumoniae_*.fasta"))
foreach (i = 1:length(files1)) %dopar%
{
  sample_i <- strsplit(files1[i],"/")[[1]];sample_i <- sample_i[length(sample_i)];sample_i <- str_replace(sample_i,".fasta","")
  print(sample_i)
  foreach (j = 1:length(files2)) %do%
  {
    sample_j <- strsplit(files2[j],"/")[[1]];sample_j <- sample_j[length(sample_j)];sample_j <- str_replace(sample_j,".fasta","")
    Command <- paste0("mummer -maxmatch -b -n -l 300 -qthreads 1 -threads 1 -F ",
                      files1[i]," ",files2[j]," > ~/HGTnew/data/processed/mummer/EscherichiaColi_KlebsiellaPneumoniae/",sample_i,"_",sample_j,".mumm")
    system(Command)
  }
}

setwd("~/HGTnew/")
files <- Sys.glob(file.path("~/HGTnew/data/processed/mummer/EscherichiaColi_KlebsiellaPneumoniae","*.mumm"))
for (i in 1:length(files)) 
{
  if (!file.exists(paste0(files[i],".h")))
  {
    file <- files[i]
    print(file)
    Command <- paste0("sed '/^>/ d' ",file," > ",file,".r")
    system(Command)
    Command <- paste0("sed -i 's/.* //' ",file,".r")
    system(Command)
    Command <- paste0("sort -n ",file,".r | uniq -c > ",file,".h")
    system(Command)
  }
}


setwd("~/HGTnew/")
files1 <- Sys.glob(file.path(".","data","processed","EscherichiaColi","*EscherichiaColi_*.fasta"))
files2 <- Sys.glob(file.path("~/HGTnew","data","processed","KlebsiellaPneumoniae","*KlebsiellaPneumoniae_*.fasta"))
for (i in 1:length(files2)) 
{
  print(i)
  fi <- read.fasta(file = files2[i])
  Li <- sum(getLength(fi))
  write(Li, paste0(files2[i],".L"),append=FALSE)
}






setwd("~/HGTnew/")
library(base)
library(foreach)
require(doParallel)
registerDoParallel(20)
library(stringr)
library(seqinr)
library(plotrix)
library(ComplexHeatmap)
Prefactor <- data.frame()
files1 <- Sys.glob(file.path(".","data","processed","EscherichiaColi","*EscherichiaColi_*.fasta"))
files2 <- Sys.glob(file.path("~/HGTnew/data/processed/KlebsiellaPneumoniae/*KlebsiellaPneumoniae_*.fasta"))
for (i in 1:length(files1)) 
{
  Li <- read.table(paste0(files1[i],".L"))$V1
  sample_i <- strsplit(files1[i],"/")[[1]];sample_i <- sample_i[length(sample_i)];sample_i <- str_replace(sample_i,".fasta","")
  country_i <- strsplit(sample_i,"_")[[1]][2]
  print(sample_i)
  for (j in 1:length(files2)) 
  {
    sample_j <- strsplit(files2[j],"/")[[1]];sample_j <- sample_j[length(sample_j)];sample_j <- str_replace(sample_j,".fasta","")
    country_j <- strsplit(sample_j,"_")[[1]][2]
    filename <- paste0("./data/processed/mummer/EscherichiaColi_KlebsiellaPneumoniae/",sample_i,"_",sample_j,".mumm.h")
    if (file.exists(filename) & file.info(filename)$size != 0)
    {
      print(paste0(sample_i," ",sample_j))
      L <- read.table(paste0("./data/processed/mummer/EscherichiaColi_KlebsiellaPneumoniae/",sample_i,"_",sample_j,".mumm.h"))
      L <- L[!is.na(as.numeric(L$V2)),]; L$V2 <- as.numeric(L$V2);  L$V1 <- as.numeric(L$V1); 
      
      if (nrow(L)>5 & max(L$V2)>300)
      {
        Lj <- read.table(paste0(files2[j],".L"))$V1
        p <- weighted.hist(x=log(L$V2),w=L$V1,breaks=10,plot=FALSE)
        r <- exp(p$mids)
        m <- p$counts/diff(exp(p$breaks))/Li/Lj
        A <- 300*sum(L$V1[L$V2>=300]*L$V2[L$V2>=300])/Li/Lj
        pdf(paste0("./plots/EscherichiaColi_KlebsiellaPneumoniae/",sample_i,"_",sample_j,".pdf"))
        # plot(r,log10(m));
        plot(log10(r),log10(m),col="black",pch=3);
        lines(log10(r),log10(A/r^3))
        title(paste0(sample_j," ",A))
        dev.off()
        Prefactor <- rbind(Prefactor,data.frame(EscherichiaColi=country_i,KlebsiellaPneumoniae=country_j,prefactor=A))
        countries <- unique(c(Prefactor$EscherichiaColi,Prefactor$KlebsiellaPneumoniae))
        PrefactorMatrix <- matrix(NA,nrow=length(countries),ncol=length(countries))
        PrefactorReciprocal <- data.frame()
        colnames(PrefactorMatrix) <- countries
        rownames(PrefactorMatrix) <- countries
        for (country1 in countries)
        {
          for (country2 in countries)
          {
            Ind <- which(Prefactor$EscherichiaColi==country1 & Prefactor$KlebsiellaPneumoniae==country2)
            if (length(Ind)>0)
            {
              PrefactorMatrix[country1,country2] <- Prefactor$prefactor[Ind]
            }
          }
        }
        if (!is.na(sd(PrefactorMatrix[!is.na(PrefactorMatrix)])))
        {
          pdf(paste0("./plots/EscherichiaColi_KlebsiellaPneumoniae/heatmap.pdf"))
          p <- Heatmap(log10(PrefactorMatrix[rowSums(is.na(PrefactorMatrix))!=ncol(PrefactorMatrix),colSums(is.na(PrefactorMatrix))!=nrow(PrefactorMatrix)]),    
                       column_names_gp = grid::gpar(fontsize = 5),
                       row_names_gp = grid::gpar(fontsize = 5),
                       show_row_names = T,
                       show_column_names = T,
                       cluster_columns = F, 
                       cluster_rows = F, 
                       row_title="Escherichia coli",
                       column_title="Klebsiella pneumoniae",
                       na_col = "black")
          print(p)
          dev.off()
        }
      }
    }
  }
}

Prefactor <- data.frame()
for (country1 in rownames(PrefactorMatrix))
{
  for (country2 in colnames(PrefactorMatrix))
  {
    if (!is.na(PrefactorMatrix[country1,country2]))
    {
      Prefactor <- rbind(Prefactor,data.frame(EscherichiaColi=country1,KlebsiellaPneumoniae=country2,A=log10(PrefactorMatrix[country1,country2]),Same=(country1==country2)))
    }
  }
}
library(ggpubr)
pdf(paste0("./plots/EscherichiaColi_KlebsiellaPneumoniae/SameDiffCountries.pdf"))
p <- ggboxplot(Prefactor, x = "Same", y = "A",
               color = "Same", palette = "jco",
               add = "jitter")
p + stat_compare_means(method = "wilcox.test")
dev.off()

SameCountriesRatios <- c()
DiffCountriesRatios <- c()
for (i in 1:nrow(Prefactor))
{
  Ind <- which(Prefactor$EscherichiaColi==Prefactor$KlebsiellaPneumoniae[i] & Prefactor$KlebsiellaPneumoniae==Prefactor$EscherichiaColi[i])
  DiffCountriesRatios <- c(DiffCountriesRatios,Prefactor$A[i]-Prefactor$A[-c(Ind,i)])
  if (length(Ind)==1 & Prefactor$EscherichiaColi[i]!=Prefactor$KlebsiellaPneumoniae[i])
  {
    SameCountriesRatios <- c(SameCountriesRatios,Prefactor$A[i]-Prefactor$A[Ind])
  }
}

DiffCountriesRatios <- sample(DiffCountriesRatios,10000,replace=FALSE)

Ratios <- data.frame(ratios=abs(c(SameCountriesRatios,DiffCountriesRatios)),Same=c(rep("same",length(SameCountriesRatios)),rep("diff",length(DiffCountriesRatios))))

library(ggpubr)
pdf(paste0("./plots/EscherichiaColi_KlebsiellaPneumoniae/Ratios.pdf"))
p <- ggviolin(Ratios, x = "Same", y = "ratios",
              draw_quantiles = 0.5,
               color = "Same", palette = "jco",
               add = "jitter")
p + stat_compare_means(method = "wilcox.test")
dev.off()



