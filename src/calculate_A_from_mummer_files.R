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
filter=args[3]

Prefactor <- data.frame()
files1 <- Sys.glob(file.path("~/HGTnew","data","processed",species1,paste0("*",species1,"_*.fasta")))
files2 <- Sys.glob(file.path("~/HGTnew","data","processed",species2,paste0("*",species2,"_*.fasta")))


species1_dir=paste("~/HGTnew/data/processed/",species1,sep="")
species2_dir=paste("~/HGTnew/data/processed/",species2,sep="")

if (filter == 1){
  suffix =  "-filtered-mash"
}else{
  suffix=''
}
system(paste0("mkdir -p ~/HGTnew/plots/",species1,"_",species2,suffix))
 
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

if (file.exists(paste0("~/HGTnew/data/processed/results/",species1,"_",species2,suffix,".RData"))==F){

  
  cc=0
  cc2=0
for (i in 1:length(files1)) {

 if (filter == 1) {  
   AllLs <- read.table(paste0(files1[i],".all.L"))$V1
   filterFile=read.table(paste("/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/data/processed/toFilterMash/",
   species1,"/",species1,"-",country1,"/toFilter.txt"),sep='')
   Li=sum(AllLs$Length[which(AllLs$ID %in% filterFile$V1)])
 }else{
   Li <- read.table(paste0(files1[i],".L"))$V1
}


  sample_i <- strsplit(files1[i],"/")[[1]];sample_i <- sample_i[length(sample_i)];sample_i <- str_replace(sample_i,".fasta","")
  country_i <- strsplit(sample_i,"_")[[1]][2]
  print(sample_i)
  for (j in 1:length(files2)) 
  {
    sample_j <- strsplit(files2[j],"/")[[1]];sample_j <- sample_j[length(sample_j)];sample_j <- str_replace(sample_j,".fasta","")
    country_j <- strsplit(sample_j,"_")[[1]][2]
    
    ###
  if (filter == 1) {
    filename <- Sys.glob(file.path("~/HGTnew/data/processed/mummer*",paste0(species1,"_",species2),
                                   paste0(species1,'-',country_i,"_",species2,'-',country_j,".mum.h-filtered-mash")))
  }else{
        filename <- Sys.glob(file.path("~/HGTnew/data/processed/mummer*",paste0(species1,"_",species2),
                                   paste0(species1,'-',country_i,"_",species2,'-',country_j,".h")))
  }
    print(filename)
    print(paste(sample_i,sample_j))
    if (length(filename)>0 && file.info(filename)$size > 0){
      print(paste0(sample_i," ",sample_j))
      L <- read.table(filename)
      L <- L[!is.na(as.numeric(L$V2)),]; L$V2 <- as.numeric(L$V2);  L$V1 <- as.numeric(L$V1); 
      
      cc2=cc2+1
      if (nrow(L)>5 & max(L$V2)>100)
      {
        cc=cc+1

	if (filter == 1) {
	   AllLs <- read.table(paste0(files1[i],".all.L"))$V1
	   filterFile=read.table(paste("/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/data/processed/toFilterMash/",
	   species1,"/",species1,"-",country1,"/toFilter.txt"),sep='')
	   Lj=sum(AllLs$Length[which(AllLs$ID %in% filterFile$V1)])
	}else{
           Lj <- read.table(paste0(files2[j],".L"))$V1
	}
        p <- weighted.hist(x=log(L$V2),w=L$V1,breaks=20,plot=FALSE)
        r <- exp(p$mids)
        m <- p$counts/diff(exp(p$breaks))/Li/Lj
        A <- 100*sum(L$V1[L$V2>=100]*L$V2[L$V2>=100])/Li/Lj
#       B <- 300*sum(L$V1[L$V2>=300]*L$V2[L$V2>=300])/Li/Lj

        
        pdf(paste0("~/HGTnew/plots/",species1,"_",species2,suffix,"/",sample_i,"_",sample_j,".pdf"))
        
        # plot(r,log10(m));
        plot(log10(r),log10(m),col="black",pch=3);
        lines(log10(r),log10(A/r^3))
#        lines(log10(r),log10(B/r^3),col=2,lty=2)
#	abline(v=log10(300),lty=3)
        title(paste0(sample_j," ",A))
        dev.off()
        Prefactor <- rbind(Prefactor,data.frame(species1=country_i,
                                                species2=country_j,prefactor=A))
      }
    }
  }
}
        
        countries <- unique(c(as.character(Prefactor$species1),as.character(Prefactor$species2)))
        PrefactorMatrix <- matrix(NA,nrow=length(countries),ncol=length(countries))
        PrefactorReciprocal <- data.frame()
        colnames(PrefactorMatrix) <- countries
        rownames(PrefactorMatrix) <- countries
        for (country1 in countries)
        {
          for (country2 in countries)
          {
            Ind <- which(Prefactor$species1==country1 & Prefactor$species2==country2)
            if (length(Ind)>0)
            {
              PrefactorMatrix[country1,country2] = Prefactor$prefactor[Ind]
            }
          }
        }
        if (!is.na(sd(PrefactorMatrix[!is.na(PrefactorMatrix)])))
        {
          pdf(paste0("~/HGTnew/plots/",species1,"_",species2,suffix,"/heatmap.pdf"))
          p <- Heatmap(log10(PrefactorMatrix[rowSums(is.na(PrefactorMatrix))!=ncol(PrefactorMatrix),colSums(is.na(PrefactorMatrix))!=nrow(PrefactorMatrix)]),    
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
}

Prefactor <- data.frame()
for (country1 in rownames(PrefactorMatrix))
{
  for (country2 in colnames(PrefactorMatrix))
  {
    if (!is.na(PrefactorMatrix[country1,country2]))
    {
      Prefactor <- rbind(Prefactor,data.frame(species1=country1,species2=country2,A=log10(PrefactorMatrix[country1,country2]),Same=(country1==country2)))
    }
  }
}
save.image(paste0("~/HGTnew/data/processed/results/",species1,"_",species2,suffix,".RData"))
} else{
  load(paste0("~/HGTnew/data/processed/results/",species1,"_",species2,suffix,".RData"))
}
  print("first")
  print(dim(Prefactor))

 if (dim(Prefactor)[1]>2){
pdf(paste0("~/HGTnew/plots/",species1,"_",species2,suffix,"/SameDiffCountries.pdf"))
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


pdf(paste0("~/HGTnew/plots/",species1,"_",species2,suffix,"/Ratios.pdf"))
p <- ggviolin(Ratios, x = "Same", y = "ratios",
              draw_quantiles = 0.5,
               color = "Same", palette = "jco")+geom_jitter(cex=0.5,aes(col=Same))
p = p + stat_compare_means(method = "wilcox.test")
print(p)
dev.off()
}


