# For github
# cd /home/m.sheinman/Development/HGTproject
# git add ./src/*.*
# git commit --all -m "CAT taxonomic annotation for metagenomic analysis"
# git push origin master


setwd("/shared/scratch/m.sheinman/assembly/")
library(base)
library(foreach)
require(doParallel)
registerDoParallel(40)
library(stringr)
library(seqinr)
files <- Sys.glob(file.path(".","*", "metaspades", "*contigs-final.fasta"))
# filter long contigs
foreach (file = files) %do%
{
  system(paste0("/home/m.sheinman/Development/Software/seqkit/seqkit seq -m 10000 ",file," > ",file,".1e4.fasta"))
}


  

# function that outputs data.frame table from mummer file
getDataFrameFromMummerOutPut <- function(mummerfile,Tax1Table,Tax2Table)
{
  Tax1Table$genus[Tax1Table$genus=="Shigella"] <- "Escherichia"
  Tax2Table$genus[Tax2Table$genus=="Shigella"] <- "Escherichia"
  dat <- data.frame()
  seqs <- readLines(mummerfile)
  for (k in 1:length(seqs))
  {
    if (substr(seqs[k],1,1)==">") 
    {
      seq1 <- strsplit(seqs[k]," ")[[1]][2]
      reverse <- ifelse(length(strsplit(seqs[k]," ")[[1]])==3,"reverse","") 
    }
    else if (substr(seqs[k],1,1)==" ") 
    {
      seq2 <- strsplit(seqs[k]," ")[[1]] 
      seq2 <- seq2[seq2!=""]
      x1 <- as.numeric(seq2[2])
      x2 <- as.numeric(seq2[3])
      r <- as.numeric(seq2[4])
      seq2 <- seq2[1]
    }
    else 
    {
      seq <- seqs[k]
      if (seq1!=seq2) 
      {
        Ind1 <- which(Tax1Table$contig==seq2)
        Ind2 <- which(Tax2Table$contig==seq1)
        if (length(Ind1)==1 & length(Ind2)==1)
        {
          L1 <- Tax1Table$contiglength[Ind1]
          L2 <- Tax2Table$contiglength[Ind2]
          if (Tax1Table$genus[Ind1] != Tax2Table$genus[Ind2] & Tax1Table$coverage[Ind1]>0.1e-10 & Tax2Table$coverage[Ind2]>0.1e-10) 
          {rank <- "genus"; taxon1 <- Tax1Table$genus[Ind1]; taxon2 <- Tax2Table$genus[Ind2]; coverage1 <- Tax1Table$coverage[Ind1]; coverage2 <- Tax2Table$coverage[Ind2]}
          else
          {rank <- "none"; taxon1 <- Tax1Table$genus[Ind1]; taxon2 <- Tax2Table$genus[Ind2]; coverage1 <- Tax1Table$coverage[Ind1]; coverage2 <- Tax2Table$coverage[Ind2]}
          datNew <- data.frame(seq1=seq1,x1=x1,L1=L1,taxon1=taxon1,coverage1=coverage1,seq2=seq2,x2=x2,L2=L2,taxon2=taxon2,coverage2=coverage2,r=r,reverse=reverse,seq=seq,rank=rank)
          dat <- rbind(dat,datNew)
        }
      }
    }
  }
  return(dat[dat$rank=="none" & dat$taxon1=="Escherichia" & dat$taxon2=="Escherichia",])
}

# outputs data.frame table from mummer file
files <- Sys.glob(file.path(".","*", "metaspades", "*contigs-final.fasta.4e4.fasta"))
foreach (i = 1:length(files)) %dopar%
{
  sample_i <- strsplit(files[i],"/")[[1]][2]
  print(sample_i)
  Tax1Table <- read.table(paste0(files[i],".blastnRefSeq.classify"), header = TRUE, sep="\t",comment.char = "")
  
  for (j in i:length(files))
  {
    sample_j <- strsplit(files[j],"/")[[1]][2]
    mummerfile <- paste0("./mummer/",sample_i,"_",sample_j,".seq.mumm")
    Tax2Table <- read.table(paste0(files[j],".blastnRefSeq.classify"), header = TRUE, sep="\t",comment.char = "")
    # system(paste0("/home/m.sheinman/Development/Software/MUMmer3.23/mummer -maxmatch -n -s -b -l 100 -F ",files[i]," ",files[j]," > ./mummer/",sample_i,"_",sample_j,".seq.mumm"))
    dat <- getDataFrameFromMummerOutPut(mummerfile, Tax1Table, Tax2Table)
    write.table(dat,
                file = paste0("./mummer/",sample_i,"_",sample_j,".seq.mumm.csv"),
                append=FALSE,row.names=FALSE,col.names=TRUE,sep = "\t",quote=FALSE)
  } 
}

calculateNormalization <- function(Tax1Table, Tax2Table)
{
  Normalization <- foreach (i = 1:nrow(Tax1Table), .combine='+', .init=0) %dopar%
  {
    NormalizationJ <- 0.0
    for (j in 1:nrow(Tax2Table))
    {
      if (Tax1Table$genus[i]==Tax2Table$genus[j] & Tax1Table$coverage[i]>0.1e-10 & Tax2Table$coverage[j]>0.1e-10 & Tax2Table$genus[j]=="Escherichia")
      {
        NormalizationJ <- NormalizationJ + as.double(Tax1Table$contiglength[i])*as.double(Tax2Table$contiglength[j])
      }
    }
    NormalizationJ
  }
  return(Normalization)
}




# make pdfs
files <- Sys.glob(file.path(".","*", "metaspades", "*contigs-final.fasta.4e4.fasta"))[1:12]
for (i in 1:length(files)) 
{
  sample_i <- strsplit(files[i],"/")[[1]][2]
  Tax1Table <- read.table(paste0(files[i],".blastnRefSeq.classify"), header = TRUE, sep="\t",comment.char = "")
  for (j in i:length(files))
  {
    sample_j <- strsplit(files[j],"/")[[1]][2]
    Tax2Table <- read.table(paste0(files[j],".blastnRefSeq.classify"), header = TRUE, sep="\t",comment.char = "")
    if (file.exists(paste0("./mummer/",sample_i,"_",sample_j,".seq.mumm.csv")) & file.size(paste0("./mummer/",sample_i,"_",sample_j,".seq.mumm.csv"))>10)
    {
      ToNorm <- calculateNormalization(Tax1Table, Tax2Table)/ifelse(i==j,2,1)
      dat <- read.table(paste0("./mummer/",sample_i,"_",sample_j,".seq.mumm.csv"), header = TRUE, sep="\t")
      # L <- dat$r[dat$genus1!=dat$genus2 & dat$genus1!="Homo sapiens (taxid 9606)" & dat$genus2!="Homo sapiens (taxid 9606)"]
      L <- dat$r
      if (length(unique(L))>2)
      {
        print(paste0(sample_i," ",sample_j))
        Ltable <- as.data.frame(table(L),stringsAsFactors=FALSE)
        Ltable$L <- as.numeric(Ltable$L)
        if (i==j) {pdf(paste0("./plots/same/",sample_i,"_",sample_j,".pdf"))
        }else {pdf(paste0("./plots/diff/",sample_i,"_",sample_j,".pdf"))}
        n <- 25
        p <- hist(log(L),n,plot=FALSE)
        while(sum(p$counts==0)>100)
        {
          n <- n - 1
          p <- hist(log(L),n,plot=FALSE)
        }
        r <- exp(p$mids)
        m <- p$counts/diff(exp(p$breaks))/ToNorm
        A <- 100*sum(Ltable$L[Ltable$L>=100]*Ltable$Freq[Ltable$L>=100])/ToNorm
        plot(log10(r),log10(m));
        lines(log10(r),log10(A/r^3))
        title(paste0(sample_i," ",A))
        dev.off()
      }
    }
  } 
}





