# For github
# cd /home/m.sheinman/Development/HGTproject
# git add ./src/*.*
# git commit --all -m "CAT taxonomic annotation for metagenomic analysis"
# git push origin master


setwd("/shared/scratch/m.sheinman/assembly/")
library(base)
library(foreach)
require(doParallel)
registerDoParallel(20)
library(stringr)
library(seqinr)
files <- Sys.glob(file.path(".","*", "metaspades", "*contigs-final.fasta"))
# filter long contigs
foreach (file = files) %do%
{
  system(paste0("/home/m.sheinman/Development/Software/seqkit/seqkit seq -m 100000 ",file," > ",file,".1e5.fasta"))
}


# copy krake2 database
# downloaded from https://benlangmead.github.io/aws-indexes/k2
# wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20210517.tar.gz
# cp /shared/scratch/m.sheinman/k2_standard/*k2d /dev/shm/
# cp /shared/scratch/m.sheinman/k2_standard/*kmer_distrib /dev/shm/

# getting kraken2 annotation for each contig in .kraken file
files <- Sys.glob(file.path(".","*", "metaspades", "*contigs-final.fasta.1e5.fasta"))
foreach (i = 1:length(files)) %dopar%
{
  print(sample_i)
  sample_i <- strsplit(files[i],"/")[[1]][2]
  Command <- paste0("/home/m.sheinman/Development/Software/kraken2/kraken2 --threads 1 --memory-mapping --use-names --report ",files[i],".kreport --db /dev/shm/ ",files[i]," --output ",files[i],".kraken")
  system(Command)
  Command <- paste0("/home/m.sheinman/Development/Software/Bracken-master/bracken -d /shared/scratch/m.sheinman/k2_standard/ -i " ,files[i], ".kreport -o ", files[i], ".bracken -l G")
  system(Command)
}

# taxonomic clasiffcation with CAT and BAT github.com/dutilh/CAT
# export PATH="/shared/scratch/m.sheinman/CAT-master/CAT_prepare_20210107/Diamond_2.0.6:$PATH"
# export PATH="/shared/scratch/m.sheinman/CAT-master/:$PATH"
setwd("/shared/scratch/m.sheinman/assembly/")
files <- Sys.glob(file.path(".","*", "metaspades", "*contigs-final.fasta.1e5.fasta"))
for (i in 1:length(files))
{
  dir <- dirname(files[i]) 
  Command <- paste0("/shared/scratch/m.sheinman/CAT-master/CAT_pack/CAT contigs -n 40 ",
                    "-c ",files[i],
  " -o ",files[i],".CAT ",
  " -d /shared/scratch/m.sheinman/CAT-master/CAT_prepare_20210107/2021-01-07_CAT_database ",
  " -t /shared/scratch/m.sheinman/CAT-master/CAT_prepare_20210107/2021-01-07_taxonomy/ ",
  " --tmpdir ",dir)
  # system(Command)
  
  Command <- paste0("/shared/scratch/m.sheinman/CAT-master/CAT_pack/CAT add_names --force --only_official ",
  "-i ",files[i],".CAT.contig2classification.txt ",
  "-o ",files[i],".CAT.contig2classification.names.txt ",
  "-t /shared/scratch/m.sheinman/CAT-master/CAT_prepare_20210107/2021-01-07_taxonomy/")
  system(Command)
  
}



  

# function that outputs data.frame table from mummer file
getDataFrameFromMummerOutPut <- function(mummerfile,kraken1,kraken2)
{
  kraken1Table <- read.table(kraken1, header = FALSE, sep="\t")
  kraken2Table <- read.table(kraken2, header = FALSE, sep="\t")
  
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
        Ind <- which(kraken1Table$V2==seq1)
        L1 <- kraken1Table$V4[Ind]
        taxon1 <- kraken1Table$V3[Ind]
        genus1 <- strsplit(taxon1," ")[[1]][1]
        genus1 <- str_replace(genus1,"\\[","");genus1 <- str_replace(genus1,"]","")
        Ind <- which(kraken2Table$V2==seq2)
        L2 <- kraken2Table$V4[Ind]
        taxon2 <- kraken2Table$V3[Ind]
        genus2 <- strsplit(taxon2," ")[[1]][1]
        genus2 <- str_replace(genus2,"\\[","");genus2 <- str_replace(genus2,"]","")
        datNew <- data.frame(seq1=seq1,x1=x1,L1=L1,taxon1=taxon1,genus1=genus1,seq2=seq2,x2=x2,L2=L2,taxon2=taxon2,genus2=genus2,r=r,reverse=reverse,seq=seq)
        dat <- rbind(dat,datNew)
      }
    }
  }
  return(dat)
}

# outputs data.frame table from mummer file
files <- Sys.glob(file.path(".","*", "metaspades", "*contigs-final.fasta.1e5.fasta"))
foreach (i = 1:length(files)) %dopar%
{
  sample_i <- strsplit(files[i],"/")[[1]][2]
  print(sample_i)
  sample_i <- strsplit(files[i],"/")[[1]][2]
  for (j in i:length(files))
  {
    sample_j <- strsplit(files[j],"/")[[1]][2]
    system(paste0("/home/m.sheinman/Development/Software/MUMmer3.23/mummer -maxmatch -n -s -b -l 100 -F ",files[i]," ",files[j]," > ./mummer/",sample_i,"_",sample_j,".seq.mumm"))
    dat <- getDataFrameFromMummerOutPut(paste0("./mummer/",sample_i,"_",sample_j,".seq.mumm"), paste0(files[j],".kraken"), paste0(files[i],".kraken"))
    write.table(dat,
                file = paste0("./mummer/",sample_i,"_",sample_j,".seq.mumm.csv"),
                append=FALSE,row.names=FALSE,col.names=TRUE,sep = "\t",quote=FALSE)
  } 
}

calculateNormalization <- function(generaTable1, generaTable2)
{
  generaUniqueTable1 <- data.frame(genus=unique(generaTable1$genus),length=0)
  rownames(generaUniqueTable1) <- generaUniqueTable1$genus
  for (genus in generaUniqueTable1$genus)
  {
    generaUniqueTable1[genus,"length"] <- sum(generaTable1$length[generaTable1$genus==genus])
  }
  generaUniqueTable2 <- data.frame(genus=unique(generaTable2$genus),length=0)
  rownames(generaUniqueTable2) <- generaUniqueTable2$genus
  for (genus in generaUniqueTable2$genus)
  {
    generaUniqueTable2[genus,"length"] <- sum(generaTable2$length[generaTable2$genus==genus])
  }
  Normalization <- 0
  for (i in 1:nrow(generaUniqueTable1))
  {
    for (j in 1:nrow(generaUniqueTable2))
    {
      if (generaUniqueTable1$genus[i]!=generaUniqueTable2$genus[j] & !generaUniqueTable1$genus[i] %in% c("Homo sapiens (taxid 9606)") & !generaUniqueTable1$genus[i] %in% c("Homo sapiens (taxid 9606)"))
      {
        Normalization <- Normalization + generaUniqueTable1$length[i]*generaUniqueTable2$length[j]
      }
    }
  }
  return(Normalization)
}


# make pdfs
files <- Sys.glob(file.path(".","*", "metaspades", "*contigs-final.fasta.4e4.fasta"))
for (i in 1:length(files)) 
{
  sample_i <- strsplit(files[i],"/")[[1]][2]
  kraken1Table <- read.table(paste0(files[i],".kraken"), header = FALSE, sep="\t")
  Genera1 <- sapply(1:nrow(kraken1Table),function(z){strsplit(kraken1Table$V3[z]," ")[[1]][1]})
  Lengths1 <- kraken1Table$V4
  generaTable1 <- data.frame(genus=Genera1,length=Lengths1)
  print(sample_i)
  sample_i <- strsplit(files[i],"/")[[1]][2]
  for (j in i:length(files))
  {
    sample_j <- strsplit(files[j],"/")[[1]][2]
    kraken2Table <- read.table(paste0(files[j],".kraken"), header = FALSE, sep="\t")
    Lengths_j <- kraken2Table$V4[!kraken2Table$V4 %in% "Homo sapiens (taxid 9606)"]
    LengthMatrix <- Lengths_i %o% Lengths_j
    Genera2 <- sapply(1:nrow(kraken2Table),function(z){strsplit(kraken2Table$V3[z]," ")[[1]][1]})
    Lengths2 <- kraken2Table$V4
    generaTable2 <- data.frame(genus=Genera2,length=Lengths2)
    ToNorm <- calculateNormalization(generaTable1, generaTable2)/ifelse(i==j,2,1)
    if (file.exists(paste0("./mummer/",sample_i,"_",sample_j,".seq.mumm.csv")) & file.size(paste0("./mummer/",sample_i,"_",sample_j,".seq.mumm.csv"))>1)
    {
      dat <- read.table(paste0("./mummer/",sample_i,"_",sample_j,".seq.mumm.csv"), header = TRUE, sep="\t")
      # L <- dat$r[dat$genus1!=dat$genus2 & dat$genus1!="Homo sapiens (taxid 9606)" & dat$genus2!="Homo sapiens (taxid 9606)"]
      L <- dat$r[dat$genus1!="Homo sapiens (taxid 9606)" & dat$genus2!="Homo sapiens (taxid 9606)"]
      if (length(unique(L))>2)
      {
        Ltable <- as.data.frame(table(L),stringsAsFactors=FALSE)
        Ltable$L <- as.numeric(Ltable$L)
        if (i==j) pdf(paste0("./plots/same/",sample_i,"_",sample_j,".pdf"))
        else pdf(paste0("./plots/diff/",sample_i,"_",sample_j,".pdf"))
        n <- 40
        p <- hist(log(L),n,plot=FALSE)
        while(sum(p$counts==0)>0)
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





