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
  system(paste0("/home/m.sheinman/Development/Software/seqkit/seqkit seq -m 10000 ",file," > ",file,".1e4.fasta"))
}


# copy krake2 database
# downloaded from https://benlangmead.github.io/aws-indexes/k2
# wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20210517.tar.gz
# cp /shared/scratch/m.sheinman/k2_standard/*k2d /dev/shm/
# cp /shared/scratch/m.sheinman/k2_standard/*kmer_distrib /dev/shm/


# taxonomic clasiffcation with CAT and BAT github.com/dutilh/CAT
# export PATH="/shared/scratch/m.sheinman/CAT-master/CAT_prepare_20210107/Diamond_2.0.6:$PATH"
# export PATH="/shared/scratch/m.sheinman/CAT-master/:$PATH"
setwd("/shared/scratch/m.sheinman/assembly/")
files <- Sys.glob(file.path(".","*", "metaspades", "*contigs-final.fasta.1e5.fasta"))
for (i in 1:length(files))
{
  dir <- dirname(files[i]) 
  Command <- paste0("/shared/scratch/m.sheinman/CAT-master/CAT_pack/CAT contigs --fraction 0.99 --force -n 40 ",
                    "-c ",files[i],
  " -o ",files[i],".CAT ",
  " -d /shared/scratch/m.sheinman/CAT-master/CAT_prepare_20210107/2021-01-07_CAT_database ",
  " -t /shared/scratch/m.sheinman/CAT-master/CAT_prepare_20210107/2021-01-07_taxonomy/ ",
  " --tmpdir ",dir)
  system(Command)
  
  Command <- paste0("/shared/scratch/m.sheinman/CAT-master/CAT_pack/CAT add_names --force --only_official ",
  "-i ",files[i],".CAT.contig2classification.txt ",
  "-o ",files[i],".CAT.contig2classification.names.txt ",
  "-t /shared/scratch/m.sheinman/CAT-master/CAT_prepare_20210107/2021-01-07_taxonomy/")
  system(Command)
  
  Command <- paste0("/shared/scratch/m.sheinman/CAT-master/CAT_pack/CAT add_names --force --only_official ",
                    "-i ",files[i],".CAT.ORF2LCA.txt ",
                    "-o ",files[i],".CAT.ORF2LCA.names.txt ",
                    "-t /shared/scratch/m.sheinman/CAT-master/CAT_prepare_20210107/2021-01-07_taxonomy/")
  system(Command)
}
# add lengths of contigs to CAT output and remove scores
for (i in 1:length(files)) 
{
  if (file.exists(paste0(files[i],".CAT.contig2classification.names.txt")))
  {
    CAT1Table <- read.table(paste0(files[i],".CAT.contig2classification.names.txt"), header = TRUE, sep="\t",comment.char = "")
    colnames(CAT1Table)[1] <- "contig"
    fi <- read.fasta(file = files[i])
    CAT1Table$length <- getLength(fi)
    CAT1Table$phylum <- sapply(1:nrow(CAT1Table),function(z){foo <- strsplit(CAT1Table$phylum[z],":")[[1]][1]; ifelse(is.na(foo),"no support",foo)})
    CAT1Table$class <- sapply(1:nrow(CAT1Table),function(z){foo <- strsplit(CAT1Table$class[z],":")[[1]][1]; ifelse(is.na(foo),"no support",foo)})
    CAT1Table$order <- sapply(1:nrow(CAT1Table),function(z){foo <- strsplit(CAT1Table$order[z],":")[[1]][1]; ifelse(is.na(foo),"no support",foo)})
    CAT1Table$family <- sapply(1:nrow(CAT1Table),function(z){foo <- strsplit(CAT1Table$family[z],":")[[1]][1]; ifelse(is.na(foo),"no support",foo)})
    CAT1Table$genus <- sapply(1:nrow(CAT1Table),function(z){foo <- strsplit(CAT1Table$genus[z],":")[[1]][1]; ifelse(is.na(foo),"no support",foo)})
    CAT1Table$species <- sapply(1:nrow(CAT1Table),function(z){foo <- strsplit(CAT1Table$species[z],":")[[1]][1]; ifelse(is.na(foo),"no support",foo)})
    write.table(CAT1Table,paste0(files[i],".CAT.contig2classification.names.lengths.txt"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
  }
  else
  {
    print(files[i])
  }
}



  

# function that outputs data.frame table from mummer file
getDataFrameFromMummerOutPut <- function(mummerfile,CAT1Table,CAT2Table)
{
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
        Ind1 <- which(CAT1Table$contig==seq2)
        Ind2 <- which(CAT2Table$contig==seq1)
        L1 <- CAT1Table$length[Ind1]
        L2 <- CAT2Table$length[Ind2]
        if (CAT1Table$genus[Ind1] != CAT2Table$genus[Ind2] & CAT1Table$genus[Ind1]!="no support" & CAT2Table$genus[Ind2]!="no support") {rank <- "genus"; taxon1 <- CAT1Table$genus[Ind1]; taxon2 <- CAT2Table$genus[Ind2]
        }else if (CAT1Table$family[Ind1] != CAT2Table$family[Ind2] & CAT1Table$family[Ind1]!="no support" & CAT2Table$family[Ind2]!="no support") {rank <- "family"; taxon1 <- CAT1Table$family[Ind1]; taxon2 <- CAT2Table$family[Ind2]
        }else if (CAT1Table$order[Ind1] != CAT2Table$order[Ind2] & CAT1Table$order[Ind1]!="no support" & CAT2Table$order[Ind2]!="no support") {rank <- "order"; taxon1 <- CAT1Table$order[Ind1]; taxon2 <- CAT2Table$order[Ind2]
        }else if (CAT1Table$class[Ind1] != CAT2Table$class[Ind2] & CAT1Table$class[Ind1]!="no support" & CAT2Table$class[Ind2]!="no support") {rank <- "class"; taxon1 <- CAT1Table$class[Ind1]; taxon2 <- CAT2Table$class[Ind2]
        }else if (CAT1Table$phylum[Ind1] != CAT2Table$phylum[Ind2] & CAT1Table$phylum[Ind1]!="no support" & CAT2Table$phylum[Ind2]!="no support") {rank <- "phylum"; taxon1 <- CAT1Table$phylum[Ind1]; taxon2 <- CAT2Table$phylum[Ind2]
        }else {rank <- "none"; taxon1 <- paste(CAT1Table$phylum[Ind1],CAT1Table$class[Ind1],CAT1Table$order[Ind1],CAT1Table$family[Ind1],CAT1Table$genus[Ind1],collapse="_"); taxon2 <- paste(CAT2Table$phylum[Ind2],CAT2Table$class[Ind2],CAT2Table$order[Ind2],CAT2Table$family[Ind2],CAT2Table$genus[Ind2],collapse="_");}
        datNew <- data.frame(seq1=seq1,x1=x1,L1=L1,taxon1=taxon1,seq2=seq2,x2=x2,L2=L2,taxon2=taxon2,r=r,reverse=reverse,seq=seq,rank=rank)
        dat <- rbind(dat,datNew)
      }
    }
  }
  return(dat[dat$rank!="none",])
}

# outputs data.frame table from mummer file
files <- Sys.glob(file.path(".","*", "metaspades", "*contigs-final.fasta.1e4.fasta"))
foreach (i = 1:length(files)) %dopar%
{
  sample_i <- strsplit(files[i],"/")[[1]][2]
  print(sample_i)
  CAT1Table <- read.table(paste0(files[i],".CAT.contig2classification.names.lengths.txt"), header = TRUE, sep="\t",comment.char = "")
  
  for (j in i:length(files))
  {
    sample_j <- strsplit(files[j],"/")[[1]][2]
    mummerfile <- paste0("./mummer/",sample_i,"_",sample_j,".seq.mumm")
    CAT2Table <- read.table(paste0(files[j],".CAT.contig2classification.names.lengths.txt"), header = TRUE, sep="\t",comment.char = "")
    # system(paste0("/home/m.sheinman/Development/Software/MUMmer3.23/mummer -maxmatch -n -s -b -l 100 -F ",files[i]," ",files[j]," > ./mummer/",sample_i,"_",sample_j,".seq.mumm"))
    dat <- getDataFrameFromMummerOutPut(mummerfile, CAT1Table, CAT2Table)
    write.table(dat,
                file = paste0("./mummer/",sample_i,"_",sample_j,".seq.mumm.csv"),
                append=FALSE,row.names=FALSE,col.names=TRUE,sep = "\t",quote=FALSE)
  } 
}

calculateNormalization <- function(CAT1Table, CAT2Table)
{
  Normalization <- 0
  for (i in 1:nrow(CAT1Table))
  {
    for (j in 1:nrow(CAT2Table))
    {
      if (CAT1Table$phylum[i]!=CAT2Table$phylum[j] | CAT1Table$class[i]!=CAT2Table$class[j] | CAT1Table$order[i]!=CAT2Table$order[j] | CAT1Table$family[i]!=CAT2Table$family[j] | CAT1Table$genus[i]!=CAT2Table$genus[j] )
      {
        Normalization <- Normalization + as.double(CAT1Table$length[i])*as.double(CAT2Table$length[j])
      }
    }
  }
  return(Normalization)
}




# make pdfs
files <- Sys.glob(file.path(".","*", "metaspades", "*contigs-final.fasta.4e4.fasta"))[1:12]
for (i in 1:length(files)) 
{
  sample_i <- strsplit(files[i],"/")[[1]][2]
  CAT1Table <- read.table(paste0(files[i],".CAT.contig2classification.names.lengths.txt"), header = TRUE, sep="\t",comment.char = "")
  for (j in i:length(files))
  {
    sample_i <- strsplit(files[i],"/")[[1]][2]
    CAT2Table <- read.table(paste0(files[j],".CAT.contig2classification.names.lengths.txt"), header = TRUE, sep="\t",comment.char = "")
    ToNorm <- calculateNormalization(CAT1Table, CAT2Table)/ifelse(i==j,2,1)
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





