##################################################### classify contigs using blastn ##################################################### 

# download RefSeq database
# perl /home/m.sheinman/Development/Software/update_blastdb.pl --blastdb_version 5 --decompress ref_prok_rep_genomes
# perl /home/m.sheinman/Development/Software/update_blastdb.pl --blastdb_version 5 --decompress nt


# blasn contig and output hits
# export BLASTDB=$BLASTDB:/shared/scratch/m.sheinman/refseq/
library(foreach)
setwd("/home/m.sheinman/Development/HGTnew/data/processed/EscherichiaColi/")
files <- Sys.glob(file.path(".","*.fasta"))
# filter long contigs
for (file in files)
{
  print(file)
  Command <- ""
  Command <- paste0(Command,"blastn -db /shared/scratch/m.sheinman/refseq/ref_prok_rep_genomes ",
" -evalue 1e-5 -word_size 64 -num_threads 40 -max_target_seqs 100000000 -task megablast ",
" -outfmt '6 qacc qlen sacc slen qstart qend sstart send evalue bitscore score length pident staxid ssciname scomname sblastname sskingdom sstrand' ", 
" -query ",file,
" -out ",file,".blastnRefSeq")
  print(Command)
  system(Command)
}

library(stringr)
files <- Sys.glob(file.path(".","*", "metaspades", "*contigs-final.fasta.4e4.fasta"))
# choose accession number with largest total length of all hits
for (file in files)
{
  print(file)
  dat <- data.frame()
  blastTable <- read.table(paste0(file,".blastnRefSeq"), header = FALSE, sep="\t", fill = TRUE)
  colnames(blastTable) <- c("qacc", "qlen", "sacc", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "score", "length", "pident", "staxid", "ssciname", "scomname", "sblastname", "sskingdom")
  contigs <- unique(blastTable$qacc)
  for (contig in contigs)
  {
    Ind1 <- which(blastTable$qacc==contig)
    saccs <- unique(blastTable$sacc[Ind1])
    TotalLength <- rep(0,length(saccs))
    for (i in 1:length(saccs))
    {
      Ind2 <- which(blastTable$sacc[Ind1]==saccs[i])
      TotalLength[i] <- sum(blastTable$length[Ind1][Ind2])
    }
    IndBest <- which.max(TotalLength)
    accs <- saccs[IndBest]
    totallength <- TotalLength[IndBest]
    IndAcc <- which(blastTable$sacc==accs)[1]
    name <- blastTable$ssciname[IndAcc]
    slen <- blastTable$slen[IndAcc]
    qlen <- blastTable$qlen[IndAcc]
    staxid <- blastTable$staxid[IndAcc]
    sskingdom <- blastTable$sskingdom[IndAcc]
    coverage <- totallength/qlen
    genus <- strsplit(name," ")[[1]][1]
    genus <- str_replace(genus,"]","")
    genus <- str_replace(genus,"\\[","")
    
    dat <- rbind(dat, data.frame(contig=contig, contiglength=qlen, blastAccNum=accs, hitlength=slen, TotalAlignmentLength=totallength,name=name,staxid=staxid,sskingdom=sskingdom,coverage=coverage,genus))
  }
  write.table(dat,paste0(file,".blastnRefSeq.classify"), row.names = FALSE, sep = "\t", quote = FALSE)
}







