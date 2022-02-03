
library(base)
library(stringr)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1 | args[1] == '--help'){
  message("usage: Rscript calculate_fasta_length.R species1  path")
  stop()
}

species1=args[1]
mypath=args[2]

files1 <- Sys.glob(file.path(mypath,species1,paste0("*",species1,"_*.fasta")))
species1_dir=paste(mypath,"/",species1,sep="")

for (i in 1:length(files1)) {
  size=file.info(files1[i])$size
  print(files1[i])
  if (length(size)>0){
    print(size)
   if (size>10){
    fi <- read.fasta(file = files1[i])
    Li <- sum(getLength(fi))
    write(Li, paste0(files1[i],".L"),append=FALSE)
   }
  }
}
