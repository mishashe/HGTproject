library(DBI)
library(restez)
library(phangorn)
library(R.utils)
library(stringr)
restez_disconnect()

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1 | args[1] == '--help'){
  message("usage: downloadGenBank.R species")
  stop()
}

species=args[1]
species_dir=args[2]
country_file_lm=args[3]
country_file_hm=args[4]



# import data from JGI
Table <- read.table(paste(file="~/HGTnew/data/external/",species,".csv",sep='')
                    ,sep="\t", row.names=NULL,header=TRUE)
restez_path_set("~/HGTnew/data/external")
restez_connect()

#species_dir=paste("~/HGTnew/data/processed/",species,sep="")
comm=paste("mkdir -p", species_dir)
system(comm)
# get species from all countries
countries <- unique(Table$Isolation.Country)
countries <- str_replace(countries," ","");countries <- str_replace(countries," ","");countries <- str_replace(countries," ","")
countries <- countries[countries!=""]
country_highmem=c()
country_lowmem=c()
for (country in countries)
{
  print(country)
  ACCN <- Table$NCBI.Biosample.Accession[Table$Isolation.Country==country & Table$High.Quality=="Yes"]
  ACCN <- ACCN[ACCN!=""]
  i <- 1
  write(data.frame(), paste0(species_dir,"/",species,"_",country,".fasta"),append=FALSE)
  while (i<=length(ACCN))
  {
    term <- paste(ACCN[i:min(i+4,length(ACCN))],collapse="[ACCN] OR ")
    IDs <- rentrez::entrez_search(db   = "biosample",term = term,retmax=101, use_history=FALSE)$ids
    linked_seq_ids <- rentrez::entrez_link(dbfrom="biosample", id=IDs, db="nuccore",retmax=40000,term="(100000[SLEN] : 1000000000000000[SLEN])", use_history=FALSE)
    
    if(is.null(linked_seq_ids$links$biosample_nuccore)==F){
     seq <- rentrez::entrez_fetch(db="nuccore", id=linked_seq_ids$links$biosample_nuccore, rettype="fasta",retmax=101)
      write(seq, paste0(species_dir,"/",species,"_",country,".fasta"),append=TRUE)
      print(paste0(i," from ",length(ACCN)))
    }
    i <- i+5
    }
    fileSize=file.info(paste0(species_dir,"/",species,"_",country,".fasta"))$size
    if (fileSize >3*10^8){
	country_highmem=c(country_highmem,country)
    }else if (fileSize >100){
	country_lowmem=c(country_lowmem,country)
    }
}
restez_disconnect()

write.table(x = country_highmem,file = country_file_hm,quote=F,row.names = F,col.names=F)
write.table(x = country_lowmem,file = country_file_lm,quote=F,row.names = F,col.names=F)








