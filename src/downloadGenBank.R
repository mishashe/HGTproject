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
country_file=args[2]


# import data from JGI
Table <- read.table(paste(file="~/HGTnew/data/external/",species,".csv",sep='')
                    ,sep="\t", row.names=NULL,header=TRUE)
restez_path_set("~/HGTnew/data/external")
restez_connect()

species_dir=paste("~/HGTnew/data/processed/",species,sep="")
comm=paste("mkdir -p", species_dir)
system(comm)
# get species from all countries
countries <- unique(Table$Isolation.Country)
countries <- str_replace(countries," ","");countries <- str_replace(countries," ","");countries <- str_replace(countries," ","")
countries <- countries[countries!=""]
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
    seq <- rentrez::entrez_fetch(db="nuccore", id=linked_seq_ids$links$biosample_nuccore, rettype="fasta",retmax=101)
    write(seq, paste0(species_dir,"/",species,"_",country,".fasta"),append=TRUE)
    print(paste0(i," from ",length(ACCN)))
    i <- i+5
    }
}
restez_disconnect()

write.table(x = countries,file = country_file,quote=F,row.names = F)







