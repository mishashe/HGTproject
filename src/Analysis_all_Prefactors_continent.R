setwd("~/HGTnew/")
library(base)
# library(foreach)
library(stringr)
library(seqinr)
library(plotrix)
library(ComplexHeatmap)
library(ggpubr)


country_per_cont=read.table("~/HGTnew/data/processed/All_countries_file.txt",head=T)
country_per_cont=country_per_cont[order(country_per_cont$Continent),]

AllContinents=unique(country_per_cont$Continent)
AllContinents=AllContinents[c(6,1,2,3,4,5,7)]

Species_list=read.table("~/HGTnew/HGTproject/src/Species_list",stringsAsFactors = F)$V1

AllPrefactors=c()
for (species1 in Species_list){
  for (species2 in Species_list){
  if (species1 != species2  &&    all(sort(c(species1,species2)) == c(species1,species2)) ){
    load(paste0("~/HGTnew/data/processed/results/",species1,"_",species2,"_perCont.RData"))
    
    if (dim(Prefactor)[1]>5){
      Prefactor2=Prefactor
      for (i in AllContinents){
        Prefactor2[i]=0
        Prefactor2[which(Prefactor2$species1 == i | Prefactor2$species2 == i ),i]=1
      }
      for (i in Species_list){
        Prefactor2[i]=0
      }
      Prefactor2[species1]=1
      Prefactor2[species2]=1
    
      AllPrefactors=rbind(AllPrefactors,Prefactor2)
        
    }else{
    print("Nothing")}
    }
  }
}    

dim(AllPrefactors)
colnames(AllPrefactors)
AllPrefactors=AllPrefactors[,-c(1,2)]

lm1=lm(data = AllPrefactors,log(-A)~NorthAmerica + Africa + Asia + Australia + Europa + Antartica + Same)
summary(lm1)

lm1=lm(data = AllPrefactors,log(-A)~.)
summary(lm1)

write.table(summary(lm1)$coefficients,file = "~/HGTnew/data/processed/results/test_all_Prefactors.txt",
            sep="\t",quote = F)
