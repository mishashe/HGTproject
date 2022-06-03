setwd("~/HGTnew/")
library(base)
# library(foreach)
library(stringr)
library(seqinr)
library(plotrix)
library(ComplexHeatmap)
library(ggpubr)
library(ggplot2)
library(ggsignif)
library(lme4)

country_per_cont=read.table("~/HGTnew/data/processed/All_countries_file.txt",head=T)
country_per_cont=country_per_cont[order(country_per_cont$Continent),]
AllContinents=unique(country_per_cont$Continent)
AllContinents=AllContinents[c(6,1,2,3,4,5,7)]


##################
## Per Country
##############

country_per_cont=read.table("~/HGTnew/data/processed/All_countries_file.txt",head=T)

AllCountries=unique(country_per_cont$Country)

Species_list=read.table("~/HGTnew/HGTproject/src/Species_list",stringsAsFactors = F)$V1

AllPrefactors=c()
for (species1 in Species_list){
  for (species2 in Species_list){
    if (species1 != species2  &&    all(sort(c(species1,species2)) == c(species1,species2)) ){
      load(paste0("~/HGTnew/data/processed/results/",species1,"_",species2,"-filtered-mash.RData"))
      
      if (dim(Prefactor)[1]>5){
        Prefactor2=Prefactor
        for (i in AllCountries){
          Prefactor2[i]=0
          Prefactor2[which(Prefactor2$species1 == i | Prefactor2$species2 == i ),i]=1
        }
        for (i in Species_list){
          Prefactor2[i]=0
        }
        Prefactor2[species1]=1
        Prefactor2[species2]=1
        for (i in AllContinents){
          Prefactor2[i]=0
        }
        
        for( ind in 1:dim(Prefactor2)[1]){
          country_sp1=as.character(Prefactor2$species1[ind])
          country_sp2=as.character(Prefactor2$species2[ind])
          tt=which( as.character(country_per_cont$Country) == country_sp1)
          contSp1=country_per_cont$Continent[tt]
          Prefactor2[ind,paste(contSp1)]=1
          tt2=which( as.character(country_per_cont$Country) == country_sp2)
          contSp2=country_per_cont$Continent[tt2]
          Prefactor2[ind,paste(contSp2)]=1
        }
        
        AllPrefactors=rbind(AllPrefactors,Prefactor2)
        
      }else{
        print("Nothing")}
    }
  }
}    

dim(AllPrefactors)
AllPrefactorsSaved = AllPrefactors

colnames(AllPrefactorsSaved)
AllPrefactors=AllPrefactorsSaved[,c(3,4,c(114:127))]

lm2=lm(data = AllPrefactors,A~.)
summary(lm2)

lm2=lm(data = AllPrefactors,A~Same)
summary(lm2)

lm2=lm(data = AllPrefactors,A~Same+EscherichiaColi+Salmonella+KlebsiellaPneumoniae+
         Vibrio+Staphylococcus+Streptococcus+Pseudomonas+Mycobacterium)
summary(lm2)

colnames(AllPrefactors)
write.table(summary(lm1)$coefficients,file = "~/HGTnew/data/processed/results/test_all_Prefactors_per_Countries.txt",
            sep="\t",quote = F)


p=ggplot(AllPrefactors,aes(x=Same,y=A,fill=Same))+
  geom_boxplot()+geom_signif()
p

colnames(AllPrefactors)
AllPrefactors2=AllPrefactors[,-c(2)]
lm1=lm(data = AllPrefactors2,A~.)
AllPrefactors$Acorrected = summary(lm1)$residuals

p=ggplot(AllPrefactors,aes(x=Same,y=Acorrected,fill=Same))+
  geom_boxplot()+geom_signif()
p

country_per_cont$Continent = as.character(country_per_cont$Continent)
AllPrefactors$Continent = NA
for ( jj in 1:dim(AllPrefactors)[1] ){
  tt = which( as.character(country_per_cont$Country) == as.character(AllPrefactorsSaved$species1[jj]) )
  AllPrefactors$Continent[jj] = country_per_cont$Continent[tt] 
}

AllPrefactorsSaved$Continent = as.factor(AllPrefactors$Continent)

p=ggplot(AllPrefactorsSaved,aes(x=,y=A,fill=Continent))+
  geom_boxplot()
# +theme(legend.position = "None")
# +geom_signif()
p 

AllPrefactors3=AllPrefactors[which(AllPrefactors$species1=="EscherichiaColi" 
                                   && AllPrefactors$species2 =="KlebsiellaPneumoniae"),]
save.image("~/HGTnew/data/processed/results/AllPrefactors-filtering-mash.Rda")


PrefactorsToPlot = data.frame("A"=c(),"Species"=c(),"Continent"=c())
for (i in Species_list){
  for (j in country_per_cont$Continent){
  AA=which(AllPrefactorsSaved[,i] == 1 & AllPrefactorsSaved$Continent == j)
  toAdd = cbind(AllPrefactorsSaved$A[AA],rep(i,length(AA)),rep(j,length(AA)))
  PrefactorsToPlot=rbind(PrefactorsToPlot,toAdd)
  }
}
head(PrefactorsToPlot)
PrefactorsToPlot$V2 = as.factor(PrefactorsToPlot$V2)
PrefactorsToPlot$V1 = as.numeric(PrefactorsToPlot$V1)
p=ggplot(PrefactorsToPlot,aes(x=V3,y=V1,fill=V3))+
  geom_boxplot(outlier.shape = NA)+theme(legend.position = "None")+
  theme(axis.text.x=element_text(angle=30))+
  facet_wrap(.~V2,nrow=2,scales="free")
#+geom_signif()
p
dim(PrefactorsToPlot)
min(which(PrefactorsToPlot$V2=="EscherichiaColi"))


################################### Correcting for Countries

colnames(AllPrefactorsSaved)
AllPrefactors2=AllPrefactorsSaved[,c(3:121)]

colnames(AllPrefactorsSaved)

AllPrefactors2=AllPrefactorsSaved[,c(3:121)]
# lm3 = lmer(data = AllPrefactors, A~Same+ (1|as.character(colnames(AllPrefactors)[-c(1,2)])))

 lmer_res = lmer(data = AllPrefactors, A~Same+ (1|Continent))
summary(lmer_res)

anova(lm3)
ranef(lmer_res)$Continent

#To be installed
#library(merTools)
