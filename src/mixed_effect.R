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
library(merTools)


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
      
      print(paste("new sp", species1,species2))
      if (dim(Prefactor)[1]>5){
        colnames(Prefactor) = c("country1","country2","A","Same")
        Prefactor2=Prefactor
        for (i in AllCountries){
          Prefactor2[i]=0
          Prefactor2[which(Prefactor2$country1 == i | Prefactor2$country2 == i ),i]=1
        }
        for (i in Species_list){
          Prefactor2[i]=0
        }
        Prefactor2$species1 = species1
        Prefactor2$species2 = species2
        
        Prefactor2[species1]=1
        Prefactor2[species2]=1
        for (i in AllContinents){
          Prefactor2[i]=0
        }
        
        for( ind in 1:dim(Prefactor2)[1]){
          country_sp1=as.character(Prefactor2$country1[ind])
          country_sp2=as.character(Prefactor2$country2[ind])
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

# AllPrefactors$A = -AllPrefactors$A

dim(AllPrefactors)
AllPrefactorsSaved = AllPrefactors

colnames(AllPrefactorsSaved)[c(114:127)]
AllPrefactors=AllPrefactorsSaved[,c(3,4,c(114:121),c(124:129))]

lm1=lm(data = AllPrefactors,A~Same)
summary(lm1)

lm2=lm(data = AllPrefactors,A~.)
summary(lm2)

lm2=lm(data = AllPrefactors,A~Same)
summary(lm2)

lm2=lm(data = AllPrefactors,A~Same+EscherichiaColi+Salmonella+KlebsiellaPneumoniae+
         Vibrio+Staphylococcus+Streptococcus+Pseudomonas+Mycobacterium)
summary(lm2)

colnames(AllPrefactors)
write.table(summary(lm2)$coefficients,file = "~/HGTnew/data/processed/results/test_all_Prefactors_per_species-filtered-mash.txt",
            sep="\t",quote = F)


p=ggplot(AllPrefactors,aes(x=Same,y=A,fill=Same))+
  geom_boxplot()+geom_jitter(cex=0.1)
  # +geom_signif(comparisons = list(c(TRUE, FALSE)))
p = p + stat_compare_means(method = "wilcox.test")
print(p)

colnames(AllPrefactors)
AllPrefactors2=AllPrefactors[,-c(2)]
lm1=lm(data = AllPrefactors2,A~.)
AllPrefactors$Acorrected = summary(lm1)$residuals

p=ggplot(AllPrefactors,aes(x=Same,y=Acorrected,fill=Same))+
  geom_boxplot()+geom_signif()
p = p + stat_compare_means(method = "wilcox.test")
print(p)


country_per_cont$Continent = as.character(country_per_cont$Continent)
AllPrefactors$Continent = NA
for ( jj in 1:dim(AllPrefactors)[1] ){
  tt = which( as.character(country_per_cont$Country) == as.character(AllPrefactorsSaved$country1[jj]) )
  AllPrefactors$Continent[jj] = country_per_cont$Continent[tt] 
}

AllPrefactorsSaved$Continent = as.factor(AllPrefactors$Continent)

p=ggplot(AllPrefactorsSaved,aes(x=Continent,y=A,fill=Continent))+
  geom_boxplot()
# +theme(legend.position = "None")
# +geom_signif()
p 

p=ggplot(AllPrefactorsSaved,aes(x=species1,y=A,fill=species1))+
  geom_boxplot()+theme(axis.text.x=element_text(angle=30))
# +theme(legend.position = "None")
# +geom_signif()
p 


AllPrefactors3=AllPrefactors[which(AllPrefactors$species1=="EscherichiaColi" 
                                   && AllPrefactors$species2 =="KlebsiellaPneumoniae"),]
save.image("~/HGTnew/data/processed/results/AllPrefactors-filtering-mash.Rda")


PrefactorsToPlot=c()
infoFact=c()
infoNum=c()
for (i in Species_list){
  for (j in country_per_cont$Continent){
  AA=which(AllPrefactorsSaved[,i] == 1 & AllPrefactorsSaved$Continent == j)
  toAdd = cbind(rep(i,length(AA)),rep(j,length(AA)))
  infoFact = rbind(infoFact,toAdd)
  toAdd = as.numeric(AllPrefactorsSaved$A[AA])
  infoNum=c(infoNum,toAdd)

  }
}


PrefactorsToPlot = data.frame("A"=infoNum)
PrefactorsToPlot$Species = infoFact[,1]
PrefactorsToPlot$Continent = infoFact[,2]
head(PrefactorsToPlot)
tail(PrefactorsToPlot)


p=ggplot(PrefactorsToPlot,aes(x=Continent,y=A,fill=Continent))+
  geom_boxplot(outlier.shape = NA)+theme(legend.position = "None")+
  theme(axis.text.x=element_text(angle=30))+
  facet_wrap(.~Species,nrow=2,scales="free")
#+geom_signif()
p

dim(PrefactorsToPlot)

################################### Correcting for Countries

colnames(AllPrefactorsSaved)
AllPrefactors2=AllPrefactorsSaved[,c(3:121)]

lm1=lm(data = AllPrefactors2,A~.)
summary(lm1)

colnames(AllPrefactorsSaved)

colnames(AllPrefactors2)
AllPrefactors3=AllPrefactors2[,-c(2)]
lm1=lm(data = AllPrefactors3,A~.)

AllPrefactors$Acorrected = summary(lm1)$residuals
summary(lm1)

p=ggplot(AllPrefactors,aes(x=Same,y=Acorrected,fill=Same))+
  geom_boxplot()+geom_signif()
p = p + stat_compare_means(method = "wilcox.test")
print(p)

lm2 = lm(data = AllPrefactors, Acorrected ~ Same) 
summary(lm2)

# lm3 = lmer(data = AllPrefactors, A~Same+ (1|as.character(colnames(AllPrefactors)[-c(1,2)])))


######## Mixed MOdels Continent

lmer_res = lmer(data = AllPrefactors, A~Same+ (1|Continent))
summary(lmer_res)


ranef(lmer_res)$Continent

#To be installed

predictInterval(lmer_res)   # for various model predictions, possibly with new data
REsim(lmer_res)
plotREsim(REsim(lmer_res))
predict(lmer_res, re.form=NA) %>% head()
predict(lmer_res) %>% head()

######## Mixed MOdels Species

colnames(AllPrefactorsSaved)

AA = AllPrefactorsSaved[,c("A","country1","country2","Same","species1","species2")]

lmer_res = lmer(data = AA, A~Same+ (1|species1) + (1|species2) + (1|country1)+ (1|country2) )
summary(lmer_res)
ranef(lmer_res)

predictInterval(lmer_res)   # for various model predictions, possibly with new data
REsim(lmer_res)
plotREsim(REsim(lmer_res))
predict(lmer_res, re.form=NA) %>% head()
predict(lmer_res) %>% head()

summary(lmer_res)
drop1(lmer_res,test = "Chisq")


## I removed "country2" because I think this is not identifiable if I keep it in
SameSame=which (as.character(AA$country1) == as.character(AA$country2))
AA$country2[SameSame]=NA
## Doesnt' work "rank defficient"

lmer_res2 = lmer(data = AA, A~Same+ (1|species1) + (1|species2) + (1|country1)+(1|country2))
lmer_res3 = lmer(data = AA, A~ (1|species1) + (1|species2) + (1|country1)+(1|country2))
summary(lmer_res2)
anova(lmer_res2,lmer_res3,test = "Chisq")
confint(lmer_res2)

plotREsim(REsim(lmer_res2))
drop1(lmer_res2,test = "Chisq")

