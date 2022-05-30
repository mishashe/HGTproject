


#allSpecies=c(#"Pseudomonas",
#"EscherichiaColi")


allSpecies=c("EscherichiaColi",
"KlebsiellaPneumoniae","Mycobacterium",
"Pseudomonas",
"Salmonella",
"Staphylococcus",
"Streptococcus",
"Vibrio")
#
for (species1 in allSpecies){

files1 <- Sys.glob(file.path("~/HGTnew","data","dist",species1,paste0("*"),"dist-all.txt"))

hh=rep(0,200)


cc=0
for (i in files1){
    a=read.table(i)
    h1=hist(a$V3,breaks = seq(0,1,0.001),plot=F)
    hh=hh+h1$counts
    cc=cc+1
    print(i)
    print(paste(cc,"/",length(files1)))
    print(range(a$V3))
    print(length(h1$counts))
    
    if (  cc / 50 == round(cc/50)){
      h1$counts=hh
      pdf(paste("~/HGTnew/plots/pairwise_distances/",species1,".pdf",sep=''))
      plot(h1,xlim=c(0,0.15),ylim=c(0,max(h1$counts[1:70])),
           main = paste(cc,"files -- ",species1),
           xlab="distance")
      dev.off()
    }
    
  }
      h1$counts=hh
      pdf(paste("~/HGTnew/plots/pairwise_distances/",species1,".pdf",sep=''))
      plot(h1,xlim=c(0,0.15),ylim=c(0,max(h1$counts[1:70])),
           main = paste(cc,"files -- ",species1),
           xlab="distance")
      dev.off()

}
