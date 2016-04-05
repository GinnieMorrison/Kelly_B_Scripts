#######################################################################
### Use this script to perform an analysis of the Golden Glow data  ###
### that mimics the Kelly, 2013 implementation of the Fisher & Ford ###
### transformation                                                  ###
#######################################################################
### Timothy M. Beissinger
### 11-1-2013
library(plyr)


clean_files<-function(filename){

   gen<-read.table(filename, sep='\t', header=T)
   gen<-gen[,c("Site.Name","Chromosome","Physical.Position","Number.of.Taxa","Major.Allele","Major.Allele.Frequency","Minor.Allele","Minor.Allele.Frequency","Gametes.Missing")]  ##This is me grabbing the fields that I want.
   cleanup<-which(gen$Major.Allele != "+" & gen$Major.Allele != '-' & gen$Minor.Allele != '+' & gen$Minor.Allele != '-' & gen$Chromosome != 0) ##No indels for allele 1 or 2
   gen<-gen[cleanup,]
   gen<-droplevels(gen)

   return(gen)
}

bind_across_gens<-function(gen1,gen2){ ##binding the two generations together based on Major allele call--so that the actual change can be measured
   #Major stays major
   
   joinMM<-join(gen1,gen2,by=c("Site.Name","Chromosome","Physical.Position","Major.Allele"),type="inner")
   names(joinMM)<-c("Site.Name","Chromosome","Physical.Position","Number.of.Taxa.Gen1","Major.Allele","Allele.Frequency.Gen1","Minor.Allele.1","Minor.Allele.Frequency.1","Gametes.Missing.Gen1","Number.of.Taxa.Gen2","Allele.Frequency.Gen2","Minor.Allele.2","Minor.Allele.Frequency.2","Gametes.Missing.Gen2")

   joinsites<-join(gen0,gen6,by=c("Site.Name","Chromosome","Physical.Position"),type="inner")
   names(joinsites)<-c("Site.Name","Chromosome","Physical.Position","Number.of.Taxa.Gen1","Major.Allele.1","Allele.Frequency.Gen1","Minor.Allele.1","Minor.Allele.Frequency.1","Gametes.Missing.Gen1","Number.of.Taxa.Gen2","Major.Allele.2","Major.Allele.Frequency.Gen2","Minor.Allele.2","Allele.Frequency.Gen2","Gametes.Missing.Gen2")
   
   joinMm<-joinsites[joinsites$Major.Allele.1==joinsites$Minor.Allele.2 & joinsites$Minor.Allele.1==joinsites$Major.Allele.2,]   

   #Take the fields of interest and apply the same names to matching fields so that the data.frames can be combined
   MM<-joinMM[,c("Site.Name","Chromosome","Physical.Position","Number.of.Taxa.Gen1","Allele.Frequency.Gen1","Gametes.Missing.Gen1","Number.of.Taxa.Gen2","Allele.Frequency.Gen2","Gametes.Missing.Gen2")]
   Mm<-joinMm[,c("Site.Name","Chromosome","Physical.Position","Number.of.Taxa.Gen1","Allele.Frequency.Gen1","Gametes.Missing.Gen1","Number.of.Taxa.Gen2","Allele.Frequency.Gen2","Gametes.Missing.Gen2")] 

   shared_sites<-rbind(MM,Mm)

   return(shared_sites)

}

### Load allele frequency data
#frequencyData <- read.table("/home/timbeissinger/Documents/GoldenGlowPE2/CalculateFrequencies/alleleFreq_filtered.txt",header=T,sep="\t",stringsAsFactors=F)
args<-commandArgs(trailingOnly=T)                                                                                                                                                        
f1<-args[1]
f2<-args[2]
outKelly_method<-args[3]
#f1<-"ZeaSyn0_default_no0_corrected3.txt"
#f2<-"ZeaSyn6_default_no0_corrected3.txt"
gen0<-clean_files(f1)
gen6<-clean_files(f2)
frequencyData<-bind_across_gens(gen0,gen6)

### Compute Fisher's angular transformation
frequencyData$xa <- 2*asin(sqrt(frequencyData$Allele.Frequency.Gen1)) #ancestral population
frequencyData$xh <- 2*asin(sqrt(frequencyData$Allele.Frequency.Gen2)) #in the kelly paper, h is for the 'high' population--equation is same for low population, we're just comparing one pop here

### The difference between x values should be normal NEED TO CHECK THIS
hist(frequencyData$xh-frequencyData$xa,freq=F,breaks=50) 

### Compute inter-quartile range and determine variance of the null model
va <- ((quantile(frequencyData$xa)[4]- quantile(frequencyData$xa)[2])/1.349)^2 ##Take interquartiles because want to get the null model. Where most of the data is. 1.349 specific to when have a normal distribution--calculate the variance of the interquartile range.
vh <- ((quantile(frequencyData$xh)[4]- quantile(frequencyData$xh)[2])/1.349)^2



### Compute snp-specific variance
frequencyData$gen0Counts<-frequencyData$Number.of.Taxa.Gen1-(frequencyData$Gametes.Missing.Gen1/2) #Just using # of individuals and not read counts 
frequencyData$gen6Counts<-frequencyData$Number.of.Taxa.Gen2-(frequencyData$Gametes.Missing.Gen2/2)

frequencyData$ua <- va + 1/frequencyData$gen0Counts
frequencyData$uh <- vh + 1/frequencyData$gen6Counts

### Try estimating uauh directly
uauh <- ((quantile(frequencyData$xh-frequencyData$xa)[4]- quantile(frequencyData$xh-frequencyData$xa)[2])/1.349)^2
##What is this value srun to find out and 
print(uauh)
### Compute single-snp test statistic
frequencyData$KellyTest1 <- (frequencyData$xh-frequencyData$xa)/sqrt(frequencyData$ua + frequencyData$uh)  # probably incorrect -- ua/uh overestimate var.

frequencyData$KellyTest2 <- (frequencyData$xh-frequencyData$xa)/sqrt(uauh)

### Compute single-snp p-values
frequencyData$Kelly.p1 <- pnorm(abs(frequencyData$KellyTest1),lower.tail=F)
frequencyData$Kelly.p2 <- pnorm(abs(frequencyData$KellyTest2),lower.tail=F)

### Identify bonferroni significant SNPs
bonferroni.p <- 0.05/nrow(frequencyData)
significantSNPs1 <- frequencyData[which(frequencyData$Kelly.p1 <= bonferroni.p),]
significantSNPs2 <- frequencyData[which(frequencyData$Kelly.p2 <= bonferroni.p),]

write.table(frequencyData,outKelly_method,sep='\t',row.names=F)
#write.table(significantSNPs1, "Kelly1_sig.txt", sep='\t',row.names=F)
#write.table(significantSNPs2, "Kelly2_sig.txt", sep='\t',row.names=F)
#####################################################################################################################
### Window approach
library(zoo)
library(GenWin)

frequencyData<-read.table("Kelly_Method_Syn0Syn6.txt",sep='\t',header=T)
### compute B for 25 SNP windows
frequencyData<-frequencyData[with(frequencyData,order(Chromosome,Physical.Position)),]
   png("Kelly_Syn0Syn6_Chr1.png")
   pos<-frequencyData[frequencyData$Chromosome==1,"Physical.Position"]
   B_chr1<-splineAnalyze(frequencyData[frequencyData$Chromosome==1,"KellyTest2"]^2,map=pos,method=3) 
   dev.off()
   Chr1_windows<-B_chr1$windowData
   Chr1_windows$Chromosome<-rep(1,dim(Chr1_windows)[1])
 
   png("Kelly_Syn0Syn6_Chr2.png")
   pos<-frequencyData[frequencyData$Chromosome==2,"Physical.Position"]
   B_chr2<-splineAnalyze(frequencyData[frequencyData$Chromosome==2,"KellyTest2"]^2,map=pos,method=3)
   dev.off()
   Chr2_windows<-B_chr2$windowData
   Chr2_windows$Chromosome<-rep(2,dim(Chr2_windows)[1])

   png("Kelly_Syn0Syn6_Chr3.png")
   pos<-frequencyData[frequencyData$Chromosome==3,"Physical.Position"]
   B_chr3<-splineAnalyze(frequencyData[frequencyData$Chromosome==3,"KellyTest2"]^2,map=pos,method=3)
   dev.off()
   Chr3_windows<-B_chr3$windowData
   Chr3_windows$Chromosome<-rep(3,dim(Chr3_windows)[1])
   
   png("Kelly_Syn0Syn6_Chr4.png")
   pos<-frequencyData[frequencyData$Chromosome==4,"Physical.Position"]
   B_chr4<-splineAnalyze(frequencyData[frequencyData$Chromosome==4,"KellyTest2"]^2,map=pos,method=3)
   dev.off()
   Chr4_windows<-B_chr4$windowData
   Chr4_windows$Chromosome<-rep(4,dim(Chr4_windows)[1])

   png("Kelly_Syn0Syn6_Chr5.png")
   pos<-frequencyData[frequencyData$Chromosome==5,"Physical.Position"]
   B_chr5<-splineAnalyze(frequencyData[frequencyData$Chromosome==5,"KellyTest2"]^2,map=pos,method=3)
   dev.off()
   Chr5_windows<-B_chr5$windowData
   Chr5_windows$Chromosome<-rep(5,dim(Chr5_windows)[1])

   png("Kelly_Syn0Syn6_Chr6.png")
   pos<-frequencyData[frequencyData$Chromosome==6,"Physical.Position"]
   B_chr6<-splineAnalyze(frequencyData[frequencyData$Chromosome==6,"KellyTest2"]^2,map=pos,method=3)
   dev.off()
   Chr6_windows<-B_chr6$windowData
   Chr6_windows$Chromosome<-rep(6,dim(Chr6_windows)[1])

   png("Kelly_Syn0Syn6_Chr7.png")
   pos<-frequencyData[frequencyData$Chromosome==7,"Physical.Position"]
   B_chr7<-splineAnalyze(frequencyData[frequencyData$Chromosome==7,"KellyTest2"]^2,map=pos,method=3)
   dev.off()
   Chr7_windows<-B_chr7$windowData
   Chr7_windows$Chromosome<-rep(7,dim(Chr7_windows)[1])

   png("Kelly_Syn0Syn6_Chr8.png")
   pos<-frequencyData[frequencyData$Chromosome==8,"Physical.Position"]
   B_chr8<-splineAnalyze(frequencyData[frequencyData$Chromosome==8,"KellyTest2"]^2,map=pos,method=3)
   dev.off()
   Chr8_windows<-B_chr8$windowData
   Chr8_windows$Chromosome<-rep(8,dim(Chr8_windows)[1])

   png("Kelly_Syn0Syn6_Chr9.png")
   pos<-frequencyData[frequencyData$Chromosome==9,"Physical.Position"]
   B_chr9<-splineAnalyze(frequencyData[frequencyData$Chromosome==9,"KellyTest2"]^2,map=pos,method=3)
   dev.off()
   Chr9_windows<-B_chr9$windowData
   Chr9_windows$Chromosome<-rep(9,dim(Chr9_windows)[1])

   png("Kelly_Syn0Syn6_Chr10.png")
   pos<-frequencyData[frequencyData$Chromosome==10,"Physical.Position"]
   B_chr10<-splineAnalyze(frequencyData[frequencyData$Chromosome==10,"KellyTest2"]^2,map=pos,method=3) 
   dev.off()
   Chr10_windows<-B_chr10$windowData
   Chr10_windows$Chromosome<-rep(10,dim(Chr10_windows)[1])

Windowed_Chrs<-rbind(Chr1_windows,Chr2_windows,Chr3_windows,Chr4_windows,Chr5_windows,Chr6_windows,Chr7_windows,Chr8_windows,Chr9_windows,Chr10_windows)
Windowed_Chrs$Number.of.Taxa.Gen1<-rep(94,dim(Windowed_Chrs)[1])
Windowed_Chrs$Number.of.Taxa.Gen2<-rep(1846,dim(Windowed_Chrs)[1])
Windowed_Chrs$Window.B<-Windowed_Chrs$SNPcount*Windowed_Chrs$MeanY
write.table(Windowed_Chrs,"Syn0Syn6_Windows_method3.txt",row.names=F,sep='\t')

Bowley.Skew <- function(data){
  skew <- as.numeric((quantile(data,na.rm=T)[4]+quantile(data,na.rm=T)[2]-2*quantile(data,na.rm=T)[3])/(quantile(data,na.rm=T)[4]-quantile(data,na.rm=T)[2]))
  return(skew)
}

skew.Bwin<- Bowley.Skew(Windowed_Chrs$Window.B)
### Write a function to determine Bowley's skewness for a chi-square distribution
Bowley.Skew.Chisq <- function(k){
 chi.skew <- (qchisq(0.75,df=k)+qchisq(0.25,df=k)-2*qchisq(0.5,df=k))/(qchisq(0.75,df=k)-qchisq(0.25,df=k))
 return(chi.skew)

}
#determine delta (the df of a chisq dist) that will give an equivalent value for Bowley's skew. This was/is a lot of trial and error. 
Bowley.Skew.Chisq(1.1274057)
delta_m4 <-1.1274057

Bowley.Skew.Chisq(1.2238045) 
delta_m3 <-1.2238045

sd.Bwin = sqrt(2*delta_m3)*((quantile(Windowed_Chrs$Window.B,na.rm=T)[4] -quantile(Windowed_Chrs$Window.B,na.rm=T)[2]) / (qchisq(0.75,df=delta_m3) - qchisq(0.25,df=delta_m3)))

Windowed_Chrs$Bstar <- delta_m3+sqrt(2*delta_m3)*( (Windowed_Chrs$Window.B - Windowed_Chrs$SNPcount) / sd.Bwin  )

Windowed_Chrs$Window.p <- pchisq(Windowed_Chrs$Bstar,df=delta_m3,lower.tail=F) # P.values  

Windowed_Chrs$Window.bonf <- p.adjust(Windowed_Chrs$Window.p,method="bonferroni")

write.table(Windowed_Chrs, "Kelly_GenWin_method3_Syn0Syn6.txt",row.name=F,sep='\t')
q()
n

