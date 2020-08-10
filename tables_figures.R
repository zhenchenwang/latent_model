library("RevoScaleR")
library(summarytools)
library(LaplacesDemon)
library(kernlab)
library(dplyr)


population <- read.delim("G:/cvdgt.txt")
#Experiments by order as presented in the paper

#Experiment I (Figure 4 A and B) comparing GT vs SY without missing data using 
#1. no. of distinct values 
#2. distribution diagram
#3. missing data rate

gtSample <- read.delim("A:/Wangz/R Scripts/Ntr/12g.txt")
synSample <- read.delim("A:/Wangz/R Scripts/Ntr/12s.txt")


synSample.nomissing<-na.omit(synSample)
gtSample.nomissing<-na.omit(gtSample)
size<-min(c(nrow(synSample.nomissing),nrow(gtSample.nomissing)),na.rm=FALSE) 
GTM.TValues<-sample_n(gtSample.nomissing, size) 
SYM.TValues<-sample_n(synSample.nomissing,  size) 
print(dfSummary(GTM.TValues) , method = "render")
print(dfSummary(SYM.TValues) , method = "render")

plot(density(GTM.TValues$age), xlab="values", ylab="probability")
plot(density(SYM.TValues$age), xlab="values", ylab="probability")

plot(density(GTM.TValues$bmi), xlab="values", ylab="probability")
plot(density(SYM.TValues$bmi), xlab="values", ylab="probability")

plot(density(GTM.TValues$choleratio), xlab="values", ylab="probability")
plot(density(SYM.TValues$choleratio), xlab="values", ylab="probability")


#Experiment IIa (Figure 4 C and D) comparing GT vs SY with missing data using latent
#1. no. of distinct values 
#2. distribution diagram
#3. missing data rate
gtSample <- read.delim("A:/Wangz/R Scripts/Ntr/11g.txt")
synSample <- read.delim("A:/Wangz/R Scripts/Ntr/11s.txt")
size<-min(c(nrow(synSample),nrow(gtSample)),na.rm=FALSE) 
GT.TValues<-sample_n(gtSample, size) 
SY.TValues<-sample_n(synSample,  size) 
print(dfSummary(GT.TValues) , method = "render")
print(dfSummary(SY.TValues) , method = "render")

plot(density(GT.TValues$age), xlab="values", ylab="probability")
plot(density(na.omit(SY.TValues$age)), xlab="values", ylab="probability")

plot(density(na.omit(GT.TValues$bmi)), xlab="values", ylab="probability")
plot(density(na.omit(SY.TValues$bmi)), xlab="values", ylab="probability")

plot(density(na.omit(GT.TValues$choleratio)), xlab="values", ylab="probability")
plot(density(na.omit(SY.TValues$choleratio)), xlab="values", ylab="probability")

#Experiment IIb (Figure 4 E ) comparing GT vs SY with missing data (missing node method) using 
#1. no. of distinct values 
#2. distribution diagram
#3. missing data rate 


plot(density(SY.TValues$age), xlab="values", ylab="probability") 

plot(density(na.omit(SY.TValues$bmi)), xlab="values", ylab="probability")
 
plot(density(na.omit(SY.TValues$choleratio)), xlab="values", ylab="probability")



#Experiment III (Table II) Chi-square test for categorical variable 

test<-chisq.test(GT.TValues$strokeha, SY.TValues$strokeha) 
test[3]$p.value
test<-chisq.test(GT.TValues$af, SY.TValues$af) 
test[3]$p.value
test<-chisq.test(GT.TValues$atyantip, SY.TValues$atyantip) 
test[3]$p.value
test<-chisq.test(GT.TValues$steroid, SY.TValues$steroid) 
test[3]$p.value
test<-chisq.test(GT.TValues$impot, SY.TValues$impot) 
test[3]$p.value
test<-chisq.test(GT.TValues$migr, SY.TValues$migr) 
test[3]$p.value
test<-chisq.test(GT.TValues$ra, SY.TValues$ra) 
test[3]$p.value
test<-chisq.test(GT.TValues$ckidney, SY.TValues$ckidney) 
test[3]$p.value
test<-chisq.test(GT.TValues$semi, SY.TValues$semi) 
test[3]$p.value
test<-chisq.test(GT.TValues$sle, SY.TValues$sle) 
test[3]$p.value
test<-chisq.test(GT.TValues$treathyp, SY.TValues$treathyp) 
test[3]$p.value
test<-chisq.test(GT.TValues$type1, SY.TValues$type1) 
test[3]$p.value
test<-chisq.test(GT.TValues$type2, SY.TValues$type2) 
test[3]$p.value
test<-chisq.test(GT.TValues$ethr, SY.TValues$ethr) 
test[3]$p.value
test<-chisq.test(GT.TValues$smoking, SY.TValues$smoking) 
test[3]$p.value
test<-chisq.test(GT.TValues$fh_cad, SY.TValues$fh_cad) 
test[3]$p.value
test<-chisq.test(GT.TValues$gender, SY.TValues$gender) 
test[3]$p.value
test<-chisq.test(GT.TValues$region, SY.TValues$region) 
test[3]$p.value

 
#Experiment IV (Table III) KS test for categorical variable 
gtSample <- read.delim("A:/Wangz/R Scripts/Ntr/12g.txt")
synSample <- read.delim("A:/Wangz/R Scripts/Ntr/12s.txt")


synSample.nomissing<-na.omit(synSample)
gtSample.nomissing<-na.omit(gtSample)
size<-1000
GT.TValues<-sample_n(gtSample.nomissing, size) 
SY.TValues<-sample_n(synSample.nomissing,  size) 

ks_test_age<- ks.test( SY.TValues$age, GT.TValues$age)
ks_test_age 


ks_test_cholest<- ks.test(GT.TValues$choleratio, SY.TValues$choleratio)
ks_test_cholest 


ks_test_bmi<- ks.test(GT.TValues$bmi,SY.TValues$bmi)
ks_test_bmi 


ks_test_sbps<- ks.test(GT.TValues$sbps, SY.TValues$sbps)
ks_test_sbps 

ks_test_sbp<- ks.test(gtSample.1$sbp, synSample.1$sbp)
ks_test_sbp 


#ks test for KLD values
KLDS <- read.csv("A:/Wangz/R Scripts/Ntr/KLDS.txt")
ks_test_kldage<- ks.test(KLDS$age_sy, KLDS$age_gt)
ks_test_kldbmi<- ks.test(KLDS$bmi_sy, KLDS$bmi_gt)
ks_test_kldchol<- ks.test(KLDS$chol_sy, KLDS$chol_gt)
ks_test_kldsbp<- ks.test(KLDS$sbp_sy, KLDS$sbp_gt)
ks_test_kldsbps<- ks.test(KLDS$sbps_sy, KLDS$sbps_gt)

#Experiment V (Table IV ) comparing GT vs SY KLD with 10 iterations of 100,000 sample patients per itration.
#base population
Q3readInI <- read.delim("G:/cvdgt.txt")

#samples
gtSample.1<-read.delim("A:/Wangz/R Scripts/Ntr/1000000/1g.txt")
synSample.1<- read.delim("A:/Wangz/R Scripts/Ntr/1000000/1s.txt")

 
gtSample.2<-read.delim("A:/Wangz/R Scripts/Ntr/1000000/2g.txt")
synSample.2<- read.delim("A:/Wangz/R Scripts/Ntr/1000000/2s.txt")

gtSample.3<-read.delim("A:/Wangz/R Scripts/Ntr/1000000/3g.txt")
synSample.3<- read.delim("A:/Wangz/R Scripts/Ntr/1000000/3s.txt")

gtSample.4<-read.delim("A:/Wangz/R Scripts/Ntr/1000000/4g.txt")
synSample.4<- read.delim("A:/Wangz/R Scripts/Ntr/1000000/4s.txt")


gtSample.5<-read.delim("A:/Wangz/R Scripts/Ntr/1000000/5g.txt")
synSample.5<- read.delim("A:/Wangz/R Scripts/Ntr/1000000/5s.txt")

gtSample.6<-read.delim("A:/Wangz/R Scripts/Ntr/1000000/6g.txt")
synSample.6<- read.delim("A:/Wangz/R Scripts/Ntr/1000000/6s.txt")

gtSample.7<-read.delim("A:/Wangz/R Scripts/Ntr/1000000/7g.txt")
synSample.7<- read.delim("A:/Wangz/R Scripts/Ntr/1000000/7s.txt")

gtSample.8<-read.delim("A:/Wangz/R Scripts/Ntr/1000000/8g.txt")
synSample.8<- read.delim("A:/Wangz/R Scripts/Ntr/1000000/8s.txt")

gtSample.9<-read.delim("A:/Wangz/R Scripts/Ntr/1000000/9g.txt")
synSample.9<- read.delim("A:/Wangz/R Scripts/Ntr/1000000/9s.txt")

gtSample.10<-read.delim("A:/Wangz/R Scripts/Ntr/1000000/10g.txt")
synSample.10<- read.delim("A:/Wangz/R Scripts/Ntr/1000000/10s.txt") 


size<-c(min(nrow(gtSample.1),nrow(synSample.1),na.rm=FALSE),
        min(nrow(gtSample.2),nrow(synSample.2),na.rm=FALSE),
        min(nrow(gtSample.3),nrow(synSample.3),na.rm=FALSE),
        min(nrow(gtSample.4),nrow(synSample.4),na.rm=FALSE),
        min(nrow(gtSample.5),nrow(synSample.5),na.rm=FALSE),
        min(nrow(gtSample.6),nrow(synSample.6),na.rm=FALSE),
        min(nrow(gtSample.7),nrow(synSample.7),na.rm=FALSE),
        min(nrow(gtSample.8),nrow(synSample.8),na.rm=FALSE),
        min(nrow(gtSample.9),nrow(synSample.9),na.rm=FALSE),
        min(nrow(gtSample.10),nrow(synSample.10),na.rm=FALSE))

GT.List<-list()
GT.List[[1]]<-gtSample.1
GT.List[[2]]<-gtSample.2
GT.List[[3]]<-gtSample.3
GT.List[[4]]<-gtSample.4
GT.List[[5]]<-gtSample.5
GT.List[[6]]<-gtSample.6
GT.List[[7]]<-gtSample.7
GT.List[[8]]<-gtSample.8
GT.List[[9]]<-gtSample.9
GT.List[[10]]<-gtSample.10

SYN.List<-list()
SYN.List[[1]]<-synSample.1
SYN.List[[2]]<-synSample.2
SYN.List[[3]]<-synSample.3
SYN.List[[4]]<-synSample.4
SYN.List[[5]]<-synSample.5
SYN.List[[6]]<-synSample.6
SYN.List[[7]]<-synSample.7
SYN.List[[8]]<-synSample.8
SYN.List[[9]]<-synSample.9
SYN.List[[10]]<-synSample.10



KLGT<-list()
KLSYN<-list()

for(i in 1:10){
  
kldlist<-list()
# for each GT conduct 10 iterative DL(gt||GT) computations
for(n in 1:10){
#the size variable determines the sample size for both  SYN and GT
  Q3readIn<-sample_n(Q3readInI,size[i])
  rapply( Q3readIn, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
  Q3readIn[is.na(Q3readIn)]<-0
  Q3readIn<-na.omit(Q3readIn)
  Q3readIn[,'age']<-as.numeric(as.character(Q3readIn[,'age']))  
  
  sizeII<-min(size[i],nrow(Q3readIn),na.rm=FALSE) 
  GT.TValuesII<-sample_n(GT.List[[i]], sizeII)
  GT.TValues<-sample_n(Q3readIn, sizeII)
   
  
  den_age<-density(na.omit(GT.TValues$age))
  den_ages<-density(na.omit(GT.TValuesII$age))
  KL<-KLD(den_age$y,den_ages$y) 
  age_m<-cbind(aKLD=c(KL$sum.KLD.px.py)) 
  den_ch<-density(na.omit(GT.TValues$choleratio))
  den_chs<-density((na.omit(GT.TValuesII$choleratio)))
  KL<-KLD(den_ch$y,den_chs$y)
  ch_m<-cbind(cKLD=c(KL$sum.KLD.px.py)) 
  den_bmi<-density(na.omit(GT.TValues$bmi))
  den_bmis<-density((na.omit(GT.TValuesII$bmi)))
  KL<-KLD(den_bmi$y,den_bmis$y)
  bmi_m<-cbind(bKLD=c(KL$sum.KLD.px.py)) 
  den_sbp<-density(na.omit(GT.TValues$sbp))
  den_sbps<-density((na.omit(GT.TValuesII$sbp)))
  KL<-KLD(den_sbp$y,den_sbps$y)
  sbp_m<-cbind(sKLD=c(KL$sum.KLD.px.py)) 
  den_sbps<-density(na.omit(GT.TValues$sbps))
  den_sbpss<-density((na.omit(GT.TValuesII$sbps)))
  KL<-KLD(den_sbps$y,den_sbpss$y)
  sbps_m<-cbind(ssdKLD=c(KL$sum.KLD.px.py))
  
  kldlist[[n]]<-cbind(age_m^2,ch_m^2,bmi_m^2,sbp_m^2,sbps_m^2)
  
}
  ite<-do.call(rbind, kldlist)
  KLGT[[i]]<-colSums(ite)/nrow(ite)


kldlist<-list()
# for each SYN conduct 10 iterative DL(syn||GT) computations
for(n in 1:10){
  Q3readIn<-sample_n(Q3readInI,size[i])
  rapply( Q3readIn, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
  Q3readIn[is.na(Q3readIn)]<-0
  Q3readIn<-na.omit(Q3readIn)
  Q3readIn[,'age']<-as.numeric(as.character(Q3readIn[,'age']))  
  
  sizeII<-min(size[i],nrow(Q3readIn),na.rm=FALSE) 
  SY.TValues<-sample_n(SYN.List[[i]], sizeII)
  GT.TValues<-sample_n(Q3readIn, sizeII)
  
  den_age<-density(na.omit(GT.TValues$age))
  den_ages<-density(na.omit(SY.TValues$age))
  KL<-KLD(den_age$y,den_ages$y) 
  age_m<-cbind(aKLD=c(KL$sum.KLD.px.py)) 
  den_ch<-density(na.omit(GT.TValues$choleratio))
  den_chs<-density((na.omit(SY.TValues$choleratio)))
  KL<-KLD(den_ch$y,den_chs$y)
  ch_m<-cbind(cKLD=c(KL$sum.KLD.px.py)) 
  den_bmi<-density(na.omit(GT.TValues$bmi))
  den_bmis<-density((na.omit(SY.TValues$bmi)))
  KL<-KLD(den_bmi$y,den_bmis$y)
  bmi_m<-cbind(bKLD=c(KL$sum.KLD.px.py)) 
  den_sbp<-density(na.omit(GT.TValues$sbp))
  den_sbps<-density((na.omit(SY.TValues$sbp)))
  KL<-KLD(den_sbp$y,den_sbps$y)
  sbp_m<-cbind(sKLD=c(KL$sum.KLD.px.py)) 
  den_sbps<-density(na.omit(GT.TValues$sbps))
  den_sbpss<-density((na.omit(SY.TValues$sbps)))
  KL<-KLD(den_sbps$y,den_sbpss$y)
  sbps_m<-cbind(ssdKLD=c(KL$sum.KLD.px.py))
  
  kldlist[[n]]<-cbind(age_m^2,ch_m^2,bmi_m^2,sbp_m^2,sbps_m^2)

}
  ite<-do.call(rbind, kldlist)
  KLSYN[[i]]<-colSums(ite)/nrow(ite)
}



table_mI<-as.data.frame(do.call(rbind, KLSYN))
table_mII<-as.data.frame(do.call(rbind, KLGT))



#Experiment V (Table VI) Maximum Mean Discrepancy (MMD),  joint distribution comparison, sample 1000 instances from ground truth and synthetic dataset respectively

for(m in 1:10){
#read gt and tag missing value
GTM.TValues<-sample_n(population,1000)
GTM.TValues[is.na(GTM.TValues)] <- -1
#read syn and tag missing value
SYM.TValues<- sample_n(read.delim(paste("A:/Wangz/R Scripts/Ntr/1000000/",m,"s.txt",sep = "")),1000)
SYM.TValues[is.na(SYM.TValues)] <- -1
#set k from 2 to 4
for(k in 2:4){
  combi<-combn(names(GTM.TValues),k,simplify=FALSE)
  for(i in combi){
    
    x <- data.matrix(GTM.TValues[,c(i)])
    y <-  data.matrix(SYM.TValues[,c(i)])
    mmdo <- kmmd(x,y)
    sink(paste('A:/Wangz/R Scripts/Ntr/1000000/mmd-',m,'-',k,'d.txt',sep = ''), append=TRUE) 
    joint<-c()
    for(j in 1:length(i)){ 
      joint<-paste(joint,c(as.character(i[j])),sep = ":")
    }
    cat(paste(joint,mmdo@H0,sep = ',')) 
    cat("\n")
    sink()
  }
}
}

 

#Experiment VI (Figure 5), PRC and ROC 
library(SuperLearner) 
library(precrec)


trainratio<-0.8
testratio<-(1-trainratio)
results<-list()
for(m in 1:10){
  
  GTM.TValues<-as.data.frame(lapply(sample_n(read.delim(paste("A:/Wangz/R Scripts/Ntr/1000000/",m,"g.txt",sep = "")),100000),as.numeric))
  GTM.TValues[is.na(GTM.TValues)] <- -1
  SYM.TValues<- as.data.frame(lapply(sample_n(read.delim(paste("A:/Wangz/R Scripts/Ntr/1000000/",m,"s.txt",sep = "")),100000),as.numeric))
  SYM.TValues[is.na(SYM.TValues)] <- -1 
  GT.TRAINO<-GTM.TValues[1:(nrow(GTM.TValues)*trainratio),]$strokeha
  GT.TRAINP<-subset(GTM.TValues[1:(nrow(GTM.TValues)*trainratio),], select = -strokeha) 
  SY.TRAINO<-SYM.TValues[1:(nrow(SYM.TValues)*trainratio),]$strokeha 
  SY.TRAINP<-subset(SYM.TValues[1:(nrow(SYM.TValues)*trainratio),], select = -strokeha) 
  GT.TESTO<-GTM.TValues[(nrow(GTM.TValues)*trainratio+1):nrow(GTM.TValues),]$strokeha
  GT.TESTP<-subset(GTM.TValues[(nrow(GTM.TValues)*trainratio+1):nrow(GTM.TValues),], select = -strokeha) 
  
  GT.TESTPNTROW<-length(which(GT.TESTO == 1))
  GT.TESTPNFROW<-length(which(GT.TESTO == 0))
  GT.TESTPNAROW<-GT.TESTPNTROW+GT.TESTPNFROW
  
  #filter equal number of positive cases
  filteredP<-filter(SYM.TValues, strokeha==1)
  SY.TESTPR<-sample_n(filteredP, GT.TESTPNTROW, replace=FALSE)
  filteredN<-filter(SYM.TValues, strokeha==0)
  SY.TESTNR<-sample_n(filteredN, GT.TESTPNFROW, replace=FALSE)
  SY.TESTR<-rbind(SY.TESTPR,SY.TESTNR)
  SY.TESTO<-SY.TESTR$strokeha
  SY.TESTP<-subset(SY.TESTR, select = -strokeha)
  
  
  GT.ALLP<-subset(GTM.TValues[1:nrow(GTM.TValues),], select = -strokeha) 
  GT.ALLO<- GTM.TValues[1:nrow(GTM.TValues),]$strokeha
  
  fit <- SuperLearner(Y = SY.TRAINO,
                      X = SY.TRAINP,
                      family ="gaussian", 
                      cvControl = list(V =10L, stratifyCV = T),
                      SL.library = c("SL.bayesglm"))
  
  SY.PREDOUTCOME<-predict(fit, newdata=SY.TESTP)#, X = SY.TESTP, Y =  SY.TESTO,onlySL = FALSE)
  GT.PREDOUTCOME<-predict(fit, newdata=GT.TESTP)#, X = GT.TESTP, Y =  GT.TESTO,onlySL = FALSE) 
  
  ####PR AND ROC 
  scores <- join_scores(as.numeric(SY.PREDOUTCOME$pred),as.numeric(GT.PREDOUTCOME$pred))
  labels <- join_labels(SY.TESTO, GT.TESTO)
  
  msmdat1 <- mmdata(scores,labels,modnames = c("SYN", "GT"),dsids = c(1, 2) )
  #mspoints <- evalmod(msmdat1)#, mode="basic")
  results[[m]]<-msmdat1
}
datapoints<-evalmod(results[[1]])
#correlation test
rocs<-attr( datapoints ,"grp_avg") [[1]]  
synroc<-rocs[[1]]
gtroc<-rocs[[2]]
#plot(synroc$x,synroc$y_avg)
#plot(gtroc$x,gtroc$y_avg)
df.d1<-data.frame(synroc$x,y=synroc$y_avg)
df.d2<-data.frame(gtroc$x,y=gtroc$y_avg)

prcs<-attr(datapoints  ,"grp_avg") [[2]]  
synpr<-prcs[[1]]
gtpr<-prcs[[2]]
#plot(synpr$x,synpr$y_avg)
#plot(gtpr$x,gtpr$y_avg)
df.d3<-data.frame(synpr$x,y=synpr$y_avg)
df.d4<-data.frame(gtpr$x,y=gtpr$y_avg)
#PRC correlation
#CPRC2<-cor(df.d3$y , df.d4$y )
PRC<-grangertest(df.d3$y,df.d4$y, order=1) 
ROCC<-grangertest(df.d1$y,df.d2$y, order=1)  

#p value
P-value<- c(PRC$`Pr(>F)`[[2]], ROCC$`Pr(>F)`[[2]]) 

##AUC test

aucs<-as.data.frame(attr( datapoints,'aucs'))

#SYN AUC of ROC and PRC
ct1<-aucs$aucs[2]  
ct2<-aucs$aucs[4]
#GT AUC of ROC and PRC
ct3<-aucs$aucs[1] 
ct4<-aucs$aucs[3]
