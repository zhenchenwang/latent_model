library(bnlearn)
library(pcalg)
library(LaplacesDemon)
library(Rgraphviz)
library(ggplot2)
library(gridExtra)
library(pracma)
library(missForest)
library(gRain)
library(cluster)
library(arules)


#ground truth
cvdgt <- read.delim("G:/cvdgt.txt")
#samples
df<-cvdgt[sample(nrow(cvdgt), 122328 ), ]
 

#variables formatting
df$strokeha<-as.factor(df$strokeha)
df$af<-as.factor(df$af)
df$atyantip<-as.factor(df$atyantip)
df$steroid<-as.factor(df$steroid)
df$impot<-as.factor(df$impot)
df$migr<-as.factor(df$migr)
df$ra<-as.factor(df$ra)
df$ckidney<-as.factor(df$ckidney)
df$semi<-as.factor(df$semi)
df$sle<-as.factor(df$sle)
df$treathyp<-as.factor(df$treathyp)
df$type1<-as.factor(df$type1)
df$type2<-as.factor(df$type2)
df$ethr<-as.factor(df$ethr)
df$smoking<-as.factor(df$smoking)
df$fh_cad<-as.factor(df$fh_cad)
df$gender<-as.factor(df$gender)
df$region<-as.factor(df$region)

df$age<-as.numeric(df$age)
df$bmi<-as.numeric(df$bmi)
df$choleratio<-as.numeric(df$choleratio)
df$sbp<-as.numeric(df$sbp)
df$sbps<-as.numeric(df$sbps)

#function that transforms the dataset into integer values (from 0..(n-1))
into.integer <- function(df){
  for (k in 1:ncol(df)){
    df[,k]<-as.integer(df[,k])-1
  }
  return(df)
}

counting.levels<-function(df){
  numLev <- c()
  for (j in 1:ncol(df)){
    numLev <- c(numLev,nlevels(df[,j]))
  }
  return(numLev)
}

#DISCRETE INDEPENDENCE TEST
discretization.equal.intervals <- function(df,n){
  #note: the df.discrete dataset need to be complete (without missing) to compute the    
  #sufficient statistic using the fci algorithm and the levels must be in the range (0..p)
  df.discrete <- df
  levels(df.discrete$gender)<-c(0:3)
  levels(df.discrete$ethr)<-c(0:8)
  #the discrete independence test need a minimum obs to compute the G^2 test, that value  depends on the levels in which the variables are discretize (more levels leads to more observations needed). I choos
  for (i in 1:length(df.discrete)){
    if (is.numeric(df.discrete[,i])){
      if (colnames(df.discrete)[i]=='bmi'){
        lev = n[1]
      } else if (colnames(df.discrete)[i]=='age'){
        lev = n[2]
      } else if (colnames(df.discrete)[i]=='choleratio'){
        lev = n[3]
      } else if (colnames(df.discrete)[i]=='sbp'){
        lev = n[3]
      } else if (colnames(df.discrete)[i]=='sbps'){
        lev = n[3]
      }
      
      df.discrete[,i]<-as.factor(
        arules::discretize(df.discrete[,i],
                           method = 'interval',
                           breaks = lev,
                           labels = c(0:(lev-1)),
                           include.lowest = TRUE,
                           right = TRUE)
      )
    }
  }
  return(df.discrete)
}



 

#latent variable experiment
##################################################################################
#discretizatin of continous variables
df.discrete <- discretization.equal.intervals(df,c(6,4,5,5))

#to run the fci algortithm I have to impute the missing values, I choose to use the random forest to do that
df.discrete.imputed<-missForest(df.discrete)$ximp

#counting the number of levels of each variable
numLev <- counting.levels(df.discrete.imputed)

#levels must starts from 0 and the variables must be integer
df.discrete.imputed<-into.integer(df.discrete.imputed)

#number of resampling
n<-10

#list in which I save the different samples
bootstrap.samples<-list()
rfci.discrete<-list()
suffStat.discrete<-list()

#bootstrap samples
for (i in 1:n){
  index <- sample(1:nrow(df.discrete.imputed),
                  size=nrow(df.discrete.imputed),
                  replace = TRUE)
  bootstrap.samples[[i]]<-df.discrete.imputed[index,]
}

for (i in 1:n){
  suffStat.discrete[[i]] <- list(dm=bootstrap.samples[[i]],
                                 nlev=numLev,
                                 adaptDF = FALSE)
  
  rfci.discrete[[i]] <- rfci(suffStat.discrete[[i]],
                             indepTest = disCItest,
                             alpha=0.9,                         #it has to be thought as tuning parameter
                             skel.method = 'stable',
                             labels=colnames(df.discrete),
                             verbose = TRUE,
                             m.max = 3)
}

#plot the PAG found
for (i in 1:length(rfci.discrete)){
  plot(rfci.discrete[[i]])
} 

df1<-df
########################################################################
for (n in c(1:10)) {


df<-df1[sample(nrow(df1), 100000), ]
gtSample<-df
write.table(df, paste("A:/Wangz/R Scripts/Ntr/",n,"g.txt",sep = ""), sep="\t",row.names = FALSE)
df$strokeha<-as.factor(df$strokeha)
df$af<-as.factor(df$af)
df$atyantip<-as.factor(df$atyantip)
df$steroid<-as.factor(df$steroid)
df$impot<-as.factor(df$impot)
df$migr<-as.factor(df$migr)
df$ra<-as.factor(df$ra)
df$ckidney<-as.factor(df$ckidney)
df$semi<-as.factor(df$semi)
df$sle<-as.factor(df$sle)
df$treathyp<-as.factor(df$treathyp)
df$type1<-as.factor(df$type1)
df$type2<-as.factor(df$type2)
df$ethr<-as.factor(df$ethr)
df$smoking<-as.factor(df$smoking)
df$fh_cad<-as.factor(df$fh_cad)
df$gender<-as.factor(df$gender)
df$region<-as.factor(df$region)

df$age<-as.numeric(df$age)
df$bmi<-as.numeric(df$bmi)
df$choleratio<-as.numeric(df$choleratio)
df$sbp<-as.numeric(df$sbp)
df$sbps<-as.numeric(df$sbps)
#structure without latent
r.without.latent<-structural.em(df,
                                maximize = "hc",
                                fit = "mle",
                                return.all = TRUE,
                                start=NULL,
                                max.iter = 5)
 
#New structure

numVariables <- ncol(df)
variables <- names(df)
#number of latent added
z<-6
#names of the latent added
z.names<-paste("L",1:z,sep="")

#number of latent added
z1<-4
#names of the latent added
z1.names<-paste("M",1:z1,sep="")


#levels of the latent added (for now, I model them as discrete with 2,3,4 levels)
states<-c(2,3,4)

#AMAT obtained using structural.EM enrichment
AMAT.enriched <- cbind(amat(r.without.latent$dag),
                       rep(0,length(variables)),
                       rep(0,length(variables)),
                       rep(0,length(variables)),
                       rep(0,length(variables)),
                       rep(0,length(variables)), 
                       rep(0,length(variables)) 
                       )
colnames(AMAT.enriched)[length(variables)+1]<-z.names[1]
colnames(AMAT.enriched)[length(variables)+2]<-z.names[2]
colnames(AMAT.enriched)[length(variables)+3]<-z.names[3]
colnames(AMAT.enriched)[length(variables)+4]<-z.names[4]
colnames(AMAT.enriched)[length(variables)+5]<-z.names[5]
colnames(AMAT.enriched)[length(variables)+6]<-z.names[6]   


AMAT.enriched<-rbind(AMAT.enriched,
                     rep(0,length(variables)+1),
                     rep(0,length(variables)+2),
                     rep(0,length(variables)+3),
                     rep(0,length(variables)+4),
                     rep(0,length(variables)+5),
                     rep(0,length(variables)+6) 
                     )
rownames(AMAT.enriched)[length(variables)+1]<-z.names[1]
rownames(AMAT.enriched)[length(variables)+2]<-z.names[2]
rownames(AMAT.enriched)[length(variables)+3]<-z.names[3]
rownames(AMAT.enriched)[length(variables)+4]<-z.names[4]
rownames(AMAT.enriched)[length(variables)+5]<-z.names[5]
rownames(AMAT.enriched)[length(variables)+6]<-z.names[6] 

# L1 -> age, af, treathyp
# L2 -> steroid, treathyp
# L3 -> impot, gender
# L4 -> migr, gender, choleratio
# L5 -> strokeha, ckidney, type2, choleratio, sbps
# L6 -> strokeha, ckidney, type2

#L1  
AMAT.enriched["L1","age"]<-1
AMAT.enriched["L1","treathyp"]<-1
AMAT.enriched["L1","af"]<-1

#L2  
AMAT.enriched["L2","steroid"]<-1
AMAT.enriched["L2","treathyp"]<-1

#L3  
AMAT.enriched["L3","impot"]<-1
AMAT.enriched["L3","gender"]<-1

#L4 
AMAT.enriched["L4","migr"]<-1
AMAT.enriched["L4","gender"]<-1
AMAT.enriched["L4","choleratio"]<-1

#L5  
AMAT.enriched["L5","strokeha"]<-1
AMAT.enriched["L5","ckidney"]<-1
AMAT.enriched["L5","type2"]<-1
AMAT.enriched["L5","choleratio"]<-1
AMAT.enriched["L5","sbps"]<-1


#L6
AMAT.enriched["L6","strokeha"]<-1
AMAT.enriched["L6","ckidney"]<-1
AMAT.enriched["L6","type2"]<-1 

 


DAG.enriched<- empty.graph(c(variables,z.names))
amat(DAG.enriched)<-AMAT.enriched
 
#########################################################################
 
##COMPARING DIFFERENT NUMBER OF STATES OF THE DISCRETE LATENT VARIABLE AND THE APPROCH IN WITH I GIVE RANDOM INITIAL VALUES
df.with.latent.clustering<-list()
df.with.latent.random<-list()

for (i in 1:length(states)){
  df.with.latent.clustering[[i]]<-df
  df.with.latent.random[[i]]<-df
  #add a new variables to the dataset and give 'clustering' values
  df.with.latent.clustering[[i]]$L1<-as.factor(
    kmeans(df.discrete.imputed[,c("age","af","treathyp")], 
           centers=states[i], 
           iter.max = 10, 
           nstart = 20)$cluster[sample(nrow(df.discrete.imputed), nrow(df))])
  
  df.with.latent.clustering[[i]]$L2<-as.factor(
    kmeans(df.discrete.imputed[,c("steroid","treathyp")],  
           centers=states[i], 
           iter.max = 10, 
           nstart = 20)$cluster[sample(nrow(df.discrete.imputed), nrow(df))])
  
  df.with.latent.clustering[[i]]$L3<-as.factor(
    kmeans(df.discrete.imputed[,c("impot","gender")],
           centers=states[i], 
           iter.max = 10, 
           nstart = 20)$cluster[sample(nrow(df.discrete.imputed), nrow(df))])
  
  df.with.latent.clustering[[i]]$L4<-as.factor(
    kmeans(df.discrete.imputed[,c("migr", "gender", "choleratio")],
           centers=states[i], 
           iter.max = 10, 
           nstart = 20)$cluster[sample(nrow(df.discrete.imputed), nrow(df))])
  
  df.with.latent.clustering[[i]]$L5<-as.factor(
    kmeans(df.discrete.imputed[,c("strokeha", "ckidney", "type2", "choleratio", "sbps")],
           centers=states[i], 
           iter.max = 10, 
           nstart = 20)$cluster[sample(nrow(df.discrete.imputed), nrow(df))])

  df.with.latent.clustering[[i]]$L6<-as.factor(
    kmeans(df.discrete.imputed[,c("strokeha", "ckidney", "type2")],
           centers=states[i], 
           iter.max = 10, 
           nstart = 20)$cluster[sample(nrow(df.discrete.imputed), nrow(df))])
  
  #add new variables to the dataset and give random values
  df.with.latent.random[[i]]$L1<-as.factor(
    sample(x=0:(states[i]-1),
           replace = TRUE,
           size=nrow(df))
  )
  df.with.latent.random[[i]]$L2<-as.factor(
    sample(x=0:(states[i]-1),
           replace = TRUE,
           size=nrow(df))
  )
  df.with.latent.random[[i]]$L3<-as.factor(
    sample(x=0:(states[i]-1),
           replace = TRUE,
           size = nrow(df))
  )
  df.with.latent.random[[i]]$L4<-as.factor(
    sample(x=0:(states[i]-1),
           replace = TRUE,
           size = nrow(df))
  )
  df.with.latent.random[[i]]$L5<-as.factor(
    sample(x=0:(states[i]-1),
           replace = TRUE,
           size = nrow(df))
  )
  df.with.latent.random[[i]]$L6<-as.factor(
    sample(x=0:(states[i]-1),
           replace = TRUE,
           size = nrow(df))
  )
  
   
  levels(df.with.latent.clustering[[i]]$L1)<-c(0:(states[i]-1))
  levels(df.with.latent.clustering[[i]]$L2)<-c(0:(states[i]-1))
  levels(df.with.latent.clustering[[i]]$L3)<-c(0:(states[i]-1))
  levels(df.with.latent.clustering[[i]]$L4)<-c(0:(states[i]-1))
  levels(df.with.latent.clustering[[i]]$L5)<-c(0:(states[i]-1))
  levels(df.with.latent.clustering[[i]]$L6)<-c(0:(states[i]-1))
  
}

 
r.random<-list()

df.synthetic.random<-list()

for (j in 1:1){
    
  r.random[[j]]<-structural.em(df.with.latent.clustering[[j]],
                               maximize = "hc", 
                               maximize.args = list(score="aic-cg"),
                               fit =  "mle",
                               fit.args = list(replace.unidentifiable =TRUE),
                               #replace.unidentifiable =TRUE ,
                               return.all = TRUE,
                               start=DAG.enriched,
                               max.iter = 5,
                               debug = FALSE )
 
  xtime<-as.POSIXct( Sys.time() )
  as.integer( xtime )

  synSample<-rbn(r.random[[j]]$fitted,
                                 n=nrow(df))
  
  #formatting and biological checks
  synSample[,'bmi']<-round(synSample[,'bmi'],1)
  synSample[,'choleratio']<-round(synSample[,'choleratio'],1)
  synSample[,'sbp']<-round(synSample[,'sbp'],0)
  synSample[,'sbps']<- abs(round(synSample[,'sbps'],2))
  synSample[,'age']<-round(synSample[,'age'],0)
  
  synSample<-synSample[!(synSample$gender=="F" & synSample$`impot`==1),]
 
  
  write.table(synSample, paste("A:/Wangz/R Scripts/Ntr/",n,"s.txt",sep = ""), sep="\t",row.names = FALSE) 
  
}
}

#FUNCTION  for every variable return the probabiity (according to the graphS obtained)
#to be linked to a specific other variable of the dataset

#input: 
#-list of models obtained by the fci function FROM DIFFERENT bootstrap samples
#-variables names 

find_level_confidence <- function(list_fciAlgo, var.names){
  output<-data.frame()
  for (i in 1:length(var.names)){ #fixed variable X1
    for (j in 1:length(var.names)){ #all other variables X2 (I'm checking the match beetween the two)
      if (var.names[[i]]!=var.names[[j]]){
        conf <- 0
        for (k in 1:size(list_fciAlgo,2)){ #models learnt on bootstrap samples
          amat <- list_fciAlgo[[k]]@amat
          if(amat[i,j]==2 && amat[j,i]==2){ #bidirected edge 
            conf <- conf + 1
          } else if ((amat[i,j]==2 && amat[j,i]==1)||(amat[i,j]==1 && amat[j,i]==2)){ # X1 o--> X2 or X2 0--> X1
            conf <- conf + 0.5
          } else if (amat[i,j]==1 && amat[j,i]==1){ #o--o
            conf <- conf + 0.33
          }
        }
        partial <- data.frame(var.names[i],var.names[j],conf/size(list_fciAlgo,2))
        output<-rbind(output,partial)
      }
    }
  }
  colnames(output)[1]<-"variable1"
  colnames(output)[2]<-"variable2"
  colnames(output)[3]<-"confidenceLevel"
  return(output)
}

output<-find_level_confidence(rfci.discrete,variables)
output


#which are the most probable latent variables beetween the observed ?
plot.level.confidence<-function(){
  p<-list()
  for (i in 1:numVariables){
    treshold <- 0
    if (i==2 || i==5 || i==6 || i==7){ #for continuous variable I consider a level of threshold less than the one for discrete one
      treshold<-0.7
    }else{
      treshold<-0.9
    }
    
    
    p[[i]]<-ggplot(subset(output,subset=(variable1==variables[i]),select=-variable1), aes(x=variable2, y=confidenceLevel))+ 
      geom_bar(aes(),   # fill depends on cond2
               stat="identity",
               colour="black",
               fill="chartreuse3", # Black outline for all
               position=position_dodge())+
      ggtitle(paste("Variable 1: ",variables[i]))+
      ylim(0, 1)+
      geom_hline(yintercept=treshold, linetype="dashed", color = "red")
    
  }
  grid.arrange(p[[2]],p[[5]],p[[6]],p[[7]])
  grid.arrange(p[[1]],p[[3]],p[[4]],p[[8]])
  grid.arrange(p[[9]],p[[10]],p[[11]])
  grid.arrange(p[[12]],p[[13]])
}
