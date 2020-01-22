
#### code to produce supplement S7
##########################################################################
  setwd("C:/Users/pdvd/Online for git/RAM_forage_piscivores/Processing/")
  source("Process_timeseries_ER.R")
  library(nlme)
  library(bbmle)  

# get predators
  stid<-sort(unique(finaldata$stocklong))

# select the number of years to determine the largest decline in biomass
  min_step <- c(5:7,5:8,5:9,5:10,5:11,5:12,5:13,5:14,5:14,5:14,5:14,5:14)
  max_step <- c(rep(8,3),rep(9,4),rep(10,5),rep(11,6),rep(12,7),rep(13,8),rep(14,9),rep(15,10),rep(16,10),rep(17,10),rep(18,10),rep(19,10))
  idx <- length(min_step) # equal to length(max_step)

# create datatable S7
  Sup7tab <- as.data.frame(matrix(data=NA, nrow = 92, ncol= 5))
  Sup7tab[,1] <- min_step
  Sup7tab[,2] <- max_step
  colnames(Sup7tab) <- c("min_lengths","max_lengths","Model","R2_mult","R2_adj")
  
# estimate for each predator stock and time period
# the largest decline in predator biomass and get ER
# during the period of largest decline and fit the model
  for (q in 1:idx){
    sh <- min_step[q]-1  #  sh  = minimum nb of years to determine largest decline in stock
    ma <- max_step[q]-1  #  ma  = maximum nb of years to determine largest decline in stock 
    stockall<-as.data.frame(matrix(data=NA,ncol=5,nrow=32))
    
    for (j in 1:32){
      stock<-subset(finaldata,finaldata$stocklong == stid[j]) ### unique stock 1:32
      
      # only run if data for > 20 
      if(length(stock$year) > 20){
        stockdata<-as.data.frame(matrix(data=100,nrow=100,ncol=18))
        
        for(p in sh:ma){
          nb<-p # number of years +1 to determine the difference
          lentime<-length(stock$year)-nb # lentime describes how many times the biomass ratio is determined
          for (i in 1:lentime){
            stockdata[i,(p-(sh-1))]<-stock$bio[i+nb]/stock$bio[i]
          }}
        nb<-which(stockdata == min(stockdata)) [1]
        cl<- seq(100,100*18,length.out = 18)
        ycol<-(which(cl-nb > 0)[1])
        if (ycol > 1){
          xrow<-abs(cl-nb)[ycol-1]
        }
        if (ycol == 1){
          xrow<-nb
        }
        
        # get the data per predator on largest decline in time series
        stockall[j,1]<-xrow
        stockall[j,2]<-ycol+(sh-1)
        stockall[j,3]<-stock$year[xrow]
        stockall[j,4]<-as.character(stock$stocklong[1])
        stockall[j,5]<-as.character(stock$stockid_prey[1])
      }}
    
    # data for all stocks, remove two predators that are increasing  
    stockall<-cbind(c(1:32),stockall)
    stockall<-stockall[complete.cases(stockall[,2:3]),]
    stockall<-stockall[-7,] ### remove blue whiting north atlantic (increasing)
    stockall<-stockall[-12,] ### remove kinglip south africa (increasing)
    colnames(stockall)<- c("nb","start_dec","lenght_dec","year_dec","predator","prey")
    
    # estimate ER, biomass difference predator
    stck<-unique(stockall[,1])
    for (j in 1:nrow(stockall)){
      stock<-subset(finaldata,finaldata$stocklong == stid[stck[j]]) ### unique stock
      yr<-stockall$start_dec[j]
      yrend<-stockall$start_dec[j]+stockall$lenght_dec[j]
      stockall$ratiodif[j]<-stock$bio[yrend]/stock$bio[yr]
      stockall$ERpred[j]<-mean(stock$ER[(yr):(yrend)])
      stockall$ERprey[j]<-mean(stock$ER_prey[(yr):(yrend)])
    }
  
    # remove blue whiting Chile -- outlier
    stockall <- stockall[-6,] 
  
    Int    <-lm(log(ratiodif)~(ERpred)+(ERprey),data=stockall)
    Full  <-(lm(log(ratiodif)~(ERpred)*(ERprey),data=stockall))
    Simple <-lm(log(ratiodif)~(ERpred),data=stockall)
    
    # use r.squared or adj.r.squared
    dat <- matrix(data=c(summary(Int)$r.squared,summary(Int)$adj.r.squared,
                         summary(Full)$r.squared,summary(Full)$adj.r.squared,
                         summary(Simple)$r.squared,summary(Simple)$adj.r.squared),ncol=3,nrow=2)
    colnames(dat) <- c("Int","Full","Simple")
    
    datAIC <-  as.data.frame(AIC(Int,Full,Simple)) 
    
    bestmod <-   ifelse((((datAIC[2,2]+2) < datAIC[1,2]) & ((datAIC[2,2]+2) < datAIC[3,2])), "Full", ifelse(((datAIC[1,2]+2) < datAIC[3,2]), "Int", "Simple"))
    Sup7tab[q,3] <- bestmod
    Sup7tab[q,4] <- dat[1,paste(bestmod)]
    Sup7tab[q,5] <- dat[2,paste(bestmod)]
  }
    
  setwd("C:/Users/pdvd/Online for git/RAM_forage_piscivores/Output/")
  write.csv(Sup7tab,file="supplement_table_S7.csv",row.names = F)
  
