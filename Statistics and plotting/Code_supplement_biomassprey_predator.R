
#### code to produce appendix figure biomass prey versus predator
##########################################################################
  setwd("C:/Users/pdvd/Online for git/RAM_forage_piscivores/Processing/")
  source("Process_timeseries_ER.R")
  library(latex2exp)

# get predators
  stid<-sort(unique(finaldata$stocklong))

# select the number of years to determine the largest decline in biomass
  sh <- 4 #  1 + 4  = minimum nb of years to determine largest decline in stock
  ma<-14  #  1 + 14 = maximum nb of years to determine largest decline in stock 

  
  # get data per predator on largest decline in time series
  stockall<-as.data.frame(matrix(data=NA,ncol=5,nrow=length(stid)))
  
  for (j in 1:length(stid)){
    stock<-subset(finaldata,finaldata$stocklong == stid[j]) ### select the predator 1:32
    
    # only run if data for > 20 
    if(length(stock$year) > 20){
      
      stockdata<-as.data.frame(matrix(data=100,nrow=100,ncol=18))
      for(p in sh:ma){
        nb <- p # number of years +1 to determine the difference
        lentime<-length(stock$year)-nb # lentime describes how many times the biomass ratio can be determined
        for (i in 1:lentime){
          stockdata[i,(p-(sh-1))]<-stock$bio[i+nb]/stock$bio[i]
        }
      }
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
  
# estimate ER, biomass difference predator and biomass predator and prey
  stck<-unique(stockall[,1])
  for (j in 1:nrow(stockall)){
    stock<-subset(finaldata,finaldata$stocklong == stid[stck[j]]) ### unique stock
    yr<-stockall$start_dec[j]
    yrend<-stockall$start_dec[j]+stockall$lenght_dec[j]
    stockall$ratiodif[j]<-stock$bio[yrend]/stock$bio[yr]
    stockall$ERpred[j]<-mean(stock$ER[(yr):(yrend)])
    stockall$ERprey[j]<-mean(stock$ER_prey[(yr):(yrend)])
    stockall$Biopred[j]<-mean(stock$bio[(yr):(yrend)])   # get mean of biomass predator (TB and when unavailable SSB)
    stockall$Bioprey[j]<-mean(stock$TB_prey[(yr):(yrend)]) # get mean of biomass predator (TB and when unavailable SSB)
    stockall$Bioprey[j] <- ifelse(is.na(stockall$Bioprey[j]),mean(stock$SSB_prey[(yr):(yrend)]),stockall$Bioprey[j])
  }

# plot the average biomass predator and prey over the time series
  setwd("C:/Users/pdvd/Online for git/RAM_forage_piscivores/Output")
  pdf("Biomass_predator_prey.pdf",width=5,height=4.5)
  plot(log10(stockall$Biopred)~log10(stockall$Bioprey),xlim = c(3,7.5),ylim=c(3,7.5),xaxt="n",yaxt="n",
       xlab="Average biomass prey (metric tons)", ylab="Average biomass predator (metric tons)")
  axis(1,c(3,5,7),c(TeX("10$^3$"),TeX("10$^5$"),TeX("10$^7$")))
  axis(2,c(3,5,7),c(TeX("10$^3$"),TeX("10$^5$"),TeX("10$^7$")),las=1)
  lines(x=c(0,10),y=c(0,10),lwd=1, lty=3)  

  dev.off()
