
#### code to produce: main manuscript data figure + table; 
#### appendix cook distance, appendix time series fish predators
##########################################################################
  
  setwd("C:/Users/pdvd/Online for git/RAM_forage_piscivores/Processing/")
  source("Process_timeseries_ER.R")
  library(nlme)
  library(bbmle)  

# make a time-series plot for each predator
# include smoother for illustration
# select the years with the largest decline in biomass 

  setwd("C:/Users/pdvd/Online for git/RAM_forage_piscivores/Output/")
  pdf("Timeseries_predators.pdf",width=8.3,height=11.27)
  par(mfrow=c(4,4))

# get predators
  stid<-sort(unique(finaldata$stocklong))
 
# select the number of years to determine the largest decline in biomass
  sh <- 4 #  1 + 4  = minimum nb of years to determine largest decline in stock
  ma<-14  #  1 + 14 = maximum nb of years to determine largest decline in stock 

# get correct names for plotting
  plname<-c("NA","NA","NA","Arrowtooth flounder \n Pacific Coast",
  "Atlantic cod \n Baltic areas 25-32",
  "Atlantic cod \n Iceland","NA","NA",
  "Atlantic cod \n North Sea",
  "Atlantic cod \n West of Scotland",
  "*** Southern blue whiting \n Chile",
  "** Blue whiting \n North East Atlantic",
  "Deep-water cape hake \n South Africa",
  "Haddock \n Iceland",
  "Haddock \n ICES IIIa and North Sea",
  "Haddock \n West of Scotland",
  "Hake \n Northeast Atlantic South",
  "** Kingklip \n South Africa",
  "Megrim \n ICES VIIIc-IXa","NA",
  "Pacific cod \n W. Coast of Vancouver I.",
  "Pacific hake \n Pacific coast",
  "Petrale sole \n Pacific coast",
  "Pink cusk-eel \n Chile",
  "Pollock \n ICES IIIa, VI and North Sea",
  "Pollock \n Iceland",
  "South hake \n Chile",
  "South pacific hake \n Chile",
  "Spotted spiny dogfish \n Pacific coast",
  "Whiting \n ICES VIa",
  "Yelloweye rockfish \n Pacific coast",
  "Yellownose skate \n Chile")

# get data per predator on largest decline in time series
  stockall<-as.data.frame(matrix(data=NA,ncol=5,nrow=length(stid)))
  
  for (j in 1:length(stid)){
    stock<-subset(finaldata,finaldata$stocklong == stid[j]) ### select the predator 1:32
    
    # only run if data for > 20 
    if(length(stock$year) > 20){
    
      # add smoother for illustration
      TBsmooth  <- smooth.spline(stock$year, stock$bio, cv = TRUE, nknots=round(length(stock$year)/3))
      Yfine <-seq(min(stock$year),max(stock$year),1)
      biosmooth<-predict(TBsmooth, Yfine) 
      stock<-cbind(stock,biosmooth$y)
      colnames(stock)[18]<-"biosmooth"
        
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
      
      # plot all data
      plot(stock$bio~stock$year,main=plname[j],xlab="Year",ylab="Total biomass",cex.lab=1,cex.main=1)
      points(stock$bio[xrow:(xrow+ycol+(sh-1))]~stock$year[xrow:(xrow+ycol+(sh-1))],col="red",pch=16)
      lines(stock$biosmooth~stock$year)
      
      # get the data per predator on largest decline in time series
      stockall[j,1]<-xrow
      stockall[j,2]<-ycol+(sh-1)
      stockall[j,3]<-stock$year[xrow]
      stockall[j,4]<-as.character(stock$stocklong[1])
      stockall[j,5]<-as.character(stock$stockid_prey[1])
    }}
  
  dev.off()
  
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

# fit the interaction model
  mod  <- lm(log(ratiodif)~(ERpred)+(ERprey),data=stockall)
  mod1 <- lm(log(ratiodif)~(ERpred)*(ERprey),data=stockall)
  mod2 <- lm(log(ratiodif)~(ERpred),data=stockall)
  #AICtab(mod,mod1,mod2)
  
# cooks distance plot
  pdf("Cook distance plot.pdf",width=5,height=4.5)
  barplot(sort(cooks.distance(mod1)),xaxt="n",xlab="Data points (stocks)",ylab="Cook's distance",
          ylim=c(0,1.2),yaxt="n")
  axis(2,c(0,0.6,1.2),las=1)
  box()
  dev.off()

# remove blue whiting Chile -- outlier
  stockall <- stockall[-6,] 
  
# check the three models (table main document)
  mod  <- lm(log(ratiodif)~(ERpred)+(ERprey),data=stockall)
  mod1 <- lm(log(ratiodif)~(ERpred)*(ERprey),data=stockall)
  mod2 <- lm(log(ratiodif)~(ERpred),data=stockall)
  #AICtab(mod,mod1,mod2)

# make figure contour
  pdf("Contour plot.pdf",width=8,height=5)
  ERpred<-seq(0.023,0.46,0.001)
  ERprey<-seq(0.078,0.268,0.001)
  tabz<-merge(ERpred,ERprey)
  pred <-data.frame(tabz)
  colnames(pred)<-c("ERpred","ERprey")
  pred$fit <-predict(mod1,newdata=pred,type="response")
  pred$Biomass_ratio <- exp(pred$fit)
  
  library(ggplot2)
  cont<-ggplot(pred) + 
    aes(x = ERpred, y = ERprey, z = Biomass_ratio, fill = Biomass_ratio) + 
    geom_tile() +
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) +
    stat_contour(breaks=c(0.15,0.2,0.3,0.4,0.6),col="white")+
    scale_fill_gradientn(colours=c("#2066ac", "#cde3ef","#f9c1a5","#de725a","#d05647","#b2182b"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            axis.text=element_text(size=15),axis.title=element_text(size=15)) +
            labs( x = "Piscivore fishing mortality",y = "Forage-fish fishing mortality")
  
  new <- data.frame(test1= stockall$ERpred, test2= stockall$ERprey)
  
  cont+ geom_point(data=new, aes(x=test1, y=test2 ,fill=NULL, z=NULL),colour="black") 

  dev.off()