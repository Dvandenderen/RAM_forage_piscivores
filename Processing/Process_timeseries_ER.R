
##### processing predator - prey timeseries based on area overlap ##### 
#######################################################################

# select data of predator and prey with estimated overlap
# data has been pre-processed to identify "predators" and "prey" (see method section)
  setwd("C:/Users/pdvd/Online for git/RAM_forage_piscivores/Data")
  ppoverlap <- readRDS("overlap_areas_predator_prey.rds")

# select all areas where more than 50 percent is overlapping (both predator with prey and prey with predator)
  ppoverlap<-subset(ppoverlap,ppoverlap[,13]>0.5)
  ppoverlap<-subset(ppoverlap,ppoverlap[,14]>0.5)

# now select first all predators with data available on ER (exploitation rate--catch/biomass) 
# or F (fishing mortality) and TB or SSB
  stocks<-read.table("timeseries_vals.txt",sep="\t",header=T)

  idpred<-unique(ppoverlap$stockid_pred)
  mysum <- function(x) if (all(is.na(x))) NA else sum(x,na.rm=T)
  datfin<-c()
  for (i in 1:length(idpred)){
    fishP <- subset(stocks,stocks$stockid == as.character(idpred[i]))
    fishP <- subset (fishP,fishP$TB >= 0 | fishP$SSB >= 0) ## select all years with TB or SSB or both
    fishP <- subset (fishP,fishP$F >= 0 | fishP$ER >= 0) ## select all years with F or ER or both
    fishP$bio<-apply(fishP[5:6],1,mysum) 
    fishP$exp<-apply(fishP[11:12],1,mysum)
    fishP <- subset (fishP, fishP$exp >= 0 &  fishP$bio >= 0) ## select all years with F/ER and TB/SSB
    fishP <- fishP[,c(1,2,3,4,5,6,11,12,15)]  
    datfin<- rbind(datfin,fishP)
  }

# check now for each predator, the prey (can be more than 1)
  idpred<-unique(ppoverlap$stockid_pred)
  datfinprey<-c()
  nbPrey<-c()
  for (i in 1:length(idpred)){
    fishP <- subset(ppoverlap,ppoverlap$stockid_pred == as.character(idpred[i]))
    idprey <- unique(fishP$stockid_prey)
    for (j in 1:length(idprey)){
      fishprey <- subset(stocks,stocks$stockid == as.character(idprey[j]))
      fishprey <- subset (fishprey,fishprey$TB >= 0 | fishprey$SSB >= 0) ## select all years with TB or SSB or both
      fishprey <- subset (fishprey,fishprey$F >= 0 | fishprey$ER >= 0) ## select all years with F or ER or both
      fishprey$bio<-apply(fishprey[5:6],1,mysum) 
      fishprey$exp<-apply(fishprey[11:12],1,mysum)
      fishprey <- subset (fishprey, fishprey$exp >= 0 &  fishprey$bio >= 0) ## select all years with F/ER and TB/SSB
      fishprey <- fishprey[,c(1,2,3,4,5,6,11,12,15)]  
      fishprey$idpredator<-idpred[i]
      datfinprey<-rbind(datfinprey,fishprey)
    } 
    nbPrey <- c(nbPrey,length(idprey))  # checks which predators have multiple prey
  }

# select all prey populations that are the single prey for the predator 
  single<-data.frame(idpred,nbPrey)
  datfinprey<-cbind(datfinprey,single[match(datfinprey$idpredator,single$idpred),c(2)])
  colnames(datfinprey)[11]<-"nb_prey"
  preypop<-subset(datfinprey,nb_prey ==1) ## data file with all single prey populations
  colnames(preypop)[9]<-"Fmsy"

# now combine multiple preys (one by one predator) following nbPrey
# nbPrey: 1 1 2 1 1 1 1 1 2 1 5 5 5 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 2 2 2 2 2

  # idpred[3] -> COD4TVn
  subPred<-subset(datfinprey,idpredator == as.character(idpred[3]))  # complete match between the two prey populations
  nb<-length(unique(subPred$year))
  assessid<-rep("merged",nb)
  stockid<-rep(paste(unique(subPred$stockid)[1],unique(subPred$stockid)[2],sep="_"),nb)
  stocklong<-rep("merged",nb)
  year<-unique(subPred$year)
  TB<-subPred$TB[1:29]+subPred$TB[30:58]
  SSB<-subPred$SSB[1:29]+subPred$SSB[30:58]
  F<-(subPred$F[1:29]*subPred$TB[1:29]+subPred$F[30:58]*subPred$TB[30:58])/TB
  ER<-(subPred$ER[1:29]*subPred$TB[1:29]+subPred$ER[30:58]*subPred$TB[30:58])/TB
  Fmsy<-subPred$F.Fmsy[1:29]+subPred$F.Fmsy[30:58]
  idpredator<-rep(subPred$idpredator[1],nb)
  nb_prey<-rep(subPred$nb_prey[1],nb)
  prey3<-data.frame(assessid,stockid,stocklong,year,TB,SSB,F,ER,Fmsy,idpredator,nb_prey)
  preypop<-rbind(preypop,prey3)

  # idpred[9] -> CODBA2532
  subPred<-subset(datfinprey,idpredator == as.character(idpred[9]))  ### complete match between the two prey populations
  nb<-length(unique(subPred$year))
  assessid<-rep("merged",nb)
  stockid<-rep(paste(unique(subPred$stockid)[1],unique(subPred$stockid)[2],sep="_"),nb)
  stocklong<-rep("merged",nb)
  year<-unique(subPred$year)
  TB<-subPred$TB[1:nb]+subPred$TB[(nb+1):(2*nb)]
  SSB<-subPred$SSB[1:nb]+subPred$SSB[(nb+1):(2*nb)]
  F<-(subPred$F[1:nb]*subPred$TB[1:nb]+subPred$F[(nb+1):(2*nb)]*subPred$TB[(nb+1):(2*nb)])/TB
  ER<-(subPred$ER[1:nb]*subPred$TB[1:nb]+subPred$ER[(nb+1):(2*nb)]*subPred$TB[(nb+1):(2*nb)])/TB
  Fmsy<-(subPred$F.Fmsy[1:nb]*subPred$TB[1:nb]+subPred$F.Fmsy[(nb+1):(2*nb)]*subPred$TB[(nb+1):(2*nb)])/TB
  idpredator<-rep(subPred$idpredator[1],nb)
  nb_prey<-rep(subPred$nb_prey[1],nb)
  prey9<-data.frame(assessid,stockid,stocklong,year,TB,SSB,F,ER,Fmsy,idpredator,nb_prey)
  preypop<-rbind(preypop,prey9)

  # idpred[11] -> HADNS-IIIa idpred[12] -> POLLNS-VI-IIIa idpred[13] -> CODNS (same preys for three predators)
  subPred<-subset(datfinprey,idpredator == as.character(idpred[11]))  # five prey populations
  TB_y<-as.data.frame(matrix(data=NA,nrow=length(unique(subPred$year)),ncol=1))
  rownames(TB_y)<-sort(unique(subPred$year))
  stoc1<-subset(subPred,subPred$stockid == unique(subPred$stockid)[1])
  TB_y<-cbind(TB_y, stoc1[match(rownames(TB_y), stoc1$year),c("TB","SSB","F","ER","F.Fmsy")])
  stoc2<-subset(subPred,subPred$stockid == unique(subPred$stockid)[2])
  TB_y<-cbind(TB_y, stoc2[match(rownames(TB_y), stoc2$year),c("TB","SSB","F","ER","F.Fmsy")])
  stoc3<-subset(subPred,subPred$stockid == unique(subPred$stockid)[3])
  TB_y<-cbind(TB_y, stoc3[match(rownames(TB_y), stoc3$year),c("TB","SSB","F","ER","F.Fmsy")])
  stoc4<-subset(subPred,subPred$stockid == unique(subPred$stockid)[4])
  TB_y<-cbind(TB_y, stoc4[match(rownames(TB_y), stoc4$year),c("TB","SSB","F","ER","F.Fmsy")])
  stoc5<-subset(subPred,subPred$stockid == unique(subPred$stockid)[5])
  TB_y<-cbind(TB_y, stoc5[match(rownames(TB_y), stoc5$year),c("TB","SSB","F","ER","F.Fmsy")])
  TB_y[is.na(TB_y)] <- 0

  nb<-length(unique(subPred$year))
  assessid<-rep("merged",nb)
  stockid<-rep(paste(unique(subPred$stockid)[1],unique(subPred$stockid)[2],unique(subPred$stockid)[3],unique(subPred$stockid)[4],unique(subPred$stockid)[5],sep="_"),nb)
  stocklong<-rep("merged_use_F_and_TB",nb)
  year<-sort(unique(subPred$year))
  TB<-rowSums(TB_y[,c(2,7,12,17,22)],na.rm=T)
  SSB<-rowSums(TB_y[,c(3,8,13,18,23)],na.rm=T)
  F<-(TB_y[,4]*TB_y[,2]+TB_y[,9]*TB_y[,7]+TB_y[,14]*TB_y[,12]+TB_y[,19]*TB_y[,17]+TB_y[,24]*TB_y[,22])/TB
  ER<-(TB_y[,5]*TB_y[,2]+TB_y[,10]*TB_y[,7]+TB_y[,15]*TB_y[,12]+TB_y[,20]*TB_y[,17]+TB_y[,25]*TB_y[,22])/TB
  Fmsy<-(TB_y[,6]*TB_y[,2]+TB_y[,11]*TB_y[,7]+TB_y[,16]*TB_y[,12]+TB_y[,21]*TB_y[,17]+TB_y[,26]*TB_y[,22])/TB
  idpredator<-rep(subPred$idpredator[1],nb)
  nb_prey<-rep(subPred$nb_prey[1],nb)

  prey11<-data.frame(assessid,stockid,stocklong,year,TB,SSB,F,ER,Fmsy,idpredator,nb_prey)
  preypop<-rbind(preypop,prey11)
  
  subPred<-subset(datfinprey,idpredator == as.character(idpred[12])) 
  prey12<-prey11
  prey12$idpredator<-subPred$idpredator[1]
  preypop<-rbind(preypop,prey12)
  
  subPred<-subset(datfinprey,idpredator == as.character(idpred[13]))  
  prey13<-prey11
  prey13$idpredator<-subPred$idpredator[1]
  preypop<-rbind(preypop,prey13)

  # idpred[14] -> CODICE idpred[15] -> POLLIEG idpred[16] -> HADICE (same prey for three predators) 
  # capelin is removed, unrealistic high ER, see method section)
  subPred<-subset(datfinprey,idpredator == as.character(idpred[14])) 
  subPred<-subset(subPred,subPred$stockid =="HERRIsum")   # remove capelin
  nb<-length(unique(subPred$year))
  assessid<-rep("merged",nb)
  stockid<-rep(paste(unique(subPred$stockid)[1],unique(subPred$stockid)[2],sep="_"),nb)
  stocklong<-rep("herring_only_capelin_removed",nb)
  year<-sort(unique(subPred$year))
  TB<-subPred$TB
  SSB<-subPred$SSB
  F<-subPred$F
  ER<-subPred$ER
  Fmsy<-subPred$F.Fmsy
  idpredator<-rep(subPred$idpredator[1],nb)
  nb_prey<-rep(subPred$nb_prey[1],nb)

  prey14<-data.frame(assessid,stockid,stocklong,year,TB,SSB,F,ER,Fmsy,idpredator,nb_prey)
  preypop<-rbind(preypop,prey14)

  subPred<-subset(datfinprey,idpredator == as.character(idpred[15]))
  prey15<-prey14
  prey15$idpredator<-subPred$idpredator[1]
  preypop<-rbind(preypop,prey15)

  subPred<-subset(datfinprey,idpredator == as.character(idpred[16]))
  prey16<-prey14
  prey16$idpredator<-subPred$idpredator[1]
  preypop<-rbind(preypop,prey16)

  # idpred[17] -> WHITVIa idpred[18] -> CODVIa  idpred[19] -> HADVIa (same prey for three predators) 
  # one prey stock is included in the other, hence only combined stock is included HERRVIaVIIbc
  subPred<-subset(datfinprey,idpredator == as.character(idpred[17]))

  nb<-length(unique(subPred$year))
  assessid<-rep("merged",nb)
  stockid<-rep(unique(subPred$stockid)[2],sep="_",nb)
  stocklong<-rep("oneHerring",nb)
  year<-unique(subPred$year)
  TB<-subPred$TB[(nb+1):(2*nb)] # only HERRVIaVIIbc is included 
  SSB<-subPred$SSB[(nb+1):(2*nb)] # only HERRVIaVIIbc is included 
  F<-subPred$F[(nb+1):(2*nb)] # only HERRVIaVIIbc is included 
  ER<-subPred$ER[(nb+1):(2*nb)] # only HERRVIaVIIbc is included 
  Fmsy<-subPred$F.Fmsy[(nb+1):(2*nb)] # only HERRVIaVIIbc is included 
  idpredator<-rep(subPred$idpredator[1],nb)
  nb_prey<-rep(subPred$nb_prey[1],nb)

  prey17<-data.frame(assessid,stockid,stocklong,year,TB,SSB,F,ER,Fmsy,idpredator,nb_prey)
  preypop<-rbind(preypop,prey17)
  
  subPred<-subset(datfinprey,idpredator == as.character(idpred[18])) 
  prey18<-prey17
  prey18$idpredator<-subPred$idpredator[1]
  preypop<-rbind(preypop,prey18)
  
  subPred<-subset(datfinprey,idpredator == as.character(idpred[19])) 
  prey19<-prey17
  prey19$idpredator<-subPred$idpredator[1]
  preypop<-rbind(preypop,prey19)

  # idpred[31] -> DEEPCHAKESA idpred[32] -> KINGKLIPSA (same prey for two predators)
  subPred<-subset(datfinprey,idpredator == as.character(idpred[31])) # no data on F for both species
  nb<-length(unique(subPred$year))
  assessid<-rep("merged",nb)
  stockid<-rep(paste(unique(subPred$stockid)[1],unique(subPred$stockid)[2],sep="_"),nb)
  stocklong<-rep("merged_nodata_on_F_bothspecies",nb)
  year<-unique(subPred$year)
  TB<-subPred$TB[1:nb]+subPred$TB[(nb+1):(2*nb)]
  SSB<-subPred$SSB[1:nb]+subPred$SSB[(nb+1):(2*nb)]
  F<-rep(NA,nb)
  ER<-(subPred$ER[1:nb]*subPred$TB[1:nb]+subPred$ER[(nb+1):(2*nb)]*subPred$TB[(nb+1):(2*nb)])/TB
  Fmsy<-rep(NA,nb)
  idpredator<-rep(subPred$idpredator[1],nb)
  nb_prey<-rep(subPred$nb_prey[1],nb)

  prey31<-data.frame(assessid,stockid,stocklong,year,TB,SSB,F,ER,Fmsy,idpredator,nb_prey)
  preypop<-rbind(preypop,prey31)
  
  subPred<-subset(datfinprey,idpredator == as.character(idpred[32]))  
  prey32<-prey31
  prey32$idpredator<-subPred$idpredator[1]
  preypop<-rbind(preypop,prey32)

  # idpred[40] -> PHAKEPCOAST idpred[41] -> ARFLOUNDPCOAST idpred[42]-> YEYEROCKPCOAST 
  # idpred[43] -> PSOLEPCOAST idpred[44] -> SPSDOGPCOAST (same preys for five predators)
  subPred<-subset(datfinprey,idpredator == as.character(idpred[40])) # no data on F both species and TB one species
  nb<-length(unique(subPred$year))
  assessid<-rep("merged",nb)
  stockid<-rep(paste(unique(subPred$stockid)[1],unique(subPred$stockid)[2],sep="_"),nb)
  stocklong<-rep("merged_nodata_on_TB_onespecies_and_F_bothspecies",nb)
  year<-unique(subPred$year)
  TB<-c(subPred$TB[1:52],subPred$TB[(nb+1):(106)],NA,NA)
  SSB<-c(subPred$SSB[1:52],subPred$SSB[(53):(nb-2)]+subPred$SSB[(nb+1):(106)],subPred$SSB[(nb-1):nb])
  F<-rep(NA,nb)
  ER<-c(subPred$ER[1:52],(subPred$ER[(53):(nb-2)]*subPred$SSB[(53):(nb-2)]+subPred$ER[(nb+1):(106)]*subPred$SSB[(nb+1):(106)])/SSB[53:78],subPred$ER[(nb-1):nb])
  Fmsy<-rep(NA,nb)
  idpredator<-rep(subPred$idpredator[1],nb)
  nb_prey<-rep(subPred$nb_prey[1],nb)
  
  prey40<-data.frame(assessid,stockid,stocklong,year,TB,SSB,F,ER,Fmsy,idpredator,nb_prey)
  preypop<-rbind(preypop,prey40)
  
  subPred<-subset(datfinprey,idpredator == as.character(idpred[41]))  
  prey41<-prey40
  prey41$idpredator<-subPred$idpredator[1]
  preypop<-rbind(preypop,prey41)
  
  subPred<-subset(datfinprey,idpredator == as.character(idpred[42]))  
  prey42<-prey40
  prey42$idpredator<-subPred$idpredator[1]
  preypop<-rbind(preypop,prey42)
  
  subPred<-subset(datfinprey,idpredator == as.character(idpred[43]))  
  prey43<-prey40
  prey43$idpredator<-subPred$idpredator[1]
  preypop<-rbind(preypop,prey43)
  
  subPred<-subset(datfinprey,idpredator == as.character(idpred[44]))  
  prey44<-prey40
  prey44$idpredator<-subPred$idpredator[1]
  preypop<-rbind(preypop,prey44)

# now combine the datasets
  preypop<-preypop[,c(2,4,5,6,7,8,9,10,11)]
  colnames(preypop)<-c("stockid_prey","year","TB_prey","SSB_prey","F_prey","ER_prey","Fmsy_prey","idpredator","nb_prey")
  datfin<-merge(datfin, preypop, by.x=c("stockid", "year"), by.y=c("idpredator", "year"), all.x=TRUE, all.y=FALSE)

  dattest<-datfin[complete.cases(datfin[,c("ER","ER_prey")]),] # select years with ER data
  dattest$bio<-ifelse(is.na(dattest$TB ==T),dattest$SSB,dattest$TB) # use SSB if there is no total biomass 
  dattest<-dattest[complete.cases(dattest[,c("bio")]),] # verify there is biomass data for all years

# order per stock and year
  finaldata <- dattest[ order(dattest$stockid, dattest$year), ]
  
# clean 
  rm(list=setdiff(ls(), "finaldata"))
