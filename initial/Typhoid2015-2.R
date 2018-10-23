#TyphoidWHO dataframe storing Typhoid sheet (2015)
require(qdap)
library(qdap)
library(data.table)
library(forecast)
library(lmtest)
library(plyr)
library(lmtest)
require(surveillance)

TyphoidWHO<-TyphoidWHO2015

#first cut out the first 4 rows 
TyphoidWHOclean<-TyphoidWHO[-c(1:4),]
#Filter out columns by the region name (Calculated interval from looking at sheet)
TyphoidWHOclean<-TyphoidWHOclean[,c(1,2,seq(4,260,4))]

#Remove extra colums at the end and the cumulative columns
TyphoidWHOclean<-TyphoidWHOclean[,-c(56:67)]

#Remove the total columns to get only the data by each individual district 
TyphoidWHOclean<- TyphoidWHOclean[-c(10,41,56,87,112,128,148,169,180,199,200),]

#Eliminate the rows at the very end with the summary stats 
TyphoidWHOClean<-TyphoidWHOclean[-c(211:219),]
TyphoidWHOClean<-TyphoidWHOClean[-c(190,210),]
TyphoidWHOClean<-TyphoidWHOClean[c(1:189),] 

#TyphoidWHOClean is the set that now contains our cleaned data, not the one with lower case c

#TyphoidWHOClean[is.na(TyphoidWHOClean <- TyphoidWHOClean)] <- 0

#Cast columns to int type from original factor type 
#TyphoidWHOClean$X.2<-as.numeric(as.character(TyphoidWHOClean$X.2))
for(i in c(3:ncol(TyphoidWHOClean))){
  TyphoidWHOClean[,i] <- as.numeric(as.character(TyphoidWHOClean[,i]))
}

TyphoidWHOClean[is.na(TyphoidWHOClean <- TyphoidWHOClean)] <- 0

#Create time series object for the local region 
for(u in 89:92)
{
for(k in 89:92)
{
for(i in 1:7)
{
Adamaoua_Banyo<-as.ts(as.numeric(WHOTotConcat[k,3:ncol(WHOTotConcat)])) 
Adamaoua_Tigneri<-as.ts(as.numeric(WHOTotConcat[u,3:ncol(WHOTotConcat)])) 
plot.ts(Adamaoua_Banyo)
plot.ts(Adamaoua_Tigneri)
#TS_Typhoid_B=auto.arima(Adamaoua_Banyo)
#TS_Typh=auto.arima(Adamaoua_Tigneri)
#lines(TS_Typhoid_B$fitted, col="red")
#lines(TS_Typh$fitted, col="blue")
if(k!=u)
{
print(grangertest(Adamaoua_Tigneri~Adamaoua_Banyo, order=i)) 
print(grangertest(Adamaoua_Banyo~Adamaoua_Tigneri))
print(WHOTotConcat[k,2])
print(WHOTotConcat[u,2])
}
}

}
} 
#Convert to type numeric for columns 

cols = c(3:146);    
WHOTotConcat[,cols] = apply(WHOTotConcat[,cols], 2, function(x) as.numeric(as.character(x)))
                  
NBisson<-as.ts(as.numeric(WHOTotConcat[15, 3:ncol(WHOTotConcat)])) 
NDongo<-as.ts(as.numeric(WHOTotConcat[34, 3:ncol(WHOTotConcat)]))
CrossCorr<-ccf(NBisson, NDongo) 

sts1<-sts(NBisson)
surv<-earsC(sts1) 
plot(sts1)
plot(surv) 
