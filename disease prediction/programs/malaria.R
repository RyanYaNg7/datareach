#Load in required libraries, if you have not installed these libraries, you can use the install.packages command 
library(data.table)
library(tidyverse)
library(forecast)
library(lmtest)
library(plyr)
library(lmtest)
library(ggplot2)
library(scales)
library(readr)

#read in the cleaned data file
data14_17 <- read_csv("/Users/Ryan/Desktop/my files/research/datareach/disease prediction/data file/Paludisme/tidyPaludisme14-17.csv")
data14_15 <- read_csv("/Users/Ryan/Desktop/my files/research/datareach/disease prediction/data file/Paludisme/tidyPaludisme14-15.csv")
data16_17 <- read_csv("/Users/Ryan/Desktop/my files/research/datareach/disease prediction/data file/Paludisme/tidyPaludisme16-17.csv")
data14 <- read_csv("/Users/Ryan/Desktop/my files/research/datareach/disease prediction/data file/Paludisme/tidyPaludisme14.csv")
data15 <- read_csv("/Users/Ryan/Desktop/my files/research/datareach/disease prediction/data file/Paludisme/tidyPaludisme15.csv")
data16 <- read_csv("/Users/Ryan/Desktop/my files/research/datareach/disease prediction/data file/Paludisme/tidyPaludisme16.csv")
data17 <- read_csv("/Users/Ryan/Desktop/my files/research/datareach/disease prediction/data file/Paludisme/tidyPaludisme17.csv")


#let's use this one as an example
data <- data14_17

# Further cleaning for data
# remove NA at the end
data <- data[,-c(ncol(data))]
# set 56th and 71st row to zero due to NA
data[56,] <- seq(0,0,length.out=ncol(data))
data[71,] <- seq(0,0,length.out=ncol(data))

# Convert the numerical data to numerical type
cols = c(4:ncol(data));    
data[,cols] = apply(data[,cols], 2, function(x) as.numeric(as.character(x)))

# replace NA with mean of predecessor and follower
N <- nrow(data)
data_noNA <- as_tibble(cbind(data))
for ( r in 1:N ){
  for ( c in cols ){
    if (is.na(data[r,c])){
      # find the prev
      has_prev <- FALSE
      for (i in c:4){
        if (!is.na(data[r,i])){
          has_prev <- TRUE
          prev <- data[r,i]
          break
        }
      } 
      # find post
      has_post <- FALSE
      for (j in c:ncol(data)){
        if (!is.na(data[r,j])){
          has_post <- TRUE
          post <- data[r,j]
          break
        }
      }
      data_noNA[r,c] <- 0
      if (has_prev){
        data_noNA[r,c] <- data_noNA[r,c]+prev
      }
      if (has_post){
        data_noNA[r,c] <- (data_noNA[r,c]+post)/(1+has_prev)
      }
    } 
  }
}

#show the result after removing NAs
r = 12
data_noNA[r,]
data[r,]
any(is.na(data_noNA))

#a function that take a row index, and return the a string
#in the format of "Region,Province"
index2region <- function(index){
  return( paste(data_noNA[index,3],data_noNA[index,2],sep=",") )
}

#a function that 
five_week_weighted_sum <- function(CrossCorr){
  score <- 0
  for (lag in -5:5){
    corr <- CrossCorr[lag][[1]]
    score <- score + 10.0/(abs(lag)+1)*corr
  }
  return(score)
}


### source: http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xegbohtmlnode39.html
### source: https://stackoverflow.com/questions/15551209/cross-correlation-using-ccf-in-r
calc_upperCI <- function(T, confidence){
  CI <- qnorm((1 + confidence)/2)/sqrt(T)
  return (CI)
}


# compute score matrix
score_mat <- matrix(nrow=N, ncol=N)
# the (i,j) entry of score mat correspond to ccf(i,j)
for (i in 1:N){
  for (j in 1:N){
    if (j > i){
      TS1<-as.ts(as.numeric(data_noNA[i, 4:ncol(data_noNA)])) 
      TS2<-as.ts(as.numeric(data_noNA[j, 4:ncol(data_noNA)]))
      CrossCorr<-ccf(TS1, TS2, lag = 5, pl = FALSE)
      score_mat[i,j] <- five_week_weighted_sum(CrossCorr)
    }
  }
}


process_scores<- function(score_mat){
  # summary stats of scores
  scores = c(score_mat)
  print(summary(scores))
  print ( c("Mean: ",mean(scores,na.rm = TRUE)) )
  print( c("SD: ",sd(scores,na.rm = TRUE)) )
  hist(scores,
       breaks=50, freq=FALSE,
       main="Histogram of Correlation Scores",
       xlab="Correlation Score",
       col="steelblue")
  # heatmap of score matrix
  melted_score = melt(score_mat, na.rm = TRUE) 
  colnames(melted_score) <- c("Region1", "Region2", "Score")
  ggplot(data=melted_score, aes(x=Region1, y=Region2, fill=Score)) +
    geom_tile(colour="white") +
    scale_fill_gradient2(low = "#a6611a", high = "#018571", mid = "#f5f5f5", midpoint = 0)
}



# Given two indices, plot their corrsponding time series
# if ccf_disp is TRUE, also plot the output graph of ccf function
visualize_corr <- function(index1, index2, ccf_disp=FALSE){
  TS1<-ts(as.numeric(data_noNA[index1, 4:ncol(data_noNA)])) 
  TS2<-ts(as.numeric(data_noNA[index2, 4:ncol(data_noNA)]))
  
  ts.plot(TS1, TS2,
          gpars=list(xlab="Week", ylab="Cases", col=c("red","blue")))
  legend("topleft",legend=c( paste(toString(index1),index2region(index1)), paste(toString(index2),index2region(index2))),
         col=c("red", "blue"),lty=c(1,1))
  
  CrossCorr<-ccf(TS1, TS2, lag = 5, pl = ccf_disp)
  score <- five_week_weighted_sum(CrossCorr)[[1]]
  
  title(sub = paste( "Correlation Score =",toString(score) ))
  if (ccf_disp){
    upperCI <- calc_upperCI(CrossCorr$n.used, 0.95)
    lowerCI <- -upperCI
    # verify CI
    lags <- -5:5
    my_upperCI <- rep(upperCI,11)
    points(lags, my_upperCI, col = "red")
    my_lowerCI <- rep(lowerCI,11)
    points(lags, my_lowerCI, col = "red")
  }
}


# Given a score matrix, and an interger k
# visualize_corr for top/bottom k scores(depending on the value of top)
# in the increasing order of scores
# ccf_disp dictates if the output graph of ccf is displayed
visualize_top_k <- function(score_mat, k, top=TRUE, ts_disp=TRUE, ccf_disp=FALSE){
  sorted_scores = sort(c(score_mat))
  bottom_k = sorted_scores[1:k]
  top_k = tail(sorted_scores,k)
  print(top_k)
  pairs<-c()
  if (top) {k_scores<-top_k} else {k_scores<-bottom_k}
  for (score in k_scores){
    index_pair = which(score_mat == score, arr.ind = TRUE)
    i<-index_pair[1]
    j<-index_pair[2]
    print(c(i, j))
    pairs<-c(pairs,paste(index2region(i),'X',index2region(j)))
    if (ts_disp) {
      visualize_corr(i,j,ccf_disp)
    }
  }
  print(pairs)
}


# 14-15 top10 region pairs
visualize_top_k(score_mat,10,top=TRUE,ts_disp=FALSE,ccf_disp=FALSE)


# 16-17 top10 region pairs
visualize_top_k(score_mat,10,top=TRUE,ts_disp=FALSE,ccf_disp=FALSE)


# 14-17 top10 region pairs
visualize_top_k(score_mat,10,top=TRUE,ts_disp=FALSE,ccf_disp=FALSE)



# outbreak detection thru approximated derivative
# this version use the whole time series and compute its mean and SD
# to detect any case that's out of CI
mark_outbreaks <- function(i,disp=FALSE){
  TS0<-as.ts(as.numeric(data_noNA[i, 4:ncol(data_noNA)])) #get the row in question
  LEN0<-length(TS0) # compute its length
  SD0<-sd(TS0) # compute its SD
  outbreaks = vector(,LEN0) 
  outbreaks[c(1,LEN0)]<-rep(0,2)
  outbreak_wks = c()
  for (t in 2:(LEN0-1)) {
    delta<-(TS0[t+1]-TS0[t-1])/SD0 # derivative approximation
    if (delta <= 2) {delta<-0} # threshold for outbreak: derivative > 2*SD
    else {outbreak_wks<-c(outbreak_wks,t)}
    outbreaks[t]<-delta
  }
  if (disp){
    ts.plot(TS0, gpars=list(xlab="Week", ylab="Cases", col=c("red")))
    ts.plot(outbreaks, gpars=list(xlab="Week", ylab="Delta", col=c("blue")))
  }
  return(outbreak_wks)
}


# Validation
# a region pair (i,j) with high correlation score in 2014-2015
# mark outbreaks in time series of i and j in 2016-2017
# to see if the outbreaks are concurrent
mark_outbreaks(54,disp=TRUE)
mark_outbreaks(81,disp=TRUE)


# try surveillance package on this data
library(surveillance)
disease_sts <- sts(observed = data_noNA[1, 4:ncol(data_noNA)],
                   start = c(2014, 01),
                   frequency = 52)

plot(disease_sts)
