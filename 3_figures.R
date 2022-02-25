##This is the code for the figures for the analysis of foraging children's participation in play and work. 

##Note that data for Tsimane should be requested directly from Jonathan Stieglitz (jonathan.stieglitz@gmail.com)

##Due to file size limits, model outputs are [available for download as a zip file here](https://www.dropbox.com/s/chs9eyejzoernv2/post_rv.zip?dl=0) or upon request from Sheina Lew-Levy (sheinalewlevy@gmail.com)

##IF you don't have the rethinking package, load it following the instructions here https://xcelab.net/rm/software/

library(rethinking)
#######################################
###load posterior samples and data ####
#######################################
load(file="TAwaic.rda")
load(file="post1.rda")
load(file="post2.rda")
load(file="post3.rda")
load(file="post4.rda")
load(file="post5.rda")
load(file="post3.4.rda")
load(file="post4.3.rda")
load(file="post3.5.rda")
load(file="post3.6.rda")

d<-read.csv("dataset.csv")
ds<-subset(d,Society!="Dukha")

d$NPP_z<-(d$NPP-mean(d$NPP))/sd(d$NPP)
d$nonforaged_z<-(d$nonforaged-mean(d$nonforaged))/sd(d$nonforaged)
d$temp_z<-(d$meanAnnualTemp-mean(d$meanAnnualTemp))/sd(d$meanAnnualTemp)
d$prec_z<-(d$totalAnnualPrec-mean(d$totalAnnualPrec))/sd(d$totalAnnualPrec)
d$CV_z<-(d$precSeasonality-mean(d$precSeasonality))/sd(d$precSeasonality)
d$dens<-ifelse(d$sum_density>1,1,0)
d$snake_z<-(d$snake_count-mean(d$snake_count))/sd(d$snake_count)

ds$NPP_z<-(ds$NPP-mean(ds$NPP))/sd(ds$NPP)
ds$nonforaged_z<-(ds$nonforaged-mean(ds$nonforaged))/sd(ds$nonforaged)
ds$temp_z<-(ds$meanAnnualTemp-mean(ds$meanAnnualTemp))/sd(ds$meanAnnualTemp)
ds$prec_z<-(ds$totalAnnualPrec-mean(ds$totalAnnualPrec))/sd(ds$totalAnnualPrec)

dc<-subset(d, Ado!=1)
dc$NPP_z<-(dc$NPP-mean(dc$NPP))/sd(dc$NPP)
dc$nonforaged_z<-(dc$nonforaged-mean(dc$nonforaged))/sd(dc$nonforaged)
dc$temp_z<-(dc$meanAnnualTemp-mean(dc$meanAnnualTemp))/sd(dc$meanAnnualTemp)
dc$prec_z<-(dc$totalAnnualPrec-mean(dc$totalAnnualPrec))/sd(dc$totalAnnualPrec)

##################
###WAIC results###
##################
waicTA
plot(waicTA)

###############################################
##Correlation matrix from model 1##############
###############################################
quartz(height = 6 , width = 10)
v_est <- apply(post1$v_id,2:3,median)
my_line <- function(x,y,...){
  points(x,y,...)
  abline(a = lm(y ~ x)$coefficients[1] , b = lm(y ~ x)$coefficients[2] , ...)
}
labels<-c("Childcare","Food production","Domestic work","Play")
pairs(v_est,upper.panel=my_line,lower.panel=NULL,col="black",pch = 21,labels=labels)

quartz.save(file="FigS1",type="pdf")
dev.off()
rm(v_est)
rm(my_line)
rm(labels)
######################################
###Generate predictions for Model 2###
######################################

link.TA2 <- function( data ) {
  K <- dim(post2$v_id)[3] + 1
  ns <- dim(post2$v_id)[1]
  if ( missing(data) ) stop( "BOOM: Need data argument" )
  n <- seq.length
  
  softmax2 <- function(x) {
    x <- max(x) - x
    exp(-x)/sum(exp(-x))
  }
  
  p <- list()
  
  for ( i in 1:n ) {
    p[[i]] <- sapply( 1:K , function(k) {
      if ( k < K ) {
        ptemp <- post2$a[,k] + 
          post2$b_sex[,k] * data$sex[i] + 
          post2$b_middle[,k] * data$middle[i] + 
          post2$b_ado[,k] * data$ado[i] + 
          post2$b_sex_middle[,k] * data$sex[i] * data$middle[i] + 
          post2$b_sex_ado[,k] * data$sex[i] * data$ado[i] + 
          post2$b_nonforaged[,k] * data$nonforaged[i]
        if ( data$id[i]>0 ) ptemp <- ptemp + post2$v_id[,data$id[i],k]
        if ( data$society[i]>0 ) ptemp <- ptemp + post2$v_society[,data$society[i],k]
      } else {
        ptemp <- rep(0,ns)
      }
      return(ptemp)
    })
    ## The values are converted to probabilities using the softmax function
    ## which ensures that the predicted values across categories sum to
    ## 100% probabilities.
    for ( s in 1:ns ) p[[i]][s,] <- softmax2( p[[i]][s,] )
  }
  return(p)
}

#######################################
###Make correlation plot for Model 2###
#######################################
quartz(height = 6 , width = 10)
v_est <- apply(post2$v_id,2:3,median)
my_line <- function(x,y,...){
  points(x,y,...)
  abline(a = lm(y ~ x)$coefficients[1] , b = lm(y ~ x)$coefficients[2] , ...)
}
labels<-c("Childcare","Food production","Domestic work","Play")
pairs(v_est,upper.panel=my_line,lower.panel=NULL,col="black",pch = 21,labels=labels)

quartz.save(file="FigS2",type="pdf")
dev.off()
rm(v_est)
rm(my_line)
rm(labels)
rm(link.TA2)

######################################
###Generate predictions for Model 5###
######################################

link.TA5 <- function( data ) {
  K <- dim(post5$v_id)[3] + 1
  ns <- dim(post5$v_id)[1]
  if ( missing(data) ) stop( "BOOM: Need data argument" )
  n <- seq.length
  
  softmax2 <- function(x) {
    x <- max(x) - x
    exp(-x)/sum(exp(-x))
  }
  
  p <- list()
  
  for ( i in 1:n ) {
    p[[i]] <- sapply( 1:K , function(k) {
      if ( k < K ) {
        ptemp <- post5$a[,k] + 
          post5$b_sex[,k] * data$sex[i] + 
          post5$b_middle[,k] * data$middle[i] + 
          post5$b_ado[,k] * data$ado[i] + 
          post5$b_sex_middle[,k] * data$sex[i] * data$middle[i] + 
          post5$b_sex_ado[,k] * data$sex[i] * data$ado[i] + 
          post5$b_nonforaged[,k] * data$nonforaged[i]+
          post5$b_div[,k] * data$div[i]+
          post5$b_sex_div[,k] * data$div[i] * data$sex[i]
        if ( data$id[i]>0 ) ptemp <- ptemp + post5$v_id[,data$id[i],k]
        if ( data$society[i]>0 ) ptemp <- ptemp + post5$v_society[,data$society[i],k]
      } else {
        ptemp <- rep(0,ns)
      }
      return(ptemp)
    })
    ## The values are converted to probabilities using the softmax function
    ## which ensures that the predicted values across categories sum to
    ## 100% probabilities.
    for ( s in 1:ns ) p[[i]][s,] <- softmax2( p[[i]][s,] )
  }
  return(p)
}


#####################################################################
###Make in-text Gendered Division of Food Production Labour Figure###
#####################################################################
pred_data_lowdiv<- data.frame(
  id = 0 ,
  society=0,
  sex = c(0,1),
  middle = 0,
  ado = 1,
  nonforaged=mean(d$nonforaged_z),
  div=-0.5
)

pred_data_highdiv<- data.frame(
  id = 0 ,
  society=0,
  sex = c(0,1),
  middle = 0,
  ado = 1,
  nonforaged=mean(d$nonforaged_z),
  div=1
)

seq.length<-2
p_lowdiv <- link.TA5 (pred_data_lowdiv)
p_lowdiv_mean <- sapply( 1:length(p_lowdiv) , function(i) apply(p_lowdiv[[i]],2,mean) )

p_highdiv <- link.TA5 (pred_data_highdiv)
p_highdiv_mean <- sapply( 1:length(p_highdiv) , function(i) apply(p_highdiv[[i]],2,mean) )


quartz(height = 6 , width = 18)
par(mfrow=c(1,2), mar=c(1,2,1,2) , oma=c(5,2.5,2.5,2.5))

plot( NULL , xlim=c(-0.5,3.5) , ylim=c(0,0.5) , xaxt = "n", main = " ",adj=0, ylab = "", xlab = "", cex.main = .9, cex.axis=1.2)
for ( k in 1:1 ) {
  p_lowdiv_PI <- sapply( 1:length(p_lowdiv) , function(i) PI(p_lowdiv[[i]][,k],prob=0.89) ) 
  segments(-0.1,p_lowdiv_PI[1,1],-0.1,p_lowdiv_PI[2,1],col=col.alpha("dark orange",1),lwd=3)
  segments(0.1,p_lowdiv_PI[1,2],0.1,p_lowdiv_PI[2,2],col=col.alpha("purple",1),lwd=3)
  points( c(-0.1,0.1) , p_lowdiv_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 2:2) {
  p_lowdiv_PI <- sapply( 1:length(p_lowdiv) , function(i) PI(p_lowdiv[[i]][,k],prob=0.89) ) 
  segments(0.9,p_lowdiv_PI[1,1],0.9,p_lowdiv_PI[2,1],col=col.alpha("dark orange",1),lwd=3)
  segments(1.1,p_lowdiv_PI[1,2],1.1,p_lowdiv_PI[2,2],col=col.alpha("purple",1),lwd=3)
  points( c(0.9,1.1) , p_lowdiv_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 3:3) {
  p_lowdiv_PI <- sapply( 1:length(p_lowdiv) , function(i) PI(p_lowdiv[[i]][,k],prob=0.89) ) 
  segments(1.9,p_lowdiv_PI[1,1],1.9,p_lowdiv_PI[2,1],col=col.alpha("dark orange",1),lwd=3)
  segments(2.1,p_lowdiv_PI[1,2],2.1,p_lowdiv_PI[2,2],col=col.alpha("purple",1),lwd=3)
  points( c(1.9,2.1) , p_lowdiv_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 4:4) {
  p_lowdiv_PI <- sapply( 1:length(p_lowdiv) , function(i) PI(p_lowdiv[[i]][,k],prob=0.89) ) 
  segments(2.9,p_lowdiv_PI[1,1],2.9,p_lowdiv_PI[2,1],col=col.alpha("dark orange",1),lwd=3)
  segments(3.1,p_lowdiv_PI[1,2],3.1,p_lowdiv_PI[2,2],col=col.alpha("purple",1),lwd=3)
  points( c(2.9,3.1) , p_lowdiv_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

##for ( k in 5:5) {
 ## p_lowdiv_PI <- sapply( 1:length(p_lowdiv) , function(i) PI(p_lowdiv[[i]][,k],prob=0.89) ) 
 ## segments(3.9,p_lowdiv_PI[1,1],3.9,p_lowdiv_PI[2,1],col=col.alpha("dark orange",1),lwd=3)
 ## segments(4.1,p_lowdiv_PI[1,2],4.1,p_lowdiv_PI[2,2],col=col.alpha("purple",1),lwd=3)
 ## points( c(3.9,4.1) , p_lowdiv_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
##}
labels.at<-c(0,1,2,3)
axis( side = 1, at = labels.at, labels = c("Childcare","Food Production","Domestic Work","Play"), cex.axis=1.3)

mtext(text="(A) Female-Biased Gendered Division of Food Production Labour",side=1,line=3,cex=1.4)
mtext(text="Probability (activity)",side=2,line=2, cex=1.4)


plot( NULL , xlim=c(-0.5,3.5) , ylim=c(0,0.5) , xaxt = "n",  main = " ",adj=0, ylab = "", xlab = "", cex.main = .9, cex.axis=1, col.axis = "white")
for ( k in 1:1 ) {
  p_highdiv_PI <- sapply( 1:length(p_highdiv) , function(i) PI(p_highdiv[[i]][,k],prob=0.89) ) 
  segments(-0.1,p_highdiv_PI[1,1],-0.1,p_highdiv_PI[2,1],col=col.alpha("dark orange",1),lwd=3)
  segments(0.1,p_highdiv_PI[1,2],0.1,p_highdiv_PI[2,2],col=col.alpha("purple",1),lwd=3)
  points( c(-0.1,0.1) , p_highdiv_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 2:2) {
  p_highdiv_PI <- sapply( 1:length(p_highdiv) , function(i) PI(p_highdiv[[i]][,k],prob=0.89) ) 
  segments(0.9,p_highdiv_PI[1,1],0.9,p_highdiv_PI[2,1],col=col.alpha("dark orange",1),lwd=3)
  segments(1.1,p_highdiv_PI[1,2],1.1,p_highdiv_PI[2,2],col=col.alpha("purple",1),lwd=3)
  points( c(0.9,1.1) , p_highdiv_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 3:3) {
  p_highdiv_PI <- sapply( 1:length(p_highdiv) , function(i) PI(p_highdiv[[i]][,k],prob=0.89) ) 
  segments(1.9,p_highdiv_PI[1,1],1.9,p_highdiv_PI[2,1],col=col.alpha("dark orange",1),lwd=3)
  segments(2.1,p_highdiv_PI[1,2],2.1,p_highdiv_PI[2,2],col=col.alpha("purple",1),lwd=3)
  points( c(1.9,2.1) , p_highdiv_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 4:4) {
  p_highdiv_PI <- sapply( 1:length(p_highdiv) , function(i) PI(p_highdiv[[i]][,k],prob=0.89) ) 
  segments(2.9,p_highdiv_PI[1,1],2.9,p_highdiv_PI[2,1],col=col.alpha("dark orange",1),lwd=3)
  segments(3.1,p_highdiv_PI[1,2],3.1,p_highdiv_PI[2,2],col=col.alpha("purple",1),lwd=3)
  points( c(2.9,3.1) , p_highdiv_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

##for ( k in 5:5) {
##  p_highdiv_PI <- sapply( 1:length(p_highdiv) , function(i) PI(p_highdiv[[i]][,k],prob=0.89) ) 
##  segments(3.9,p_highdiv_PI[1,1],3.9,p_highdiv_PI[2,1],col=col.alpha("dark orange",1),lwd=3)
##  segments(4.1,p_highdiv_PI[1,2],4.1,p_highdiv_PI[2,2],col=col.alpha("purple",1),lwd=3)
##  points( c(3.9,4.1) , p_highdiv_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
##}
##labels.at<-c(0,1,2,3,4)
axis( side = 1, at = labels.at, labels = c("Childcare","Food Production","Domestic Work","Play"), cex.axis=1.3)
mtext(text="(B) Male-Biased Gendered Division of Food Production Labour",side=1,line=3,cex=1.4)

legend("topright", 
       inset=0, 
       cex = 1.2, 
       c("Girls","Boys"), 
       horiz=FALSE, 
       lty=c(1,1), 
       lwd=c(2,2), 
       col=c("darkorange","purple"), 
       title="",
       bty="n",
       text.font=1.3)



quartz.save(file="fig4.pdf",type="pdf")
dev.off()
rm(p_highdiv)
rm(p_highdiv_mean)
rm(p_highdiv_PI)
rm(p_lowdiv)
rm(p_lowdiv_mean)
rm(p_lowdiv_PI)
rm(pred_data_highdiv)
rm(pred_data_lowdiv)
rm(seq.length)
rm(labels.at)

###########################################################################
###Make supplementary Gendered Division of Food Production Labour Figure###
###########################################################################

for(t in unique(d$variable)) {d[paste("r_",t,sep="")] <- ifelse(d$variable==t,1,0)}
agg <- aggregate ( cbind (r_food_production, r_znonworkobs, r_play, r_childcare, r_household) ~ Society+ div+Sex, data = d, FUN = mean)
agg_f<-subset(agg,Sex==0)
agg_m<-subset(agg,Sex==1)

seq.length <- 100
div_seq <- seq(from= min(d$div) , to= max(d$div) , length.out = seq.length)

divf <- data.frame(
  id = 0 ,
  society=0,
  sex = 0,
  middle = 1,
  ado = 0,
  nonforaged=mean(d$nonforaged_z),
  div=div_seq
)

p_divf <- link.TA5 (divf)
p_divf_mean <- sapply( 1:length(p_divf) , function(i) apply(p_divf[[i]],2,mean) )

divm <- data.frame(
  id = 0 ,
  society=0,
  sex = 1,
  middle = 1,
  ado = 0,
  nonforaged=mean(d$nonforaged_z),
  div=div_seq
)

p_divm <- link.TA5 (divm)
p_divm_mean <- sapply( 1:length(p_divm) , function(i) apply(p_divm[[i]],2,mean) )

labels.at <- c(-0.5,0,0.5,1)

quartz(height = 6 , width = 10)
par(mfrow=c(2,3), mar=c(1,2,1,2) + 0.1 , oma=c(2.5,4,2.5,2.5))

plot( NULL , xlim=c(min(d$div),max(d$div)) , ylim=c(0,0.35) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 1:1 ) {
  lines( div_seq , p_divf_mean[k,], lwd = 2 ,col=col.alpha("dark orange",1))
  p_divf_PI <- sapply( 1:length(p_divf) , function(i) PI(p_divf[[i]][,k],prob=0.89) )
  shade( p_divf_PI , div_seq,col=col.alpha("dark orange",0.1))
}
points(agg_f$div, agg_f$r_childcare,col="dark orange")
for ( k in 1:1 ) {
  lines( div_seq , p_divm_mean[k,], lwd = 2 ,col=col.alpha("purple",1))
  p_divm_PI <- sapply( 1:length(p_divm) , function(i) PI(p_divm[[i]][,k],prob=0.89) )
  shade( p_divm_PI , div_seq,col=col.alpha("purple",0.1))
}
points(agg_m$div, agg_m$r_childcare,col="purple")
axis( side = 1, at = labels.at, labels = F)
mtext(text="Probability (Childcare)",side=2,line=2)

legend("topright", 
       inset=0, 
       cex = 0.75, 
       c("Girls","Boys"), 
       horiz=TRUE, 
       lty=c(1,1), 
       lwd=c(2,2), 
       col=c("darkorange","purple"), 
       bg="grey96",
       title="",
       text.font=0.5)

plot( NULL , xlim=c(min(d$div),max(d$div)) , ylim=c(0,0.35) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 2:2 ) {
  lines( div_seq , p_divf_mean[k,], lwd = 2 ,col=col.alpha("dark orange",1))
  p_divf_PI <- sapply( 1:length(p_divf) , function(i) PI(p_divf[[i]][,k],prob=0.89) )
  shade( p_divf_PI , div_seq,col=col.alpha("dark orange",0.1))
}
points(agg_f$div, agg_f$r_food_production,col="dark orange")
for ( k in 2:2 ) {
  lines( div_seq , p_divm_mean[k,], lwd = 2 ,col=col.alpha("purple",1))
  p_divm_PI <- sapply( 1:length(p_divm) , function(i) PI(p_divm[[i]][,k],prob=0.89) )
  shade( p_divm_PI , div_seq,col=col.alpha("purple",0.1))
}
points(agg_m$div, agg_m$r_food_production,col="purple")
axis( side = 1, at = labels.at, labels = F)
mtext(text="Probability (Food Production)",side=2,line=2)

plot( NULL , xlim=c(min(d$div),max(d$div)) , ylim=c(0,0.35) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 3:3 ) {
  lines( div_seq , p_divf_mean[k,], lwd = 2 ,col=col.alpha("dark orange",1))
  p_divf_PI <- sapply( 1:length(p_divf) , function(i) PI(p_divf[[i]][,k],prob=0.89) )
  shade( p_divf_PI , div_seq,col=col.alpha("dark orange",0.1))
}
points(agg_f$div, agg_f$r_household,col="dark orange")
for ( k in 3:3 ) {
  lines( div_seq , p_divm_mean[k,], lwd = 2 ,col=col.alpha("purple",1))
  p_divm_PI <- sapply( 1:length(p_divm) , function(i) PI(p_divm[[i]][,k],prob=0.89) )
  shade( p_divm_PI , div_seq,col=col.alpha("purple",0.1))
}
points(agg_m$div, agg_m$r_household,col="purple")
axis( side = 1, at = labels.at, labels = labels.at)
mtext(text="Probability (Domestic Work)",side=2,line=2)

plot( NULL , xlim=c(min(d$div),max(d$div)) , ylim=c(0,0.35) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 4:4 ) {
  lines( div_seq , p_divf_mean[k,], lwd = 2 ,col=col.alpha("dark orange",1))
  p_divf_PI <- sapply( 1:length(p_divf) , function(i) PI(p_divf[[i]][,k],prob=0.89) )
  shade( p_divf_PI , div_seq,col=col.alpha("dark orange",0.1))
}
points(agg_f$div, agg_f$r_play,col="dark orange")
for ( k in 4:4 ) {
  lines( div_seq , p_divm_mean[k,], lwd = 2 ,col=col.alpha("purple",1))
  p_divm_PI <- sapply( 1:length(p_divm) , function(i) PI(p_divm[[i]][,k],prob=0.89) )
  shade( p_divm_PI , div_seq,col=col.alpha("purple",0.1))
}
points(agg_m$div, agg_m$r_play,col="purple")
axis( side = 1, at = labels.at, labels = labels.at)
mtext(text="Probability (Play)",side=2,line=2)

plot( NULL , xlim=c(min(d$div),max(d$div)) , ylim=c(0,0.8) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 5:5 ) {
  lines( div_seq , p_divf_mean[k,], lwd = 2 ,col=col.alpha("dark orange",1))
  p_divf_PI <- sapply( 1:length(p_divf) , function(i) PI(p_divf[[i]][,k],prob=0.89) )
  shade( p_divf_PI , div_seq,col=col.alpha("dark orange",0.1))
}
points(agg_f$div, agg_f$r_znonworkobs,col="dark orange")
for ( k in 5:5 ) {
  lines( div_seq , p_divm_mean[k,], lwd = 2 ,col=col.alpha("purple",1))
  p_divm_PI <- sapply( 1:length(p_divm) , function(i) PI(p_divm[[i]][,k],prob=0.89) )
  shade( p_divm_PI , div_seq,col=col.alpha("purple",0.1))
}
points(agg_m$div, agg_m$r_znonworkobs,col="purple")
axis( side = 1, at = labels.at, labels = labels.at)
mtext(text="Probability (Other Activities)",side=2,line=2)

mtext(text="Gendered Division of Food Production Labour",side=1,line=1,outer=TRUE)


quartz.save(file="figS7.pdf",type="pdf")
dev.off()
rm(agg)
rm(agg_f)
rm(agg_m)
rm(divm)
rm(divf)
rm(p_divf)
rm(p_divf_mean)
rm(p_divf_PI)
rm(p_divm)
rm(p_divm_mean)
rm(p_divm_PI)
rm(div_seq)
rm(labels.at)
rm(seq.length)
rm(link.TA5)
rm(t)
rm(k)

######################################
###Generate predictions for Model 3###
######################################

link.TA3 <- function( data ) {
  K <- dim(post3$v_id)[3] + 1
  ns <- dim(post3$v_id)[1]
  if ( missing(data) ) stop( "BOOM: Need data argument" )
  n <- seq.length
  
  softmax2 <- function(x) {
    x <- max(x) - x
    exp(-x)/sum(exp(-x))
  }
  
  p <- list()
  
  for ( i in 1:n ) {
    p[[i]] <- sapply( 1:K , function(k) {
      if ( k < K ) {
        ptemp <- post3$a[,k] + 
          post3$b_sex[,k] * data$sex[i] + 
          post3$b_middle[,k] * data$middle[i] + 
          post3$b_ado[,k] * data$ado[i] + 
          post3$b_sex_middle[,k] * data$sex[i] * data$middle[i] + 
          post3$b_sex_ado[,k] * data$sex[i] * data$ado[i] + 
          post3$b_nonforaged[,k] * data$nonforaged[i]+
          post3$b_NPP[,k] * data$NPP[i]+
          post3$b_temp[,k] * data$temp[i]+
          post3$b_prec[,k] * data$prec[i]
        if ( data$id[i]>0 ) ptemp <- ptemp + post3$v_id[,data$id[i],k]
        if ( data$society[i]>0 ) ptemp <- ptemp + post3$v_society[,data$society[i],k]
      } else {
        ptemp <- rep(0,ns)
      }
      return(ptemp)
    })
    ## The values are converted to probabilities using the softmax function
    ## which ensures that the predicted values across categories sum to
    ## 100% probabilities.
    for ( s in 1:ns ) p[[i]][s,] <- softmax2( p[[i]][s,] )
  }
  return(p)
}

#####################
###Make NPP Figure###
#####################
seq.length <- 100
NPP_seq <- seq(from= min(d$NPP_z) , to= max(d$NPP_z), length.out = seq.length)
agg_NPP <- aggregate ( cbind (r_food_production, r_znonworkobs, r_play, r_childcare, r_household) ~ Society+ NPP_z, data = d, FUN = mean)

NPP <- data.frame(
  id = 0 ,
  society=0,
  sex = 0,
  middle = 0,
  ado = 0,
  nonforaged=mean(d$nonforaged_z),
  NPP=NPP_seq,
  temp=mean(d$temp_z),
  prec=mean(d$prec_z)
)

p_NPP <- link.TA3 (NPP)
p_NPP_mean <- sapply( 1:length(p_NPP) , function(i) apply(p_NPP[[i]],2,mean) )
preferred.NPP<-c(150,650,1150,1650,2150,2650)
labels.at <- (preferred.NPP - mean(d$NPP))/sd(d$NPP)

quartz(height = 6 , width = 10)
par(mfrow=c(2,3), mar=c(1,2,1,2) + 0.1 , oma=c(2.5,4,2.5,2.5))
plot( NULL , xlim=c(min(d$NPP_z),max(d$NPP_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 1:1 ) {
  lines( NPP_seq , p_NPP_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_NPP_PI <- sapply( 1:length(p_NPP) , function(i) PI(p_NPP[[i]][,k],prob=0.89) )
  shade( p_NPP_PI , NPP_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = F)
points(agg_NPP$NPP_z, agg_NPP$r_childcare)
mtext(text="Probability (Childcare)",side=2,line=2)

plot( NULL , xlim=c(min(d$NPP_z),max(d$NPP_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 2:2 ) {
  lines( NPP_seq , p_NPP_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_NPP_PI <- sapply( 1:length(p_NPP) , function(i) PI(p_NPP[[i]][,k],prob=0.89) )
  shade( p_NPP_PI , NPP_seq,col=col.alpha("black",0.1))
}
points(agg_NPP$NPP_z, agg_NPP$r_food_production)

axis( side = 1, at = labels.at, labels = F)
mtext(text="Probability (Food Production)",side=2,line=2)

plot( NULL , xlim=c(min(d$NPP_z),max(d$NPP_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 3:3 ) {
  lines( NPP_seq , p_NPP_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_NPP_PI <- sapply( 1:length(p_NPP) , function(i) PI(p_NPP[[i]][,k],prob=0.89) )
  shade( p_NPP_PI , NPP_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.NPP)
points(agg_NPP$NPP_z, agg_NPP$r_household)
mtext(text="Probability (Domestic Work)",side=2,line=2)

plot( NULL , xlim=c(min(d$NPP_z),max(d$NPP_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 4:4 ) {
  lines( NPP_seq , p_NPP_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_NPP_PI <- sapply( 1:length(p_NPP) , function(i) PI(p_NPP[[i]][,k],prob=0.89) )
  shade( p_NPP_PI , NPP_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.NPP)
points(agg_NPP$NPP_z, agg_NPP$r_play)
mtext(text="Probability (Play)",side=2,line=2)

plot( NULL , xlim=c(min(d$NPP_z),max(d$NPP_z)) , ylim=c(0,0.8) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 5:5 ) {
  lines( NPP_seq , p_NPP_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_NPP_PI <- sapply( 1:length(p_NPP) , function(i) PI(p_NPP[[i]][,k],prob=0.89) )
  shade( p_NPP_PI , NPP_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.NPP)
points(agg_NPP$NPP_z, agg_NPP$r_znonworkobs)
mtext(text="Probability (Other Activities)",side=2,line=2)

mtext(text="Net Primary Productivity (NPP)",side=1,line=1,outer=TRUE)
quartz.save(file="figS3.pdf",type="pdf")
rm(agg_NPP)
rm(NPP)
rm(p_NPP)
rm(p_NPP_mean)
rm(p_NPP_PI)
rm(k)
rm(labels.at)
rm(NPP_seq)
rm(preferred.NPP)
rm(seq.length)
dev.off()

######################
###Make temp Figure###
######################
seq.length <- 100
temp_seq <- seq(from= min(d$temp_z) , to= max(d$temp_z), length.out = seq.length)
agg_temp <- aggregate ( cbind (r_food_production, r_znonworkobs, r_play, r_childcare, r_household) ~ Society+ temp_z, data = d, FUN = mean)

temp <- data.frame(
  id = 0 ,
  society=0,
  sex = 0,
  middle = 0,
  ado = 0,
  nonforaged=mean(d$nonforaged_z),
  NPP=mean(d$NPP_z),
  temp=temp_seq,
  prec=mean(d$prec_z)
)

p_temp <- link.TA3 (temp)
p_temp_mean <- sapply( 1:length(p_temp) , function(i) apply(p_temp[[i]],2,mean) )
preferred.temp<-c(-5,5,15,25)
labels.at <- (preferred.temp - mean(d$meanAnnualTemp))/sd(d$meanAnnualTemp)

quartz(height = 6 , width = 10)
par(mfrow=c(2,3), mar=c(1,2,1,2) + 0.1 , oma=c(2.5,4,2.5,2.5))
plot( NULL , xlim=c(min(d$temp_z),max(d$temp_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 1:1 ) {
  lines( temp_seq , p_temp_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_PI <- sapply( 1:length(p_temp) , function(i) PI(p_temp[[i]][,k],prob=0.89) )
  shade( p_temp_PI , temp_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = F)
points(agg_temp$temp_z, agg_temp$r_childcare)
mtext(text="Probability (Childcare)",side=2,line=2)

plot( NULL , xlim=c(min(d$temp_z),max(d$temp_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 2:2 ) {
  lines( temp_seq , p_temp_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_PI <- sapply( 1:length(p_temp) , function(i) PI(p_temp[[i]][,k],prob=0.89) )
  shade( p_temp_PI , temp_seq,col=col.alpha("black",0.1))
}
points(agg_temp$temp_z, agg_temp$r_food_production)

axis( side = 1, at = labels.at, labels = F)
mtext(text="Probability (Food Production)",side=2,line=2)

plot( NULL , xlim=c(min(d$temp_z),max(d$temp_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 3:3 ) {
  lines( temp_seq , p_temp_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_PI <- sapply( 1:length(p_temp) , function(i) PI(p_temp[[i]][,k],prob=0.89) )
  shade( p_temp_PI , temp_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.temp)
points(agg_temp$temp_z, agg_temp$r_household)
mtext(text="Probability (Domestic Work)",side=2,line=2)

plot( NULL , xlim=c(min(d$temp_z),max(d$temp_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 4:4 ) {
  lines( temp_seq , p_temp_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_PI <- sapply( 1:length(p_temp) , function(i) PI(p_temp[[i]][,k],prob=0.89) )
  shade( p_temp_PI , temp_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.temp)
points(agg_temp$temp_z, agg_temp$r_play)
mtext(text="Probability (Play)",side=2,line=2)

plot( NULL , xlim=c(min(d$temp_z),max(d$temp_z)) , ylim=c(0,1) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 5:5) {
  lines( temp_seq , p_temp_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_PI <- sapply( 1:length(p_temp) , function(i) PI(p_temp[[i]][,k],prob=0.89) )
  shade( p_temp_PI , temp_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.temp)
points(agg_temp$temp_z, agg_temp$r_znonworkobs)
mtext(text="Probability (Other Activities)",side=2,line=2)

mtext(text="Annual Mean Temperature (°C)",side=1,line=1,outer=TRUE)

quartz.save(file="figS4.pdf",type="pdf")
rm(agg_temp)
rm(p_temp)
rm(p_temp_mean)
rm(p_temp_PI)
rm(temp)
rm(k)
rm(labels.at)
rm(preferred.temp)
rm(seq.length)
rm(temp_seq)
dev.off()

######################
###Make prec Figure###
######################
seq.length <- 100
prec_seq <- seq(from= min(d$prec_z) , to= max(d$prec_z), length.out = seq.length)
agg_prec <- aggregate ( cbind (r_food_production, r_znonworkobs, r_play, r_childcare, r_household) ~ Society+ prec_z, data = d, FUN = mean)

prec <- data.frame(
  id = 0 ,
  society=0,
  sex = 0,
  middle = 0,
  ado = 0,
  nonforaged=mean(d$nonforaged_z),
  NPP=mean(d$NPP_z),
  temp=mean(d$temp_z),
  prec=prec_seq
)

p_prec <- link.TA3 (prec)
p_prec_mean <- sapply( 1:length(p_prec) , function(i) apply(p_prec[[i]],2,mean) )

preferred.prec<-c(500,1000,1500,2000,2500)
labels.at <- (preferred.prec - mean(d$totalAnnualPrec))/sd(d$totalAnnualPrec)

quartz(height = 6 , width = 10)
par(mfrow=c(2,3), mar=c(1,2,1,2) + 0.1 , oma=c(2.5,4,2.5,2.5))
plot( NULL , xlim=c(min(d$prec_z),max(d$prec_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 1:1 ) {
  lines( prec_seq , p_prec_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_prec_PI <- sapply( 1:length(p_prec) , function(i) PI(p_prec[[i]][,k],prob=0.89) )
  shade( p_prec_PI , prec_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = F)
points(agg_prec$prec_z, agg_prec$r_childcare)
mtext(text="Probability (Childcare)",side=2,line=2)

plot( NULL , xlim=c(min(d$prec_z),max(d$prec_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 2:2 ) {
  lines( prec_seq , p_prec_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_prec_PI <- sapply( 1:length(p_prec) , function(i) PI(p_prec[[i]][,k],prob=0.89) )
  shade( p_prec_PI , prec_seq,col=col.alpha("black",0.1))
}
points(agg_prec$prec_z, agg_prec$r_food_production)

axis( side = 1, at = labels.at, labels = F)
mtext(text="Probability (Food Production)",side=2,line=2)

plot( NULL , xlim=c(min(d$prec_z),max(d$prec_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 3:3 ) {
  lines( prec_seq , p_prec_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_prec_PI <- sapply( 1:length(p_prec) , function(i) PI(p_prec[[i]][,k],prob=0.89) )
  shade( p_prec_PI , prec_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.prec)
points(agg_prec$prec_z, agg_prec$r_household)
mtext(text="Probability (Domestic Work)",side=2,line=2)

plot( NULL , xlim=c(min(d$prec_z),max(d$prec_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 4:4 ) {
  lines( prec_seq , p_prec_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_prec_PI <- sapply( 1:length(p_prec) , function(i) PI(p_prec[[i]][,k],prob=0.89) )
  shade( p_prec_PI , prec_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.prec)
points(agg_prec$prec_z, agg_prec$r_play)
mtext(text="Probability (Play)",side=2,line=2)

plot( NULL , xlim=c(min(d$prec_z),max(d$prec_z)) , ylim=c(0,0.8) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 5:5) {
  lines( prec_seq , p_prec_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_prec_PI <- sapply( 1:length(p_prec) , function(i) PI(p_prec[[i]][,k],prob=0.89) )
  shade( p_prec_PI , prec_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.prec)
points(agg_prec$prec_z, agg_prec$r_znonworkobs)
mtext(text="Probability (Other Activities)",side=2,line=2)
mtext(text="Annual Precipitation (mm)",side=1,line=1,outer=TRUE)

quartz.save(file="figS5.pdf",type="pdf")
rm(agg_prec)
rm(p_prec)
rm(p_prec_mean)
rm(p_prec_PI)
rm(prec)
rm(k)
rm(labels.at)
rm(prec_seq)
rm(preferred.prec)
rm(seq.length)
dev.off()
rm(link.TA3)

######################################
###Generate predictions for Model 4###
######################################

link.TA4 <- function( data ) {
  K <- dim(post4$v_id)[3] + 1
  ns <- dim(post4$v_id)[1]
  if ( missing(data) ) stop( "BOOM: Need data argument" )
  n <- seq.length
  
  softmax2 <- function(x) {
    x <- max(x) - x
    exp(-x)/sum(exp(-x))
  }
  
  p <- list()
  
  for ( i in 1:n ) {
    p[[i]] <- sapply( 1:K , function(k) {
      if ( k < K ) {
        ptemp <- post4$a[,k] + 
          post4$b_sex[,k] * data$sex[i] + 
          post4$b_middle[,k] * data$middle[i] + 
          post4$b_ado[,k] * data$ado[i] + 
          post4$b_sex_middle[,k] * data$sex[i] * data$middle[i] + 
          post4$b_sex_ado[,k] * data$sex[i] * data$ado[i] + 
          post4$b_nonforaged[,k] * data$nonforaged[i]+
          post4$b_water[,k] * data$water[i]+
          post4$b_dens[,k] * data$dens[i]
        if ( data$id[i]>0 ) ptemp <- ptemp + post4$v_id[,data$id[i],k]
        if ( data$society[i]>0 ) ptemp <- ptemp + post4$v_society[,data$society[i],k]
      } else {
        ptemp <- rep(0,ns)
      }
      return(ptemp)
    })
    ## The values are converted to probabilities using the softmax function
    ## which ensures that the predicted values across categories sum to
    ## 100% probabilities.
    for ( s in 1:ns ) p[[i]][s,] <- softmax2( p[[i]][s,] )
  }
  return(p)
}

################################
###Make Mammal Density Figure###
################################
dens<- data.frame(
  id = 0 ,
  society=0,
  sex = 0,
  middle = 1,
  ado = 0,
  water=1,
  dens=c(0,1),
  nonforaged=mean(d$nonforaged_z)
)

seq.length<-2
p_dens <- link.TA4 (dens)
p_dens_mean <- sapply( 1:length(p_dens) , function(i) apply(p_dens[[i]],2,mean) )

labels.at <- c(1,2,3,4)

##Make figure
quartz(height = 6 , width = 18)
par(mfrow=c(1,2), mar=c(1,2,1,2) , oma=c(5,2.5,2.5,2.5))

##
plot( NULL , xlim=c(0.5,4.5) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9, cex.axis=1.2)
for ( k in 1:1 ) {
  p_dens_PI <- sapply( 1:length(p_dens) , function(i) PI(p_dens[[i]][,k],prob=0.89) ) 
  segments(0.9,p_dens_PI[1,1],0.9,p_dens_PI[2,1],col=col.alpha("red",1),lwd=3)
  segments(1.1,p_dens_PI[1,2],1.1,p_dens_PI[2,2],col=col.alpha("purple",1),lwd=3)
  
  points( c(0.9,1.1) , p_dens_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

##B Food Production
for ( k in 2:2) {
  p_dens_PI <- sapply( 1:length(p_dens) , function(i) PI(p_dens[[i]][,k],prob=0.89) ) 
  segments(1.9,p_dens_PI[1,1],1.9,p_dens_PI[2,1],col=col.alpha("red",1),lwd=3)
  segments(2.1,p_dens_PI[1,2],2.1,p_dens_PI[2,2],col=col.alpha("purple",1),lwd=3)
  
  points( c(1.9,2.1) , p_dens_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 3:3) {
  p_dens_PI <- sapply( 1:length(p_dens) , function(i) PI(p_dens[[i]][,k],prob=0.89) ) 
  segments(2.9,p_dens_PI[1,1],2.9,p_dens_PI[2,1],col=col.alpha("red",1),lwd=3)
  segments(3.1,p_dens_PI[1,2],3.1,p_dens_PI[2,2],col=col.alpha("purple",1),lwd=3)
  
  points( c(2.9,3.1) , p_dens_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 4:4) {
  p_dens_PI <- sapply( 1:length(p_dens) , function(i) PI(p_dens[[i]][,k],prob=0.89) ) 
  segments(3.9,p_dens_PI[1,1],3.9,p_dens_PI[2,1],col=col.alpha("red",1),lwd=3)
  segments(4.1,p_dens_PI[1,2],4.1,p_dens_PI[2,2],col=col.alpha("purple",1),lwd=3)
  
  points( c(3.9,4.1) , p_dens_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}


axis( side = 1, at = labels.at, labels = c("Childcare","Food Production","Domestic Work","Play"), cex.axis=1.3)
mtext(text="Probability (Activity)",side=2,line=2, cex=1.4)
mtext(text="(A) Activity Participation by Dangerous Mammal Density",side=1,line=3,cex=1.4)

legend("topleft", 
       inset=0, 
       cex = 1.2, 
       c("<1 n/km2",">10 n/km2"), 
       horiz=FALSE, 
       lty=c(1,1), 
       lwd=c(2,2), 
       col=c("red","purple"), 
       bg="",
       bty="n",
       title="Dangerous Mammal Density",
       text.font=1.3)

########################################
###Make Water Quality/Quantity Figure###
########################################

Water<- data.frame(
  id = 0 ,
  society=0,
  sex = 0,
  middle = 1,
  ado = 0,
  water=c(1,2,4),
  dens=0,
  nonforaged=mean(d$nonforaged_z)
)

seq.length<-3
p_Water<-link.TA4 (Water)
p_Water_mean <- sapply( 1:length(p_Water) , function(i) apply(p_Water[[i]],2,mean) )

##Make Figure for effect of water on children's activities
labels.at <- c(1,2,3,4)

plot( NULL , xlim=c(0.5,4.5) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9, cex.axis=1, col.axis = "white")
for ( k in 1:1 ) {
  p_Water_PI <- sapply( 1:length(p_Water) , function(i) PI(p_Water[[i]][,k],prob=0.89) ) 
  segments(0.75,p_Water_PI[1,1],0.75,p_Water_PI[2,1],col=col.alpha("lightblue",1),lwd=3)
  segments(1,p_Water_PI[1,2],1,p_Water_PI[2,2],col=col.alpha("navy",1),lwd=3)
  segments(1.25,p_Water_PI[1,3],1.25,p_Water_PI[2,3],col=col.alpha("brown",1),lwd=3)
  
  points( c(0.75,1,1.25) , p_Water_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 2:2 ) {
  p_Water_PI <- sapply( 1:length(p_Water) , function(i) PI(p_Water[[i]][,k],prob=0.89) ) 
  segments(1.75,p_Water_PI[1,1],1.75,p_Water_PI[2,1],col=col.alpha("lightblue",1),lwd=3)
  segments(2,p_Water_PI[1,2],2,p_Water_PI[2,2],col=col.alpha("navy",1),lwd=3)
  segments(2.25,p_Water_PI[1,3],2.25,p_Water_PI[2,3],col=col.alpha("brown",1),lwd=3)
  
  points( c(1.75,2,2.25) , p_Water_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 3:3 ) {
  p_Water_PI <- sapply( 1:length(p_Water) , function(i) PI(p_Water[[i]][,k],prob=0.89) ) 
  segments(2.75,p_Water_PI[1,1],2.75,p_Water_PI[2,1],col=col.alpha("lightblue",1),lwd=3)
  segments(3,p_Water_PI[1,2],3,p_Water_PI[2,2],col=col.alpha("navy",1),lwd=3)
  segments(3.25,p_Water_PI[1,3],3.25,p_Water_PI[2,3],col=col.alpha("brown",1),lwd=3)
  
  points( c(2.75,3,3.25) , p_Water_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 4:4 ) {
  p_Water_PI <- sapply( 1:length(p_Water) , function(i) PI(p_Water[[i]][,k],prob=0.89) ) 
  segments(3.75,p_Water_PI[1,1],3.75,p_Water_PI[2,1],col=col.alpha("lightblue",1),lwd=3)
  segments(4,p_Water_PI[1,2],4,p_Water_PI[2,2],col=col.alpha("navy",1),lwd=3)
  segments(4.25,p_Water_PI[1,3],4.25,p_Water_PI[2,3],col=col.alpha("brown",1),lwd=3)
  
  points( c(3.75,4,4.25) , p_Water_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}


axis( side = 1, at = labels.at, labels = c("Childcare","Food Production","Domestic Work","Play"), cex.axis=1.3)

mtext(text="(B) Activity Participation by Water Quality/Quantity",side=1,line=3,cex=1.4)
legend("topleft", 
       inset=0, 
       cex = 1.2, 
       c("High Quality/Quantity","Low Quality/High Quantity","Low Quality/Quantity"), 
       horiz=FALSE, 
       lty=c(1,1), 
       lwd=c(2,2), 
       col=c("lightblue","navy","brown"), 
       bg="",
       bty="n",
       title="Water",
       text.font=1.3)

quartz.save(file="fig3.pdf",type="pdf")


#####################################################################
###Same figure but with the Other Activities category for the supp###
#####################################################################

labels.at <- c(1,2,3,4,5)

##Make figure
quartz(height = 6 , width = 11)
par(mfrow=c(1,1), mar=c(1,2,1,2) + 2 , oma=c(2.5,4,2.5,2.5))

##
plot( NULL , xlim=c(0.5,5.5) , ylim=c(0,0.8) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 1:1 ) {
  p_dens_PI <- sapply( 1:length(p_dens) , function(i) PI(p_dens[[i]][,k],prob=0.89) ) 
  segments(0.9,p_dens_PI[1,1],0.9,p_dens_PI[2,1],col=col.alpha("red",1),lwd=3)
  segments(1.1,p_dens_PI[1,2],1.1,p_dens_PI[2,2],col=col.alpha("purple",1),lwd=3)
  
  points( c(0.9,1.1) , p_dens_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

##B Food Production
for ( k in 2:2) {
  p_dens_PI <- sapply( 1:length(p_dens) , function(i) PI(p_dens[[i]][,k],prob=0.89) ) 
  segments(1.9,p_dens_PI[1,1],1.9,p_dens_PI[2,1],col=col.alpha("red",1),lwd=3)
  segments(2.1,p_dens_PI[1,2],2.1,p_dens_PI[2,2],col=col.alpha("purple",1),lwd=3)
  
  points( c(1.9,2.1) , p_dens_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 3:3) {
  p_dens_PI <- sapply( 1:length(p_dens) , function(i) PI(p_dens[[i]][,k],prob=0.89) ) 
  segments(2.9,p_dens_PI[1,1],2.9,p_dens_PI[2,1],col=col.alpha("red",1),lwd=3)
  segments(3.1,p_dens_PI[1,2],3.1,p_dens_PI[2,2],col=col.alpha("purple",1),lwd=3)
  
  points( c(2.9,3.1) , p_dens_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 4:4) {
  p_dens_PI <- sapply( 1:length(p_dens) , function(i) PI(p_dens[[i]][,k],prob=0.89) ) 
  segments(3.9,p_dens_PI[1,1],3.9,p_dens_PI[2,1],col=col.alpha("red",1),lwd=3)
  segments(4.1,p_dens_PI[1,2],4.1,p_dens_PI[2,2],col=col.alpha("purple",1),lwd=3)
  
  points( c(3.9,4.1) , p_dens_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 5:5) {
 p_dens_PI <- sapply( 1:length(p_dens) , function(i) PI(p_dens[[i]][,k],prob=0.89) ) 
 segments(4.75,p_dens_PI[1,1],4.75,p_dens_PI[2,1],col=col.alpha("red",1),lwd=3)
 segments(5.25,p_dens_PI[1,2],5.25,p_dens_PI[2,2],col=col.alpha("purple",1),lwd=3)

 points( c(4.75,5.25) , p_dens_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}



axis( side = 1, at = labels.at, labels = c("Childcare","Food Production","Domestic Work","Play","Other Activities"))
mtext(text="Probability (Activity)",side=2,line=2)
mtext(text="(A) Activity Participation by Dangerous Mammal Density",side=1,line=3,font=2)
legend("topleft", 
       inset=0, 
       cex = 0.75, 
       c("<1",">10"), 
       horiz=FALSE, 
       lty=c(1,1), 
       lwd=c(2,2), 
       col=c("red","purple"), 
       bg="grey96",
       title="Dangerous Mammal Density (n/km2)",
       text.font=0.5)

quartz.save(file="figS6a.pdf",type="pdf")
rm(dens)
rm(p_dens)
rm(p_dens_mean)
rm(p_dens_PI)
rm(k)
rm(labels.at)
rm(seq.length)
dev.off()


#######################################################################
###Make same figure but for the supplement with the Other Activities###
#######################################################################
labels.at <- c(1,2,3,4,5)

quartz(height = 6 , width =11)
par(mfrow=c(1,1), mar=c(1,2,1,2)+2 , oma=c(5,2.5,2.5,2.5))
plot( NULL , xlim=c(0.5,5.5) , ylim=c(0,0.8) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 1:1 ) {
  p_Water_PI <- sapply( 1:length(p_Water) , function(i) PI(p_Water[[i]][,k],prob=0.89) ) 
  segments(0.75,p_Water_PI[1,1],0.75,p_Water_PI[2,1],col=col.alpha("lightblue",1),lwd=3)
  segments(1,p_Water_PI[1,2],1,p_Water_PI[2,2],col=col.alpha("navy",1),lwd=3)
  segments(1.25,p_Water_PI[1,3],1.25,p_Water_PI[2,3],col=col.alpha("brown",1),lwd=3)
  
  points( c(0.75,1,1.25) , p_Water_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 2:2 ) {
  p_Water_PI <- sapply( 1:length(p_Water) , function(i) PI(p_Water[[i]][,k],prob=0.89) ) 
  segments(1.75,p_Water_PI[1,1],1.75,p_Water_PI[2,1],col=col.alpha("lightblue",1),lwd=3)
  segments(2,p_Water_PI[1,2],2,p_Water_PI[2,2],col=col.alpha("navy",1),lwd=3)
  segments(2.25,p_Water_PI[1,3],2.25,p_Water_PI[2,3],col=col.alpha("brown",1),lwd=3)
  
  points( c(1.75,2,2.25) , p_Water_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 3:3 ) {
  p_Water_PI <- sapply( 1:length(p_Water) , function(i) PI(p_Water[[i]][,k],prob=0.89) ) 
  segments(2.75,p_Water_PI[1,1],2.75,p_Water_PI[2,1],col=col.alpha("lightblue",1),lwd=3)
  segments(3,p_Water_PI[1,2],3,p_Water_PI[2,2],col=col.alpha("navy",1),lwd=3)
  segments(3.25,p_Water_PI[1,3],3.25,p_Water_PI[2,3],col=col.alpha("brown",1),lwd=3)
  
  points( c(2.75,3,3.25) , p_Water_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 4:4 ) {
  p_Water_PI <- sapply( 1:length(p_Water) , function(i) PI(p_Water[[i]][,k],prob=0.89) ) 
  segments(3.75,p_Water_PI[1,1],3.75,p_Water_PI[2,1],col=col.alpha("lightblue",1),lwd=3)
  segments(4,p_Water_PI[1,2],4,p_Water_PI[2,2],col=col.alpha("navy",1),lwd=3)
  segments(4.25,p_Water_PI[1,3],4.25,p_Water_PI[2,3],col=col.alpha("brown",1),lwd=3)
  
  points( c(3.75,4,4.25) , p_Water_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

for ( k in 5:5 ) {
  p_Water_PI <- sapply( 1:length(p_Water) , function(i) PI(p_Water[[i]][,k],prob=0.89) ) 
  segments(4.75,p_Water_PI[1,1],4.75,p_Water_PI[2,1],col=col.alpha("lightblue",1),lwd=3)
  segments(5,p_Water_PI[1,2],5,p_Water_PI[2,2],col=col.alpha("navy",1),lwd=3)
  segments(5.25,p_Water_PI[1,3],5.25,p_Water_PI[2,3],col=col.alpha("brown",1),lwd=3)

  points( c(4.75,5,5.25) , p_Water_mean[k,], pch = 21, bg = "white",col=col.alpha("black",1))
}

axis( side = 1, at = labels.at, labels = c("Childcare","Food Production","Domestic Work","Play","Other Activities"))
mtext(text="Probability (Activity)",side=2,line=2)
mtext(text="(B) Activity Participation by Water Quality/Quantity",side=1,line=3,font=2)
legend("topleft", 
       inset=0, 
       cex = 0.75, 
       c("High Quality/Quantity","Low Quality/High Quantity","Low Quality/Quantity"), 
       horiz=FALSE, 
       lty=c(1,1), 
       lwd=c(2,2), 
       col=c("lightblue","navy","brown"), 
       bg="grey96",
       title="Water",
       text.font=0.5)

quartz.save(file="figS6b.pdf",type="pdf")
dev.off()
rm(p_Water)
rm(p_Water_mean)
rm(p_Water_PI)
rm(Water)
rm(k)
rm(labels.at)
rm(link.TA4)
dev.off()

########################################
###Generate predictions for Model 3.5###
########################################
link.TA3.5 <- function( data ) {
  K <- dim(post3.5$v_id)[3] + 1
  ns <- dim(post3.5$v_id)[1]
  if ( missing(data) ) stop( "BOOM: Need data argument" )
  n <- seq.length
  
  softmax2 <- function(x) {
    x <- max(x) - x
    exp(-x)/sum(exp(-x))
  }
  
  p <- list()
  
  for ( i in 1:n ) {
    p[[i]] <- sapply( 1:K , function(k) {
      if ( k < K ) {
        ptemp <- post3.5$a[,k] + 
          post3.5$b_sex[,k] * data$sex[i] + 
          post3.5$b_middle[,k] * data$middle[i] + 
          post3.5$b_ado[,k] * data$ado[i] + 
          post3.5$b_sex_middle[,k] * data$sex[i] * data$middle[i] + 
          post3.5$b_sex_ado[,k] * data$sex[i] * data$ado[i] + 
          post3.5$b_nonforaged[,k] * data$nonforaged[i]+
          post3.5$b_NPP[,k] * data$NPP[i]+
          post3.5$b_temp[,k] * data$temp[i]+
          post3.5$b_prec[,k] * data$prec[i]
        if ( data$id[i]>0 ) ptemp <- ptemp + post3.5$v_id[,data$id[i],k]
        if ( data$society[i]>0 ) ptemp <- ptemp + post3.5$v_society[,data$society[i],k]
      } else {
        ptemp <- rep(0,ns)
      }
      return(ptemp)
    })
    ## The values are converted to probabilities using the softmax function
    ## which ensures that the predicted values across categories sum to
    ## 100% probabilities.
    for ( s in 1:ns ) p[[i]][s,] <- softmax2( p[[i]][s,] )
  }
  return(p)
}

####################################
###Make temp Figure without Dukha###
####################################
seq.length <- 100
for(t in unique(ds$variable)) {ds[paste("r_",t,sep="")] <- ifelse(ds$variable==t,1,0)}
temp_seq_ds <- seq(from= min(ds$temp_z) , to= max(ds$temp_z), length.out = seq.length)
agg_temp_ds <- aggregate ( cbind (r_food_production, r_znonworkobs, r_play, r_childcare, r_household) ~ Society+ temp_z, data = ds, FUN = mean)

temp_ds <- data.frame(
  id = 0 ,
  society=0,
  sex = 0,
  middle = 0,
  ado = 0,
  nonforaged=mean(ds$nonforaged_z),
  NPP=mean(ds$NPP_z),
  temp=temp_seq_ds,
  prec=mean(ds$prec_z)
)

p_temp_ds <- link.TA3.5 (temp_ds)
p_temp_mean_ds <- sapply( 1:length(p_temp_ds) , function(i) apply(p_temp_ds[[i]],2,mean) )
preferred.temp<-c(17,20,23,26)
labels.at <- (preferred.temp - mean(ds$meanAnnualTemp))/sd(ds$meanAnnualTemp)

quartz(height = 6 , width = 10)
par(mfrow=c(2,3), mar=c(1,2,1,2) + 0.1 , oma=c(2.5,4,2.5,2.5))
plot( NULL , xlim=c(min(ds$temp_z),max(ds$temp_z)) , ylim=c(0,1) , xaxt = "n", main = " ",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 1:1 ) {
  lines( temp_seq_ds , p_temp_mean_ds[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_ds_PI <- sapply( 1:length(p_temp_ds) , function(i) PI(p_temp_ds[[i]][,k],prob=0.89) )
  shade( p_temp_ds_PI , temp_seq_ds,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = F)
points(agg_temp_ds$temp_z, agg_temp_ds$r_childcare)
mtext(text="Probability (Childcare)",side=2,line=2)

plot( NULL , xlim=c(min(ds$temp_z),max(ds$temp_z)) , ylim=c(0,1) , xaxt = "n", main = " ",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 2:2 ) {
  lines( temp_seq_ds , p_temp_mean_ds[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_ds_PI <- sapply( 1:length(p_temp_ds) , function(i) PI(p_temp_ds[[i]][,k],prob=0.89) )
  shade( p_temp_ds_PI , temp_seq_ds,col=col.alpha("black",0.1))
}
points(agg_temp_ds$temp_z, agg_temp_ds$r_food_production)

axis( side = 1, at = labels.at, labels = F)
mtext(text="Probability (Food Production)",side=2,line=2)

plot( NULL , xlim=c(min(ds$temp_z),max(ds$temp_z)) , ylim=c(0,1) , xaxt = "n", main = " ",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 3:3 ) {
  lines( temp_seq_ds , p_temp_mean_ds[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_ds_PI <- sapply( 1:length(p_temp_ds) , function(i) PI(p_temp_ds[[i]][,k],prob=0.89) )
  shade( p_temp_ds_PI , temp_seq_ds,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.temp)
points(agg_temp_ds$temp_z, agg_temp_ds$r_household)
mtext(text="Probability (Domestic Work)",side=2,line=2)

plot( NULL , xlim=c(min(ds$temp_z),max(ds$temp_z)) , ylim=c(0,1) , xaxt = "n", main = " ",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 4:4 ) {
  lines( temp_seq_ds , p_temp_mean_ds[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_ds_PI <- sapply( 1:length(p_temp_ds) , function(i) PI(p_temp_ds[[i]][,k],prob=0.89) )
  shade( p_temp_ds_PI , temp_seq_ds,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.temp)
points(agg_temp_ds$temp_z, agg_temp_ds$r_play)
mtext(text="Probability (Play)",side=2,line=2)

plot( NULL , xlim=c(min(ds$temp_z),max(ds$temp_z)) , ylim=c(0,1) , xaxt = "n", main = " ",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 5:5 ) {
  lines( temp_seq_ds , p_temp_mean_ds[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_ds_PI <- sapply( 1:length(p_temp_ds) , function(i) PI(p_temp_ds[[i]][,k],prob=0.89) )
  shade( p_temp_ds_PI , temp_seq_ds,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.temp)
points(agg_temp_ds$temp_z, agg_temp_ds$r_znonworkobs)
mtext(text="Probability (Other Activities)",side=2,line=2)

mtext(text="Mean Annual Temperature (°C)",side=1,line=1,outer=TRUE)

quartz.save(file="figs8.pdf",type="pdf")
dev.off()
rm(agg_temp_ds)
rm(p_temp_ds)
rm(p_temp_ds_PI)
rm(p_temp_mean_ds)
rm(temp_ds)
rm(k)
rm(labels.at)
rm(preferred.temp)
rm(seq.length)
rm(t)
rm(temp_seq_ds)
rm(link.TA3.5)

########################################
###Generate predictions for Model 3.4###
########################################
link.TA3.4 <- function( data ) {
  K <- dim(post3.4$v_id)[3] + 1
  ns <- dim(post3.4$v_id)[1]
  if ( missing(data) ) stop( "BOOM: Need data argument" )
  n <- seq.length
  
  softmax2 <- function(x) {
    x <- max(x) - x
    exp(-x)/sum(exp(-x))
  }
  
  p <- list()
  
  for ( i in 1:n ) {
    p[[i]] <- sapply( 1:K , function(k) {
      if ( k < K ) {
        ptemp <- post3.4$a[,k] + 
          post3.4$b_sex[,k] * data$sex[i] + 
          post3.4$b_middle[,k] * data$middle[i] + 
          post3.4$b_ado[,k] * data$ado[i] + 
          post3.4$b_sex_middle[,k] * data$sex[i] * data$middle[i] + 
          post3.4$b_sex_ado[,k] * data$sex[i] * data$ado[i] + 
          post3.4$b_nonforaged[,k] * data$nonforaged[i]+
          post3.4$b_CV[,k] * data$CV[i]
        if ( data$id[i]>0 ) ptemp <- ptemp + post3.4$v_id[,data$id[i],k]
        if ( data$society[i]>0 ) ptemp <- ptemp + post3.4$v_society[,data$society[i],k]
      } else {
        ptemp <- rep(0,ns)
      }
      return(ptemp)
    })
    ## The values are converted to probabilities using the softmax function
    ## which ensures that the predicted values across categories sum to
    ## 100% probabilities.
    for ( s in 1:ns ) p[[i]][s,] <- softmax2( p[[i]][s,] )
  }
  return(p)
}


####################
###Make CV Figure###
####################
seq.length <- 100
CV_seq <- seq(from= min(d$CV_z) , to= max(d$CV_z), length.out = seq.length)
for(t in unique(d$variable)) {d[paste("r_",t,sep="")] <- ifelse(d$variable==t,1,0)}
agg_CV <- aggregate ( cbind (r_food_production, r_znonworkobs, r_play, r_childcare, r_household) ~ Society+ CV_z, data = d, FUN = mean)

CV_d <- data.frame(
  id = 0 ,
  society=0,
  sex = 0,
  middle = 0,
  ado = 0,
  nonforaged=mean(d$nonforaged_z),
  CV=CV_seq
)


p_CV <- link.TA3.4 (CV_d)
p_CV_mean <- sapply( 1:length(p_CV) , function(i) apply(p_CV[[i]],2,mean) )
preferred.CV<-c(50,75,100,125)
labels.at <- (preferred.CV - mean(d$precSeasonality))/sd(d$precSeasonality)

quartz(height = 6 , width = 10)
par(mfrow=c(2,3), mar=c(1,2,1,2) + 0.1 , oma=c(2.5,4,2.5,2.5))
plot( NULL , xlim=c(min(d$CV_z),max(d$CV_z)) , ylim=c(0,1) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 1:1 ) {
  lines( CV_seq , p_CV_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_CV_PI <- sapply( 1:length(p_CV) , function(i) PI(p_CV[[i]][,k],prob=0.89) )
  shade( p_CV_PI , CV_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = F)
points(agg_CV$CV_z, agg_CV$r_childcare)
mtext(text="Probability (Childcare)",side=2,line=2)

plot( NULL , xlim=c(min(d$CV_z),max(d$CV_z)) , ylim=c(0,1) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 2:2 ) {
  lines( CV_seq , p_CV_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_CV_PI <- sapply( 1:length(p_CV) , function(i) PI(p_CV[[i]][,k],prob=0.89) )
  shade( p_CV_PI , CV_seq,col=col.alpha("black",0.1))
}
points(agg_CV$CV_z, agg_CV$r_food_production)

axis( side = 1, at = labels.at, labels = F)
mtext(text="Probability (Food Production)",side=2,line=2)

plot( NULL , xlim=c(min(d$CV_z),max(d$CV_z)) , ylim=c(0,1) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 3:3 ) {
  lines( CV_seq , p_CV_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_CV_PI <- sapply( 1:length(p_CV) , function(i) PI(p_CV[[i]][,k],prob=0.89) )
  shade( p_CV_PI , CV_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.CV)
points(agg_CV$CV_z, agg_CV$r_household)
mtext(text="Probability (Domestic Work)",side=2,line=2)

plot( NULL , xlim=c(min(d$CV_z),max(d$CV_z)) , ylim=c(0,1) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 4:4 ) {
  lines( CV_seq , p_CV_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_CV_PI <- sapply( 1:length(p_CV) , function(i) PI(p_CV[[i]][,k],prob=0.89) )
  shade( p_CV_PI , CV_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.CV)
points(agg_CV$CV_z, agg_CV$r_play)
mtext(text="Probability (Play)",side=2,line=2)

plot( NULL , xlim=c(min(d$CV_z),max(d$CV_z)) , ylim=c(0,1) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 5:5 ) {
  lines( CV_seq , p_CV_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_CV_PI <- sapply( 1:length(p_CV) , function(i) PI(p_CV[[i]][,k],prob=0.89) )
  shade( p_CV_PI , CV_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.CV)
points(agg_CV$CV_z, agg_CV$r_znonworkobs)
mtext(text="Probability (Other Activities)",side=2,line=2)


mtext(text="Coefficient of Variation of Monthly Precipitation (%) ",side=1,line=1,outer=TRUE)

quartz.save(file="figs9.pdf",type="pdf")
dev.off()
rm(agg_CV)
rm(CV_d)
rm(p_CV)
rm(p_CV_PI)
rm(p_CV_mean)
rm(k)
rm(labels.at)
rm(preferred.CV)
rm(seq.length)
rm(t)
rm(CV_seq)
rm(link.TA3.4)

########################################
###Generate predictions for Model 4.3###
########################################
link.TA4.3 <- function( data ) {
  K <- dim(post4.3$v_id)[3] + 1
  ns <- dim(post4.3$v_id)[1]
  if ( missing(data) ) stop( "BOOM: Need data argument" )
  n <- seq.length
  
  softmax2 <- function(x) {
    x <- max(x) - x
    exp(-x)/sum(exp(-x))
  }
  
  p <- list()
  
  for ( i in 1:n ) {
    p[[i]] <- sapply( 1:K , function(k) {
      if ( k < K ) {
        ptemp <- post4.3$a[,k] + 
          post4.3$b_sex[,k] * data$sex[i] + 
          post4.3$b_middle[,k] * data$middle[i] + 
          post4.3$b_ado[,k] * data$ado[i] + 
          post4.3$b_sex_middle[,k] * data$sex[i] * data$middle[i] + 
          post4.3$b_sex_ado[,k] * data$sex[i] * data$ado[i] + 
          post4.3$b_nonforaged[,k] * data$nonforaged[i]+
          post4.3$b_snake[,k] * data$snake[i]
        if ( data$id[i]>0 ) ptemp <- ptemp + post4.3$v_id[,data$id[i],k]
        if ( data$society[i]>0 ) ptemp <- ptemp + post4.3$v_society[,data$society[i],k]
      } else {
        ptemp <- rep(0,ns)
      }
      return(ptemp)
    })
    ## The values are converted to probabilities using the softmax function
    ## which ensures that the predicted values across categories sum to
    ## 100% probabilities.
    for ( s in 1:ns ) p[[i]][s,] <- softmax2( p[[i]][s,] )
  }
  return(p)
}
#############################
###Make Snake Count Figure###
#############################

seq.length <- 100
snake_seq <- seq(from= min(d$snake_z) , to= max(d$snake_z), length.out = seq.length)

snake<- data.frame(
  id = 0 ,
  society=0,
  sex = 0,
  middle = 0,
  ado = 0,
  snake=snake_seq,
  nonforaged=mean(d$nonforaged_z)
)

p_snake <- link.TA4.3 (snake)
p_snake_mean <- sapply( 1:length(p_snake) , function(i) apply(p_snake[[i]],2,mean) )
preferred.snake<-c(0,3,6,9,12)
labels.at <- (preferred.snake - mean(d$snake_count))/sd(d$snake_count)

##Aggregate raw data
for(t in unique(d$variable)) {d[paste("r_",t,sep="")] <- ifelse(d$variable==t,1,0)}
agg_snake <- aggregate ( cbind (r_food_production, r_znonworkobs, r_play, r_childcare, r_household) ~ Society+ snake_z, data = d, FUN = mean)

##Make figure
quartz(height = 6 , width = 10)
par(mfrow=c(2,3), mar=c(1,2,1,2) + 0.1 , oma=c(2.5,4,2.5,2.5))

##A Childcare
plot( NULL , xlim=c(min(d$snake_z),max(d$snake_z)) , ylim=c(0,0.8) , xaxt = "n", main = " ",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 1:1 ) {
  lines( snake_seq , p_snake_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_snake_PI <- sapply( 1:length(p_snake) , function(i) PI(p_snake[[i]][,k],prob=0.89) )
  shade( p_snake_PI , snake_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = F)
points(agg_snake$snake_z, agg_snake$r_childcare)
mtext(text="Probability (Childcare)",side=2,line=2)

##B Food Production
plot( NULL , xlim=c(min(d$snake_z),max(d$snake_z)) , ylim=c(0,0.8) , xaxt = "n", main = " ",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 2:2) {
  lines( snake_seq , p_snake_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_snake_PI <- sapply( 1:length(p_snake) , function(i) PI(p_snake[[i]][,k],prob=0.89) )
  shade( p_snake_PI , snake_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = F)
points(agg_snake$snake_z, agg_snake$r_food_production)
mtext(text="Probability (Food Production)",side=2,line=2)

##C Domestic Work
plot( NULL , xlim=c(min(d$snake_z),max(d$snake_z)) , ylim=c(0,0.8) , xaxt = "n", main = " ",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 3:3) {
  lines( snake_seq , p_snake_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_snake_PI <- sapply( 1:length(p_snake) , function(i) PI(p_snake[[i]][,k],prob=0.89) )
  shade( p_snake_PI , snake_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.snake)
points(agg_snake$snake_z, agg_snake$r_household)
mtext(text="Probability (Domestic Work)",side=2,line=2)

##D Play
plot( NULL , xlim=c(min(d$snake_z),max(d$snake_z)) , ylim=c(0,0.8) , xaxt = "n", main = " ",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 4:4) {
  lines( snake_seq , p_snake_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_snake_PI <- sapply( 1:length(p_snake) , function(i) PI(p_snake[[i]][,k],prob=0.89) )
  shade( p_snake_PI , snake_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.snake)
points(agg_snake$snake_z, agg_snake$r_play)
mtext(text="Probability (Play)",side=2,line=2)

plot( NULL , xlim=c(min(d$snake_z),max(d$snake_z)) , ylim=c(0,0.8) , xaxt = "n", main = " ",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 5:5) {
  lines( snake_seq , p_snake_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_snake_PI <- sapply( 1:length(p_snake) , function(i) PI(p_snake[[i]][,k],prob=0.89) )
  shade( p_snake_PI , snake_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.snake)
points(agg_snake$snake_z, agg_snake$r_znonworkobs)
mtext(text="Probability (Other Activities)",side=2,line=2)

mtext(text="Medically Important Venomous Snake Species (count)",side=1,line=1,outer=TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

quartz.save(file="figS10.pdf",type="pdf")
rm(agg_snake)
rm(p_snake)
rm(p_snake_mean)
rm(p_snake_PI)
rm(snake)
rm(k)
rm(labels.at)
rm(preferred.snake)
rm(seq.length)
rm(snake_seq)
rm(t)
rm(link.TA4.3)
dev.off()

########################################
###Generate predictions for Model 3.6###
########################################

link.TA3.6 <- function( data ) {
  K <- dim(post3.6$v_id)[3] + 1
  ns <- dim(post3.6$v_id)[1]
  if ( missing(data) ) stop( "BOOM: Need data argument" )
  n <- seq.length
  
  softmax2 <- function(x) {
    x <- max(x) - x
    exp(-x)/sum(exp(-x))
  }
  
  p <- list()
  
  for ( i in 1:n ) {
    p[[i]] <- sapply( 1:K , function(k) {
      if ( k < K ) {
        ptemp <- post3.6$a[,k] + 
          post3.6$b_sex[,k] * data$sex[i] + 
          post3.6$b_middle[,k] * data$middle[i] + 
          post3.6$b_sex_middle[,k] * data$sex[i] * data$middle[i] + 
          post3.6$b_nonforaged[,k] * data$nonforaged[i]+
          post3.6$b_NPP[,k] * data$NPP[i]+
          post3.6$b_temp[,k] * data$temp[i]+
          post3.6$b_prec[,k] * data$prec[i]
        if ( data$id[i]>0 ) ptemp <- ptemp + post3.6$v_id[,data$id[i],k]
        if ( data$society[i]>0 ) ptemp <- ptemp + post3.6$v_society[,data$society[i],k]
      } else {
        ptemp <- rep(0,ns)
      }
      return(ptemp)
    })
    ## The values are converted to probabilities using the softmax function
    ## which ensures that the predicted values across categories sum to
    ## 100% probabilities.
    for ( s in 1:ns ) p[[i]][s,] <- softmax2( p[[i]][s,] )
  }
  return(p)
}

#####################
###Make NPP Figure###
#####################
seq.length <- 100

NPP_seq <- seq(from= min(dc$NPP_z) , to= max(dc$NPP_z), length.out = seq.length)
agg_NPP <- aggregate ( cbind (r_food_production, r_znonworkobs, r_play, r_childcare, r_household) ~ Society+ NPP_z, data = dc, FUN = mean)

NPP <- data.frame(
  id = 0 ,
  society=0,
  sex = 0,
  middle = 0,
  nonforaged=mean(dc$nonforaged_z),
  NPP=NPP_seq,
  temp=mean(dc$temp_z),
  prec=mean(dc$prec_z)
)

p_NPP <- link.TA3.6 (NPP)
p_NPP_mean <- sapply( 1:length(p_NPP) , function(i) apply(p_NPP[[i]],2,mean) )
preferred.NPP<-c(145,645,1145,1645,2145,2645)
labels.at <- (preferred.NPP - mean(dc$NPP))/sd(dc$NPP)

quartz(height = 6 , width = 10)
par(mfrow=c(2,3), mar=c(1,2,1,2) + 0.1 , oma=c(2.5,4,2.5,2.5))
plot( NULL , xlim=c(min(dc$NPP_z),max(dc$NPP_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 1:1 ) {
  lines( NPP_seq , p_NPP_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_NPP_PI <- sapply( 1:length(p_NPP) , function(i) PI(p_NPP[[i]][,k],prob=0.89) )
  shade( p_NPP_PI , NPP_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = F)
points(agg_NPP$NPP_z, agg_NPP$r_childcare)
mtext(text="Probability (Childcare)",side=2,line=2)

plot( NULL , xlim=c(min(dc$NPP_z),max(dc$NPP_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 2:2 ) {
  lines( NPP_seq , p_NPP_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_NPP_PI <- sapply( 1:length(p_NPP) , function(i) PI(p_NPP[[i]][,k],prob=0.89) )
  shade( p_NPP_PI , NPP_seq,col=col.alpha("black",0.1))
}
points(agg_NPP$NPP_z, agg_NPP$r_food_production)

axis( side = 1, at = labels.at, labels = F)
mtext(text="Probability (Food Production)",side=2,line=2)

plot( NULL , xlim=c(min(dc$NPP_z),max(dc$NPP_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 3:3 ) {
  lines( NPP_seq , p_NPP_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_NPP_PI <- sapply( 1:length(p_NPP) , function(i) PI(p_NPP[[i]][,k],prob=0.89) )
  shade( p_NPP_PI , NPP_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.NPP)
points(agg_NPP$NPP_z, agg_NPP$r_household)
mtext(text="Probability (Domestic Work)",side=2,line=2)

plot( NULL , xlim=c(min(dc$NPP_z),max(dc$NPP_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 4:4 ) {
  lines( NPP_seq , p_NPP_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_NPP_PI <- sapply( 1:length(p_NPP) , function(i) PI(p_NPP[[i]][,k],prob=0.89) )
  shade( p_NPP_PI , NPP_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.NPP)
points(agg_NPP$NPP_z, agg_NPP$r_play)
mtext(text="Probability (Play)",side=2,line=2)

plot( NULL , xlim=c(min(dc$NPP_z),max(dc$NPP_z)) , ylim=c(0,0.8) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 5:5 ) {
  lines( NPP_seq , p_NPP_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_NPP_PI <- sapply( 1:length(p_NPP) , function(i) PI(p_NPP[[i]][,k],prob=0.89) )
  shade( p_NPP_PI , NPP_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.NPP)
points(agg_NPP$NPP_z, agg_NPP$r_znonworkobs)
mtext(text="Probability (Other Activities)",side=2,line=2)

mtext(text="Net Primary Productivity (NPP)",side=1,line=1,outer=TRUE)
quartz.save(file="figS11.pdf",type="pdf")
rm(agg_NPP)
rm(NPP)
rm(p_NPP)
rm(p_NPP_mean)
rm(p_NPP_PI)
rm(k)
rm(labels.at)
rm(NPP_seq)
rm(preferred.NPP)
rm(seq.length)
dev.off()

######################
###Make temp Figure###
######################
seq.length <- 100
temp_seq <- seq(from= min(dc$temp_z) , to= max(dc$temp_z), length.out = seq.length)
agg_temp <- aggregate ( cbind (r_food_production, r_znonworkobs, r_play, r_childcare, r_household) ~ Society+ temp_z, data = dc, FUN = mean)

temp <- data.frame(
  id = 0 ,
  society=0,
  sex = 0,
  middle = 0,
  nonforaged=mean(dc$nonforaged_z),
  NPP=mean(dc$NPP_z),
  temp=temp_seq,
  prec=mean(dc$prec_z)
)

p_temp <- link.TA3.6 (temp)
p_temp_mean <- sapply( 1:length(p_temp) , function(i) apply(p_temp[[i]],2,mean) )
preferred.temp<-c(-5,5,15,25)
labels.at <- (preferred.temp - mean(dc$meanAnnualTemp))/sd(dc$meanAnnualTemp)

quartz(height = 6 , width = 10)
par(mfrow=c(2,3), mar=c(1,2,1,2) + 0.1 , oma=c(2.5,4,2.5,2.5))
plot( NULL , xlim=c(min(dc$temp_z),max(dc$temp_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 1:1 ) {
  lines( temp_seq , p_temp_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_PI <- sapply( 1:length(p_temp) , function(i) PI(p_temp[[i]][,k],prob=0.89) )
  shade( p_temp_PI , temp_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = F)
points(agg_temp$temp_z, agg_temp$r_childcare)
mtext(text="Probability (Childcare)",side=2,line=2)

plot( NULL , xlim=c(min(dc$temp_z),max(dc$temp_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 2:2 ) {
  lines( temp_seq , p_temp_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_PI <- sapply( 1:length(p_temp) , function(i) PI(p_temp[[i]][,k],prob=0.89) )
  shade( p_temp_PI , temp_seq,col=col.alpha("black",0.1))
}
points(agg_temp$temp_z, agg_temp$r_food_production)

axis( side = 1, at = labels.at, labels = F)
mtext(text="Probability (Food Production)",side=2,line=2)

plot( NULL , xlim=c(min(dc$temp_z),max(dc$temp_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 3:3 ) {
  lines( temp_seq , p_temp_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_PI <- sapply( 1:length(p_temp) , function(i) PI(p_temp[[i]][,k],prob=0.89) )
  shade( p_temp_PI , temp_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.temp)
points(agg_temp$temp_z, agg_temp$r_household)
mtext(text="Probability (Domestic Work)",side=2,line=2)

plot( NULL , xlim=c(min(dc$temp_z),max(dc$temp_z)) , ylim=c(0,0.4) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 4:4 ) {
  lines( temp_seq , p_temp_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_PI <- sapply( 1:length(p_temp) , function(i) PI(p_temp[[i]][,k],prob=0.89) )
  shade( p_temp_PI , temp_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.temp)
points(agg_temp$temp_z, agg_temp$r_play)
mtext(text="Probability (Play)",side=2,line=2)

plot( NULL , xlim=c(min(dc$temp_z),max(dc$temp_z)) , ylim=c(0,1) , xaxt = "n", main = "",adj=0, ylab = "", xlab = "", cex.main = .9)
for ( k in 5:5) {
  lines( temp_seq , p_temp_mean[k,], lwd = 2 ,col=col.alpha("black",1))
  p_temp_PI <- sapply( 1:length(p_temp) , function(i) PI(p_temp[[i]][,k],prob=0.89) )
  shade( p_temp_PI , temp_seq,col=col.alpha("black",0.1))
}
axis( side = 1, at = labels.at, labels = preferred.temp)
points(agg_temp$temp_z, agg_temp$r_znonworkobs)
mtext(text="Probability (Other Activities)",side=2,line=2)

mtext(text="Annual Mean Temperature (°C)",side=1,line=1,outer=TRUE)

quartz.save(file="figS12.pdf",type="pdf")
rm(agg_temp)
rm(p_temp)
rm(p_temp_mean)
rm(p_temp_PI)
rm(temp)
rm(k)
rm(labels.at)
rm(preferred.temp)
rm(seq.length)
rm(temp_seq)
dev.off()

rm(d)
rm(ds)
