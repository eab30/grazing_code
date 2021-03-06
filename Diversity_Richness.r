## Melica, Emily & Abby attempt code

install.packages("googlesheets")
library(googlesheets)

##install metafor
install.packages("metafor")
require(metafor)

setwd("C:\\Users\\alr569\\Desktop")

#Authorize access to google sheets
gs_auth(new_user=TRUE)

##Download sheet from google
my.sheet <- gs_title("Data Extraction")
my.data <- data.frame(gs_read_csv(ss=my.sheet, ws="Data", is.na(TRUE)))
head(my.data)

##Subset to ritchness data
richness<- subset(my.data,Response=="Richness")
head(richness)

#size of dataset
dim(richness)

#Calculate SD from Var
richness$sd_cont
richness$sd_cont<-(sqrt(richness$CntrlVar))
head(richness$sd_cont)

richness$sd_treat
richness$sd_treat<-(sqrt(richness$TrtVar))
head(richness$sd_treat)

#Remove NAs
richness <- richness[!is.na(richness$CntrlVar),]
richness <- richness[!is.na(richness$TrtMean),]
richness <- richness[!is.na(richness$n.Trt),]
nrow(richness)


#Take standared mean difference and SMD variance
rich<-escalc(measure="SMD",m1i=TrtMean, m2i=CntrlMean,sd1i=sd_treat
              ,sd2i=sd_cont,n1i=n.Trt,n2i=n.Cntrl, data=richness,
              var.names=c("SMD","SMD_var"),digits=4)

##Histogram of headge's d
hist(rich$SMD, breaks=15, xlab="hedge's d", col=2, main="Species Richness")
abline(v=0,col=4,lty=3,lwd=5)

##Forest plot
forest(rich$SMD,rich$SMD_var,slab=rich$Taxa,pch=19)

##Fixed effects model of species richness with NA omitted
fixef.model <- rma(SMD, SMD_var, data=rich, method = "FE")

summary(fixef.model)
plot(fixef.model)
forest(fixef.model,slab=rich$Taxa)


##Random effects model of species richness with NA omitted
ranef.model <- rma(SMD, SMD_var, data=rich, method = "HE")

summary(ranef.model)
plot(ranef.model)
forest(ranef.model,slab=rich$Taxa)

