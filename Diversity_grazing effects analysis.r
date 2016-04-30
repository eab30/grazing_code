## Melica, Emily & Abby attempt code
?rma
install.packages("googlesheets")
library(googlesheets)

##install metafor
install.packages("metafor")
require(metafor)

setwd("Users\\Strongbuck\\Desktop")
setwd("/Users/Strongbuck/Desktop/in progress/MetaAnalysis")
#Authorize access to google sheets
gs_auth(new_user=TRUE)

##Download sheet from google
my.sheet <- gs_title("Data Extraction")
my.data <- data.frame(gs_read_csv(ss=my.sheet, ws="Data", is.na(TRUE)))
head(my.data)

##Subset to diversity data
diversity<-rbind(my.data[which(my.data$Response=="Shannon"),],
                 my.data[which(my.data$Response=="Diversity"),],
                 my.data[which(my.data$Response=="Shannon diversity"),])

unique(diversity$Response)

#size of dataset
dim(diversity)

#Calculate SD from Var
diversity$sd_cont
diversity$sd_cont<-(sqrt(diversity$CntrlVar))
head(diversity$sd_cont)

diversity$sd_treat
diversity$sd_treat<-(sqrt(diversity$TrtVar))
head(diversity$sd_treat)

#Remove NAs
diversity <- diversity[!is.na(diversity$CntrlVar),]
diversity <- diversity[!is.na(diversity$TrtMean),]
diversity <- diversity[!is.na(diversity$nTrt),]
nrow(diversity)


#Take standared mean difference and SMD variance
div<-escalc(measure="SMD",m1i=TrtMean, m2i=CntrlMean,sd1i=sd_treat
              ,sd2i=sd_cont,n1i=nTrt,n2i=nCntrl, data=diversity,
              var.names=c("SMD","SMD_var"),digits=4)

##Histogram of headge's d
hist(div$SMD, breaks=15, xlab="hedge's d", col=2, main="Species diversity",
     xlim=c(-1.0,0.1))
abline(v=0,col=4,lty=3,lwd=5)

##Forest plot
forest(div$SMD,div$SMD_var,slab=div$Taxa,pch=19)

##Fixed effects model of species diversity with NA omitted
Divfixef.model <- rma(SMD, SMD_var, data=div, method = "FE")

summary(Divfixef.model)
plot(Divfixef.model)
forest(Divfixef.model,slab=div$Taxa,main="Diversity Fixed Eff")


##Random effects model of species diversity with NA omitted
Divranef.model <- rma(SMD, SMD_var, data=div, method = "HE")

summary(Divranef.model)
plot(Divranef.model)
forest(Divranef.model,slab=div$Taxa,main="Diversity Random Eff")

##Fail-Safe-N and Trim and Fill Funnel Plot
fsn(div$SMD, div$SMD_var, div$TrtSE, data=div, type="Rosenthal", 
    alpha=.05, digits=4)
fsn(div$SMD, div$SMD_var, div$TrtSE, data=div, type="Rosenberg", 
    alpha=.05, digits=4)

funnel(Divfixef.model,main="Diversity Fixed Effects")

Divtaf <- trimfill(Divfixef.model)
funnel(Divtaf,main="Diversity Trim & Fill Fixed Eff")

## Examining Effect of Grazing Regime and Agency Owner
#not enough levels available at this time to examine this
