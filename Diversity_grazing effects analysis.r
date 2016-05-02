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
plot(Divfixef.model,main="Diversity Fixed Eff")
forest(Divfixef.model,slab=div$Taxa,main="Diversity Fixed Eff")

### Fixed effects models with various moderators:
Divfixef.model1<- rma(SMD~Taxa+PaperID,SMD_var, data=div, method = "FE")
summary(Divfixef.model1)
plot(Divfixef.model1)
forest(Divfixef.model1,slab=div$Taxa,main="Diversity Fixed Eff with Taxa")

## with entire length
Divfixef.model2<- rma(SMD~EntireLength,SMD_var, data=div, method = "FE")
summary(Divfixef.model2)
plot(Divfixef.model2)
forest(Divfixef.model2,slab=div$Taxa,main="Diversity Fixed Eff with Entire Length")

## with study length
Divfixef.model3<- rma(SMD~StudyLength,SMD_var, data=div, method = "FE")
summary(Divfixef.model3)
plot(Divfixef.model3)
forest(Divfixef.model3,slab=div$Taxa,main="Diversity Fixed Eff with Study Length")

## with study year
Divfixef.model4<- rma(SMD~StudyYear,SMD_var, data=div, method = "FE")
summary(Divfixef.model4)
plot(Divfixef.model4)
forest(Divfixef.model4,slab=div$Taxa,main="Diversity Fixed Eff with Study Year")

## with sample design
Divfixef.model5<- rma(SMD~SampleDesign,SMD_var, data=div, method = "FE")
summary(Divfixef.model5)
plot(Divfixef.model5)
forest(Divfixef.model5,slab=div$Taxa,main="Diversity Fixed Eff with Sample Design")
boxplot(div$SMD~div$SampleDesign)


##Random effects model of species diversity with NA omitted, method = ME
Divranef.model <- rma(SMD, SMD_var, data=div, method = "HE")
summary(Divranef.model)
plot(Divranef.model)
forest(Divranef.model,slab=div$Taxa,main="Diversity Random Eff")

Divranef.model <- rma(SMD, SMD_var, data=div, method = "REML")
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
is.factor(div$Taxa)
##Random Effects Models with various moderators
##taxa and paper ID only: drops paper ID, taxa non-significant
Divranef.model1<- rma(SMD~Taxa+PaperID,SMD_var, data=div, method = "REML")
summary(Divranef.model1)
plot(Divranef.model1)
forest(Divranef.model1,slab=div$Taxa,main="Diversity Random Eff with Taxa")

Divranef.model1b<- rma(SMD~Taxa,SMD_var, data=div, method = "REML")
summary(Divranef.model1b)
plot(Divranef.model1b)
forest(Divranef.model1b,slab=div$Taxa,main="Diversity Random Eff with Taxa")

## with entire length
Divranef.model2<- rma(SMD~EntireLength,SMD_var, data=div, method = "REML")

summary(Divranef.model2)
plot(Divranef.model2)
forest(Divranef.model2,slab=div$Taxa,main="Diversity Random Eff with Entire Length")

## with study length
Divranef.model3<- rma(SMD~StudyLength,SMD_var, data=div, method = "REML")
summary(Divranef.model3)
plot(Divranef.model3)
forest(Divranef.model3,slab=div$Taxa,main="Diversity Random Eff with Study Length")

## with study year
Divranef.model4<- rma(SMD~StudyYear,SMD_var, data=div, method = "REML")
summary(Divranef.model4)
plot(Divranef.model4)
forest(Divranef.model4,slab=div$Taxa,main="Diversity Random Eff with Study Year")

## with sample design
Divranef.model5<- rma(SMD~SampleDesign,SMD_var, data=div, method = "REML")
summary(Divranef.model5)
plot(Divranef.model5)
forest(Divranef.model5,slab=div$Taxa,main="Diversity Random Eff with Sample Design")
boxplot(rich$SMD~rich$SampleDesign)
boxplot(rich$SMD~rich$SampleDesign)

## try to specify paper ID as random effect...no success yet
Divranef.model6<- rma(SMD~(1|as.numeric("PaperID")),SMD_var, data=div, method = "REML")
summary(Divranef.model6)
plot(Divranef.model6)
forest(Divranef.model6,slab=div$Taxa,main="Diversity Random Eff of Paper ID")


Divranef.modelfull<- rma(SMD~Taxa,SMD_var, data=div, method = "REML")
summary(Divranef.modelfull)
plot(Divranef.modelfull)
forest(Divranef.modelfull,slab=div$Taxa,main="Diversity Random Eff with Taxa")

Divranef.modelreduced <- rma(SMD, SMD_var, data=div, method = "REML")
summary(Divranef.modelreduced )
plot(Divranef.modelreduced )
forest(Divranef.modelreduced ,slab=div$Taxa,main="Diversity Random Eff")