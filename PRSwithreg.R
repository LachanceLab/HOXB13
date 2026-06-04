library(data.table)

# pop <- West
# pop <- South
# pop <- All


#fread or read.table
prsfile <- read.table("/home/rohinijanivara/Documents/ChrisHaimanPRS/Jun16HaimanPRS_withdosageimputation.txt",sep=" ",header=TRUE)
PRS <- prsfile$PRS

#phen file has PrCa status, age, PCs, sites
phen <- read.table("/home/rohinijanivara/Documents/ChrisHaimanPRS/finalcovfile",sep=",",header=TRUE)

# IID <- phen$IID
# PrCa <- phen$PrCa
# PRSAge <- phen$PRSAge
# PC1  <- phen$PC1
# PC2  <- phen$PC2
# PC3  <- phen$PC3
# PC4  <- phen$PC4
# PC5  <- phen$PC5
# PC6  <- phen$PC6
# PC7  <- phen$PC7
# PC8  <- phen$PC8
# PC9  <- phen$PC9
# PC10  <- phen$PC10
# Agecat <- phen$Agecat
# imputeAge <- phen$missingyoungcaseisenrollment
# dxAge <- phen$dxAge

# HOGGY  <- phen$HOGGY
# KBTH  <- phen$KBTH
# Military  <- phen$Military
# SUN  <- phen$SUN
# UCH  <- phen$UCH
# UATH  <- phen$UATH


#print(PRSAge[phen$PrCa==1 & phen$HOGGY==1])

#Continuous

PRS.std<-scale(PRS)
cont <- glm(PrCa ~ PRS.std + PRSAge + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=phen, family=binomial)
cont.s <- summary(cont)
print("Association of PRS with overall prostate cancer risk (continuous)")
print(cont.s)
print ("----------------------------------------------------------------")

#Categorical 1

PRS.categories <-c(0, .1, .2, .3, .4, .6, .7, .8, .9,  1)

#Overall

PRS.ref <- "40% - 60%"
PRS.labels <- paste(paste(100*PRS.categories[1:(length(PRS.categories)-1)], 100*PRS.categories[2:(length(PRS.categories))], sep="% - "), "%", sep="")
PRS.labels.s <- PRS.labels[!PRS.labels==PRS.ref]
q.v <- quantile(PRS[phen$PrCa==0], probs=PRS.categories)
q.v[1] <- min(q.v[1], min(PRS)) # if min(PRS) not in controls
q.v[length(q.v)] <- max(q.v[length(q.v)], max(PRS)) # if max(PRS) is not in controls
PRS.c <- cut(PRS, breaks=q.v, right=T, include.lowest=T, labels=PRS.labels)
phen$PRS.c <- relevel(PRS.c, ref = PRS.ref)

overall1 <- glm(PrCa ~ PRS.c + PRSAge + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=phen, family = binomial)
overall1.s <- summary(overall1)
print("Association of PRS with overall prostate cancer risk (Category 1)")
print(overall1.s)
print ("----------------------------------------------------------------")

# #AgeStrat

#Below 55


below55 <- read.table("/home/rohinijanivara/Documents/ChrisHaimanPRS/Jun20/below55.csv",sep=",",header=TRUE)

below55PRS <- below55$PRS

PRS.ref <- "40% - 60%"
PRS.labels <- paste(paste(100*PRS.categories[1:(length(PRS.categories)-1)], 100*PRS.categories[2:(length(PRS.categories))], sep="% - "), "%", sep="")
PRS.labels.s <- PRS.labels[!PRS.labels==PRS.ref]
q.v <- quantile(below55PRS[below55$PrCa==0], probs=PRS.categories)
q.v[1] <- min(q.v[1], min(below55PRS)) # if min(PRS) not in controls
q.v[length(q.v)] <- max(q.v[length(q.v)], max(below55PRS)) # if max(PRS) is not in controls
below55PRS.c <- cut(below55PRS, breaks=q.v, right=T, include.lowest=T, labels=PRS.labels)
below55$PRS.c <- relevel(below55PRS.c, ref = PRS.ref)


agestrat1young <- glm(PrCa ~ PRS.c + PRSAge + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=below55, family=binomial)
agestrat1young.s <- summary(agestrat1young)
print("Association of PRS with overall prostate cancer risk for cases below 55 (Category 1)")
print("PRS decile cutoffs:")
print(q.v)
print("PRS control counts per decile (age below 55)")
print(table(below55$PRS.c[below55$PrCa==0]))
print("PRS case counts per decile (age below 55)")
print(table(below55$PRS.c[below55$PrCa==1]))
print(agestrat1young.s)

PRS.std<-scale(below55PRS)
cont <- glm(PrCa ~ PRS.std + PRSAge + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=below55, family=binomial)
cont.s <- summary(cont)
print("Association of PRS with overall prostate cancer risk for cases below 55 (continuous)")
print(cont.s)
print("----------------------------------------------------------------")


#Above 55
above55 <- read.table("/home/rohinijanivara/Documents/ChrisHaimanPRS/Jun20/above55.csv",sep=",",header=TRUE)
above55PRS <- above55$PRS

PRS.ref <- "40% - 60%"
PRS.labels <- paste(paste(100*PRS.categories[1:(length(PRS.categories)-1)], 100*PRS.categories[2:(length(PRS.categories))], sep="% - "), "%", sep="")
PRS.labels.s <- PRS.labels[!PRS.labels==PRS.ref]
q.v <- quantile(above55PRS[above55$PrCa==0], probs=PRS.categories)
q.v[1] <- min(q.v[1], min(above55PRS)) # if min(PRS) not in controls
q.v[length(q.v)] <- max(q.v[length(q.v)], max(above55PRS)) # if max(PRS) is not in controls
above55PRS.c <- cut(above55PRS, breaks=q.v, right=T, include.lowest=T, labels=PRS.labels)
above55$PRS.c <- relevel(above55PRS.c, ref = PRS.ref)


agestrat1old <- glm(PrCa ~ PRS.c + PRSAge + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=above55, family=binomial)
agestrat1old.s <- summary(agestrat1old)
print("Association of PRS with overall prostate cancer risk for cases above 55 (Category 1)")
print("PRS decile cutoffs:")
print(q.v)
print("PRS control counts per decile (age above 55)")
print(table(above55$PRS.c[above55$PrCa==0]))
print("PRS case counts per decile (age above 55)")
print(table(above55$PRS.c[above55$PrCa==1]))
print(agestrat1old.s)


PRS.std<-scale(above55PRS)
cont <- glm(PrCa ~ PRS.std + PRSAge + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=above55, family=binomial)
cont.s <- summary(cont)
print("Association of PRS with overall prostate cancer risk for cases above 55 (continuous)")
print(cont.s)
print("----------------------------------------------------------------")

#AgeAss

caseswdx <-read.table("/home/rohinijanivara/Documents/ChrisHaimanPRS/Jun20/caseswithdxage.csv",sep=",",header=TRUE)
caseswdxPRS <- caseswdx$PRS

PRS.ref <- "0% - 10%"
PRS.labels <- paste(paste(100*PRS.categories[1:(length(PRS.categories)-1)], 100*PRS.categories[2:(length(PRS.categories))], sep="% - "), "%", sep="")
PRS.labels.s <- PRS.labels[!PRS.labels==PRS.ref]
q.v <- quantile(PRS[phen$PrCa==0], probs=PRS.categories)
q.v[1] <- min(q.v[1], min(PRS)) # if min(PRS) not in controls
q.v[length(q.v)] <- max(q.v[length(q.v)], max(PRS)) # if max(PRS) is not in controls
caseswdxPRS.c <- cut(caseswdxPRS, breaks=q.v, right=T, include.lowest=T, labels=PRS.labels)
caseswdx$PRS.c <- relevel(caseswdxPRS.c, ref = PRS.ref)

ageass1 <- glm(PRSAge ~ PRS.c+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=caseswdx, family=gaussian)
ageass1.s <- summary(ageass1)
print("Association of PRS with age at diagnosis (Category 1)")
print("PRS decile cutoffs:")
print(q.v)
print("PRS control counts per decile (age above 55)")
print(table(caseswdx$PRS.c[caseswdx$PrCa==0]))
print("PRS case counts per decile (age above 55)")
print(table(caseswdx$PRS.c[caseswdx$PrCa==1]))
print(ageass1.s)


# PRS.std<-scale(caseswdxPRS)
# cont <- glm(PrCa ~ PRS.std + PRSAge + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=caseswdx, family=binomial)
# cont.s <- summary(cont)
# print("Association of PRS with age at diagnosis (continuous)")
# print(cont.s)
print ("----------------------------------------------------------------")



# Categorical 2
PRS.categories <-c(0, .1, .2, .3, .4, .6, .7, .8, .9, .99, 1)

#Overall

PRS.ref <- "40% - 60%"
PRS.labels <- paste(paste(100*PRS.categories[1:(length(PRS.categories)-1)], 100*PRS.categories[2:(length(PRS.categories))], sep="% - "), "%", sep="")
PRS.labels.s <- PRS.labels[!PRS.labels==PRS.ref]
q.v <- quantile(PRS[phen$PrCa==0], probs=PRS.categories)
q.v[1] <- min(q.v[1], min(PRS)) # if min(PRS) not in controls
q.v[length(q.v)] <- max(q.v[length(q.v)], max(PRS)) # if max(PRS) is not in controls
PRS.c <- cut(PRS, breaks=q.v, right=T, include.lowest=T, labels=PRS.labels)
phen$PRS.c <- relevel(PRS.c, ref = PRS.ref)

overall1 <- glm(PrCa ~ PRS.c + PRSAge + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=phen, family = binomial)
overall1.s <- summary(overall1)
print("Association of PRS with overall prostate cancer risk (Category 2)")
print(overall1.s)
print ("----------------------------------------------------------------")

# #AgeStrat

#Below 55

# below55 <- read.table("/home/rohinijanivara/Documents/ChrisHaimanPRS/Jun20/below55.csv",sep=",",header=TRUE)
# below55PRS <- below55$PRS

PRS.ref <- "40% - 60%"
PRS.labels <- paste(paste(100*PRS.categories[1:(length(PRS.categories)-1)], 100*PRS.categories[2:(length(PRS.categories))], sep="% - "), "%", sep="")
PRS.labels.s <- PRS.labels[!PRS.labels==PRS.ref]
q.v <- quantile(below55PRS[below55$PrCa==0], probs=PRS.categories)
q.v[1] <- min(q.v[1], min(below55PRS)) # if min(PRS) not in controls
q.v[length(q.v)] <- max(q.v[length(q.v)], max(below55PRS)) # if max(PRS) is not in controls
below55PRS.c <- cut(below55PRS, breaks=q.v, right=T, include.lowest=T, labels=PRS.labels)
below55$PRS.c <- relevel(below55PRS.c, ref = PRS.ref)


agestrat1young <- glm(PrCa ~ PRS.c + PRSAge + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=below55, family=binomial)
agestrat1young.s <- summary(agestrat1young)
print("Association of PRS with overall prostate cancer risk for cases below 55 (Category 2)")
print("PRS decile cutoffs:")
print(q.v)
print("PRS control counts per decile (age below 55)")
print(table(below55$PRS.c[below55$PrCa==0]))
print("PRS case counts per decile (age below 55)")
print(table(below55$PRS.c[below55$PrCa==1]))
print(agestrat1young.s)

PRS.std<-scale(below55PRS)
cont <- glm(PrCa ~ PRS.std + PRSAge + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=below55, family=binomial)
cont.s <- summary(cont)
print("Association of PRS with overall prostate cancer risk for cases below 55 (continuous)")
print(cont.s)
print("----------------------------------------------------------------")


#Above 55
# above55 <- read.table("/home/rohinijanivara/Documents/ChrisHaimanPRS/Jun20/above55.csv",sep=",",header=TRUE)
# above55PRS <- above55$PRS

PRS.ref <- "40% - 60%"
PRS.labels <- paste(paste(100*PRS.categories[1:(length(PRS.categories)-1)], 100*PRS.categories[2:(length(PRS.categories))], sep="% - "), "%", sep="")
PRS.labels.s <- PRS.labels[!PRS.labels==PRS.ref]
q.v <- quantile(above55PRS[above55$PrCa==0], probs=PRS.categories)
q.v[1] <- min(q.v[1], min(above55PRS)) # if min(PRS) not in controls
q.v[length(q.v)] <- max(q.v[length(q.v)], max(above55PRS)) # if max(PRS) is not in controls
above55PRS.c <- cut(above55PRS, breaks=q.v, right=T, include.lowest=T, labels=PRS.labels)
above55$PRS.c <- relevel(above55PRS.c, ref = PRS.ref)


agestrat1old <- glm(PrCa ~ PRS.c + PRSAge + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=above55, family=binomial)
agestrat1old.s <- summary(agestrat1old)
print("Association of PRS with overall prostate cancer risk for cases above 55 (Category 2)")
print("PRS decile cutoffs:")
print(q.v)
print("PRS control counts per decile (age above 55)")
print(table(above55$PRS.c[above55$PrCa==0]))
print("PRS case counts per decile (age above 55)")
print(table(above55$PRS.c[above55$PrCa==1]))
print(agestrat1old.s)


PRS.std<-scale(above55PRS)
cont <- glm(PrCa ~ PRS.std + PRSAge + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=above55, family=binomial)
cont.s <- summary(cont)
print("Association of PRS with overall prostate cancer risk for cases above 55 (continuous)")
print(cont.s)
print("----------------------------------------------------------------")

#AgeAss

# caseswdx <-read.table("/home/rohinijanivara/Documents/ChrisHaimanPRS/Jun20/caseswithdxage.csv",sep=",",header=TRUE)
# caseswdxPRS <- caseswdx$PRS

PRS.ref <- "0% - 10%"
PRS.labels <- paste(paste(100*PRS.categories[1:(length(PRS.categories)-1)], 100*PRS.categories[2:(length(PRS.categories))], sep="% - "), "%", sep="")
PRS.labels.s <- PRS.labels[!PRS.labels==PRS.ref]
q.v <- quantile(PRS[phen$PrCa==0], probs=PRS.categories)
q.v[1] <- min(q.v[1], min(PRS)) # if min(PRS) not in controls
q.v[length(q.v)] <- max(q.v[length(q.v)], max(PRS)) # if max(PRS) is not in controls
caseswdxPRS.c <- cut(caseswdxPRS, breaks=q.v, right=T, include.lowest=T, labels=PRS.labels)
caseswdx$PRS.c <- relevel(caseswdxPRS.c, ref = PRS.ref)

ageass1 <- glm(PRSAge ~ PRS.c+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=caseswdx, family=gaussian)
ageass1.s <- summary(ageass1)
print("Association of PRS with age at diagnosis (Category 2)")
print("PRS decile cutoffs:")
print(q.v)
print("PRS control counts per decile (age above 55)")
print(table(caseswdx$PRS.c[caseswdx$PrCa==0]))
print("PRS case counts per decile (age above 55)")
print(table(caseswdx$PRS.c[caseswdx$PrCa==1]))
print(ageass1.s)


# PRS.std<-scale(caseswdxPRS)
# cont <- glm(PrCa ~ PRS.std + PRSAge + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + HOGGY + KBTH + Military + SUN + UCH + UATH, data=caseswdx, family=binomial)
# cont.s <- summary(cont)
# print("Association of PRS with age at diagnosis (continuous)")
# print(cont.s)
# print ("----------------------------------------------------------------")