setwd("C:\\Users\\pdhakal\\Dropbox\\Thesis\\")

aim3y.orig<-read.csv(file="C:\\Users\\pdhakal\\Dropbox\\Thesis\\clinicalinfowithunit.csv", ,row.names = 1,header=TRUE, stringsAsFactors = F,na.strings = "n.a.")
aim3y <- aim3y.orig[aim3y.orig$State %in% c("controls","carcinoma"),]

for(i in 1:18) {print(colnames(aim3y)[i]);aim3y[,i] <- as.numeric(aim3y[,i])}

for(i in 19:34) {print(colnames(aim3y)[i]);aim3y[,i] <- as.factor(aim3y[,i])}
aim3y<- aim3y[,-13]

aim3y$State = abs(as.numeric(aim3y$State)-2)

model <- glm(State ~ ., data = aim3y, family="binomial")

summary(model)
write.table(summary.glm(model)$coefficients,file="cancer_vs_control_regression.csv",quote=F)

##library(nnet)
#multi.fit <- multinom(State ~ ., data = aim3y)

#result <- summary(multi.fit)
#result$