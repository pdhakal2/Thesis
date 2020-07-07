setwd("C:\\Users\\pdhakal\\Dropbox\\Thesis\\")


aim3x<-read.csv(file="C:\\Users\\pdhakal\\Dropbox\\Thesis\\ERRsampleInfo.csv", header=TRUE, stringsAsFactors = F)
aim3y<-read.csv(file="C:\\Users\\pdhakal\\Dropbox\\Thesis\\clinicalinfo.csv", header=TRUE, stringsAsFactors = F)
alias.index <- match(aim3y$Alias, aim3x$Alias)
aim3x <- aim3x[alias.index,]

controlIndex <- grep("controls",aim3x$title)
carcinomaIndex <- grep("carcinoma",aim3x$title)
aaIndex <- grep("advanced adenoma",aim3x$title)

aim3x.control <- aim3x[controlIndex,]
aim3x.carcinoma <- aim3x[carcinomaIndex,]
aim3x.aa<- aim3x[aaIndex,]

# adding sample id in aim3y
ind <- match(aim3y$Alias, aim3x$Alias)
sample.ids <- aim3x[ind,1]
aim3y <- cbind(sample.ids,aim3y)


h1class<-read.csv(file="C:\\Users\\pdhakal\\Dropbox\\Thesis\\Aim2\\carcinoma sample\\class\\class_level_carcinoma_sample.csv", 
                  header=TRUE, stringsAsFactors = F, na.strings = "NaN",row.names = 1)
h2class <-read.csv(file="C:\\Users\\pdhakal\\Dropbox\\Thesis\\Aim2\\advanced adenoma sample\\class\\class_level_advance_adenoma_sample.csv", 
                   header=TRUE, stringsAsFactors = F, na.strings = "NaN",row.names = 1)

h3class <- read.csv(file="C:\\Users\\pdhakal\\Dropbox\\Thesis\\Aim2\\control sample\\class\\class_level_control_sample.csv", 
                    header=TRUE, stringsAsFactors = F, na.strings = "NaN",row.names = 1)

#h <- cbind(h1class,h2class,h3class)


abundance_table_carcinoma <- t(h1class) #for healthy samples
abundance_table_adenoma <- t(h2class) #for advance adenoma sample
abundance_table_healthy <- t(h3class) #for carcinoma sample

for(i in 3:20) {print(colnames(aim3y)[i]);aim3y[,i] <- as.numeric(aim3y[,i])}

for(i in 21:35) {print(colnames(aim3y)[i]);aim3y[,i] <- as.factor(aim3y[,i])}
#summary(aim3y)

#library(GGally)
#ggpairs(aim3y[,3:20])

#removing column 7 (GPT(ALT))and 15 (total meat)
aim3y <- aim3y[,c(-15)]
#ggpairs(aim3y[,3:18])


#####regression carcinoma

print( colnames(abundance_table_carcinoma))

bacteria_list_carcinoma <- colnames(abundance_table_carcinoma)



abundance_table_carcinoma[is.na(abundance_table_carcinoma)] <- 0



ind <- match(aim3x.carcinoma$Alias,aim3y$Alias)
aim3y.carcinoma <- aim3y[ind,]
for(boi in bacteria_list_carcinoma) {
  print(boi)
  col_number = which(colnames(abundance_table_carcinoma) == boi)
  
  
  indd <- match(aim3y.carcinoma$sample.ids, rownames(abundance_table_carcinoma))
  small_abd_table_carcinoma <- abundance_table_carcinoma[indd,col_number]
  
  big.matrix_carcinoma <- cbind("bacteria"= as.numeric(small_abd_table_carcinoma), aim3y.carcinoma[,-1:-2])
  
  
  fit <- lm(bacteria ~ ., data = big.matrix_carcinoma)
  reg.result_carcinoma <- summary(fit)$coefficients
  reg.result_carcinoma <- cbind("factor" = rownames(reg.result_carcinoma),reg.result_carcinoma)
  reg.result.05_carcinoma <- reg.result_carcinoma[which(as.numeric(reg.result_carcinoma[,5]) < 0.05),, drop=F]
  
  boi <- gsub("\\|",".",boi)
  #write.csv(file = paste(boi,".regression_carcinoma.csv",sep=""),reg.result_carcinoma,row.names=F)
  if(nrow(reg.result.05_carcinoma)>0) write.csv(file = paste("Aim3/class_level_regression/carcinoma/",boi,".regression.05_carcinoma.csv",sep=""),reg.result.05_carcinoma,row.names=F)
  
}




################ regression advance adenoma

ind <- match(aim3x.aa$Alias, aim3y$Alias)
aim3y.aa <- aim3y[ind,]
print( colnames(abundance_table_adenoma))
abundance_table_adenoma[is.na(abundance_table_adenoma)] <- 0

bacteria_list_adenoma <- colnames(abundance_table_adenoma)

for(boi in bacteria_list_adenoma) {
  print(boi)
  col_number = which(colnames(abundance_table_adenoma) == boi)
  
  
  indd <- match(aim3y.aa$sample.ids, rownames(abundance_table_adenoma))
  small_abd_table_adenoma <- abundance_table_adenoma[indd,col_number]
  
  big.matrix_adenoma <- cbind("bacteria"= as.numeric(small_abd_table_adenoma), aim3y.aa[,-1:-2])
  
  
  fit <- lm(bacteria ~ ., data = big.matrix_adenoma)
  reg.result_adenoma <- summary(fit)$coefficients
  reg.result_adenoma <- cbind("factor" = rownames(reg.result_adenoma),reg.result_adenoma)
  reg.result.05_adenoma <- reg.result_adenoma[which(as.numeric(reg.result_adenoma[,5]) < 0.05),, drop=F]
  
  boi <- gsub("\\|",".",boi)
 # write.csv(file = paste(boi,".regression_adenoma.csv",sep=""),reg.result_adenoma,row.names=F)
  if(nrow(reg.result.05_adenoma)>0) write.csv(file = paste("Aim3/class_level_regression/adenoma/",boi,".regression.05_adenoma.csv",sep=""),reg.result.05_adenoma,row.names=F)
  
}


###########regression healthy########

ind <- match( aim3x.control$Alias,aim3y$Alias)
aim3y.control <- aim3y[ind,]
print( colnames(abundance_table_healthy))

abundance_table_healthy[is.na(abundance_table_healthy)] <- 0
bacteria_list_healthy <- colnames(abundance_table_healthy)


for(boi in bacteria_list_healthy[1]) {
  print(boi)
  col_number = which(colnames(abundance_table_healthy) == boi)
  
  
  indd <- match(aim3y.control$sample.ids, rownames(abundance_table_healthy))
  small_abd_table_healthy <- abundance_table_healthy[indd,col_number]
  
  big.matrix_healthy <- cbind("bacteria"= as.numeric(small_abd_table_healthy), aim3y.control[,-1:-2])
  
  
  fit <- lm(bacteria ~ ., data = big.matrix_healthy)
  reg.result_healthy <- summary(fit)$coefficients
  reg.result_healthy <- cbind("factor" = rownames(reg.result_healthy),reg.result_healthy)
  reg.result.05_healthy <- reg.result_healthy[which(as.numeric(reg.result_healthy[,5]) < 0.05),, drop=F]
  
  boi <- gsub("\\|",".",boi)
  #write.csv(file = paste(boi,".regression_healthy.csv",sep=""),reg.result_healthy,row.names=F)
  if(nrow(reg.result.05_healthy)>0) write.csv(file = paste("Aim3/class_level_regression/control/",boi,".regression.05_healthy.csv",sep=""),reg.result.05_healthy,row.names=F)
  
}

