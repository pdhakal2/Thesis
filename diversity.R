#library(phyloseq)
library(vegan)
setwd("../Dropbox/Thesis/")
species.control.matrix  <- data.matrix(read.csv("Aim2/control sample/species/species_level_control_sample.txt",row.names=1,sep="\t",na.strings = "NaN"))
species.control.matrix[is.na(species.control.matrix)] <- 0
control.diversity <- diversity(t(species.control.matrix), index = "shannon")
control.spa <- specaccum(t(species.control.matrix))

species.carcinoma.matrix  <- data.matrix(read.csv("Aim2/carcinoma sample/species/species_level_carcinoma_sample.txt",row.names=1,sep="\t",na.strings = "NaN"))
species.carcinoma.matrix[is.na(species.carcinoma.matrix)] <- 0
carcinoma.diversity <- diversity(t(species.carcinoma.matrix), index = "shannon")
carcinoma.spa <- specaccum(t(species.carcinoma.matrix))

species.adenoma.matrix  <- data.matrix(read.csv("Aim2/advanced adenoma sample/species/species_level_advance_adenoma_sample.txt",row.names=1,sep="\t",na.strings = "NaN"))
species.adenoma.matrix[is.na(species.adenoma.matrix)] <- 0
adenoma.diversity <- diversity(t(species.adenoma.matrix), index = "shannon")
adenoma.spa <- specaccum(t(species.adenoma.matrix))

summary(control.diversity)
summary(adenoma.diversity)
summary(carcinoma.diversity)
boxplot(control.diversity, adenoma.diversity, carcinoma.diversity, names=c("control","adenoma","carcinoma"),main="species richness")

plot(control.spa,col="black",main="Species Accumulation Curve",ylab="Richness",xlab="Samples")
lines(adenoma.spa,col="blue")
lines(carcinoma.spa,col="red")
legend("bottomright",legend=c("control","adenoma","carcinoma"),pch=c(16,16,16),col=c("black","blue","red"))


