# Stefan van Dam
# Janna Jilesen

# Exercise 7
# 13-1-2015

# libraries
library(raster)

# loading data
load("data/GewataB1.rda")
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")
load("data/GewataB5.rda")
load("data/GewataB7.rda")
load("data/vcfGewata.rda")

# removing the unrealistic values from the vcf file (which contains tree cover)
vcfGewata[vcfGewata > 100] <- NA

# putting the bands together in a brick
gewata <- brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7, vcfGewata)
gewata_novcf <- brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7)

# plotting the bands in pairs to show correlations
pairs(gewata)

# plotting the bands with the highest correlations seperately
pairs(brick(GewataB7, vcfGewata))
pairs(brick(GewataB3, vcfGewata))

valuetable <- as.data.frame(getValues(gewata))

# creating a linear model, predicting the VCF with the reflections of the different bands
# in the summary a column named 'Pr(>|t|)' is given. Small values mean that the probability that this correlation results from coincidence is really small
# so the bands with the small p-values are important: these are all bands except for band 7
modelLM <- lm(formula = vcf2000Gewata ~ gewataB1 + gewataB2 + gewataB3 + gewataB4 + gewataB5 + gewataB7, data = valuetable)
summary(modelLM)

# making the predicted tree cover raster and plotting it
predLC <- predict(gewata_novcf, model=modelLM, na.rm=TRUE)
predLC[predLC < 0] <- NA
plot(predLC, main="Predicted forest coverage", zlim=c(0,100))
plot(gewata$vcf2000Gewata, main="Actual forest coverage", zlim=c(0,100))

# calculating RMSE by using rasterlayers
sqdif <- (gewata$vcf2000Gewata-predLC)^2
sqdifdf <- as.data.frame(sqdif)
RMSE_rasterlayer<-sqrt(mean(sqdifdf$layer, na.rm=TRUE))

# calculating RMSE by using CellStats
RMSE_cellstats <- sqrt(cellStats(((gewata$vcf2000Gewata-predLC)^2), stat='mean'))

# calculating difference between the actual and predicted tree coverage per class (RMSE)
differenceTC<-(gewata$vcf2000Gewata-predLC)^2

load("data/trainingPoly.rda")
trainingPoly@data$Code <- as.numeric(trainingPoly@data$Class)
ndvi <- overlay(GewataB4, GewataB3, fun=function(x,y){(x-y)/(x+y)})
classes <- rasterize(trainingPoly, ndvi, field='Code', progress='text')

averagediff<-zonal(differenceTC, classes)
RMSEClasses<-sqrt(averagediff[,2])

# this results in an RMSE=8.87 for cropland, an RMSE=4.99 for forest and an RMSE=11.91 for wetland
# this means that the difference between the predicted and actual tree cover is the biggest for wetland.