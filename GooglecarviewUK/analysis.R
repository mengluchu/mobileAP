a = read.csv("~/Downloads/measurements-3.csv") #first clip of 2019, yelds 10-24 -10-28
head(a)
tail(a)
a1 = st_as_sf(a,coords = c("longitude", "latitude"))
plot(a1['value'])
summary(as.POSIXlt(a1$utc,format="%Y-%m-%dT%H:%M"))
head(a$utc)
library(devtools)
install_github("mengluchu/APMtools") 
library(APMtools)
 
help(package = APMtools)
