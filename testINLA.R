library(ranger)
library(rasterVis)
library(RColorBrewer)
library(INLA)
source("~/Documents/mobileAP/INLA_functions.R")

library(sf)
library(dplyr)
library(terra)
library(APMtools)
library(xgboost)
library(ggplot2)
library(raster)
covnames = c("b0", "nightlight_450", "population_1000", "population_3000",
             "road_class_1_5000", "road_class_2_100", "road_class_3_300",  
             "trop_mean_filt", "road_class_1_100")
formula = as.formula(paste0('y ~ 0 + ', paste0(covnames, collapse = '+'), " + f(s, model = spde)"))

#scale covariates 


#data preparation 
a = read_sf("~/Documents/GitHub/mobileAP/AMS_MixedModels_25May2021.shp")
a = st_transform(a, 4326)
a = a["Mixed_NO2"]%>%filter(Mixed_NO2>1)
samp = st_sample(a, 10000)
samp = samp [!st_is_empty(samp)]
# the most costly, but takes not so long
asp = as_Spatial(a)
 

ams = read_sf("~/Documents/GitHub/AP_AMS/AMS_polygon/")
ams=st_transform(ams, 4326)

s = rast(list.files("~/Documents/GitHub/global-mapping/denl/ams", full.names = T))
s = crop(s, ams)

rs_template=  rast("~/Documents/GitHub/global-mapping/denl/ams/road_class_2_25.tif")

no2 = terra::rasterize(vect(a), rs_template,  field = "Mixed_NO2", fun = mean)
no2 =terra::crop(no2, s[[1]])


add(s)<-no2

r225 = s[["road_class_2_25"]]
s= stack(no2,s)

nlayers(s)
samp = st_cast(samp, "POINT")
samp_sp =  as(samp, "Spatial")

# sample all traffic points
samp_r= terra::extract(s, vect(samp_sp), xy = T)
sf_r225 = st_as_sf(samp_r, coords = c("x","y"),crs = 4326)
quantile(sf_r225$road_class_2_25, c(0.25, 0.5, 0.92, 0.95)) # 9 percent of the data, from 100,000
traffic2 = sf_r225%>%filter(road_class_2_25>0.1&road_class_1_25<0.1)
traffic1 = sf_r225%>%filter(road_class_1_25>0.1&road_class_2_25<0.1)
traffic = rbind(traffic1, traffic2)
nrow(traffic2)
#notraffic 
samp_r= terra::extract(s,vect(samp_sp), xy= T)
sf_r = st_as_sf(samp_r , coords = c("x","y"),crs = 4326)
no_traf = sf_r%>%filter(road_class_2_500<0.1&road_class_1_500<0.1)

ggplot(traffic1)+geom_sf()
 
ggplot(traffic2)+geom_sf()
 
ggplot(traffic)+geom_sf()
 
ggplot(no_traf)+geom_sf()
 
create_dataset = function(traffic = traffic, no_traffic=no_traf,ras_stack = s2, n_traffic= 300, n_no_traffic= 50)
{
  sample_traffic = st_sample(traffic, n_traffic)%>%st_cast("POINT")
  st_crs(sample_traffic) = 4326
  
  sample_no_traffic = st_sample(no_traffic, n_no_traffic)%>%st_cast("POINT")
  st_crs(sample_no_traffic) = 4326
  p = ggplot()+geom_sf(data = sample_traffic, col = "red") +geom_sf(data =sample_no_traffic,color= "green")
  ggsave(paste0("~/Downloads/sample_", i, ".png"), p)
  
  
  tn = c(sample_no_traffic, sample_traffic)
  sptn = as(tn, "Spatial")
  terra::extract(ras_stack, vect(sptn), xy= T)%>%data.frame()
  
}

data0  = create_dataset (traffic = traffic, no_traffic=no_traf,ras_stack = s, n_traffic= s1[i], n_no_traffic= s2[i])

d2 = data0
d2$b0 = 1 # intercept

d2 = d2%>%rename(coox = x, cooy=y, y = Mixed_NO2)
d2$real = d2$y
nrow(d2)
head(d2)
d2= d2%>%dplyr::select(covnames, real, coox, cooy, y)%>%data.frame

lres <- fnFitModelINLA(d2, dp = d2, covnames, formula = formula, TFPOSTERIORSAMPLES = TRUE, family = "gaussian")
lres[[1]]$summary.fixed
plot(lres[[1]]$marginals.fix[[1]], type='l', xlab=expression(Intercept), ylab="Density", cex.lab=1.6, cex.axis=1.4)
# Precision for Gaussian observations
plot(lres[[1]]$marginals.hy[[1]], type='l',xlab=expression(tau[y]), ylab="Density", cex.lab=1.6, cex.axis=1.4)

# Marginal for tau of the spatial random field
plot(lres[[1]]$marginals.hy[[2]], type='l',xlab=expression(tau[x]), ylab="Density", cex.lab=1.6, cex.axis=1.4)

# Marginal for kappa of the spatial random field
plot(lres[[1]]$marginals.hy[[3]], type='l',xlab=expression(kappa), ylab="Density", cex.lab=1, cex.axis=1)


#=========================================
#     Get predicted data on grid
#=========================================

index.pred <- inla.stack.index(lres[[2]], "pred")$data

pred_mean <- lres[[1]]$summary.fitted.values[index.pred, "mean"]
pred_ll <- lres[[1]]$summary.fitted.values[index.pred, "0.025quant"]
pred_ul <- lres[[1]]$summary.fitted.values[index.pred, "0.975quant"]

#==============================
#     Plot of predictions
#==============================

library(rgdal)
library(viridis)
#shapefile <- readOGR(dsn = "~/Documents/GitHub/AP_AMS/AMS_polygon/")
#sjer_aoi_WGS84 <- spTransform(shapefile)
data_pred = data.frame(pred_mean, pred_ll, pred_ul, lres[[5]][, 1], lres[[5]][, 2])
colnames(data_pred) <- c("pred_mean", "pred_ll", "pred_ul", "Longitude", "Latitude")
#==================================
# predictions,  
#=================================
resdf = data.frame(observation =d2$y, prediction_mean = data_pred$pred_mean, 
                   prediction_upper = data_pred$pred_ul,
                   prediction_lower = data_pred$pred_ll, lat =d2$cooy, 
                   lon = d2$coox)

gres = gather(resdf, "key", "value",-lat, -lon)
head(gressf)
gressf = st_as_sf(gres, coords = c("lon", "lat"), crs = 4326)
ggplot() + geom_sf(data = ams)+
  geom_sf(data =gressf,aes(colour=value)) +
  facet_wrap(~key)+
  ggtitle("Observation vs. INLA prediction") +
  theme(plot.title = element_text(hjust = 0, size = 10),
        strip.text.x = element_text(size = 15, colour = "black", angle = 0))+
  scale_color_gradientn(name = expression(paste(NO[2],~mu, g/m^{3})),
                        colours = c(viridis(100, begin = 0.3, end = 0.9),rev( magma(100, begin = 0.3))), limits = c(0,48))

#geom_point(aes(colour=value)) +facet_wrap(~key)+

#geom_polygon(data = shapefile, colour = "black",aes(x = long, y = lat, group = group), fill = NA) + 
#ggsave("~/Documents/GitHub/uncertainty/pred.png",height = 10, width = 10)


##============================================
## differences between prediction and observations, 
##=============================================

res_dif = data.frame(obs_predmean = d2$real- data_pred$pred_mean,  
                     obs_predul = d2$real- data_pred$pred_ul,
                     obs_predll = d2$real- data_pred$pred_ll,
                     lat =d2$cooy, lon = d2$coox)
gres = gather(res_dif, "key", "value",-lat, -lon)
gressf = st_as_sf(gres, coords = c("lon", "lat"), crs = 4326)


ggplot() + geom_sf(data = ams)+
  geom_sf(data =gressf,aes(colour=value)) +
  facet_wrap(~key)+
  ggtitle("Observation - INLA prediction") +
  theme(plot.title = element_text(hjust = 0, size = 10),
        strip.text.x = element_text(size = 15, colour = "black", angle = 0))+
  scale_color_gradientn(name = expression(paste(NO[2],~mu, g/m^{3})),
                        colours = rev(brewer.pal(10,"Spectral")))

ggplot(gres, aes(lon,lat)) + 
  ggtitle("Observation - INLA prediction") +
  theme(plot.title = element_text(hjust = 0, size = 10),
        strip.text.x = element_text(size = 15, colour = "black", angle = 0))+
  geom_point(aes(colour=value)) +facet_wrap(~key)+ 
  geom_polygon(data = shapefile, colour = "black",aes(x = long, y = lat, group = group), fill = NA) + 
  scale_color_gradientn(name = expression(paste(NO[2],~mu, g/m^{3})),
                        colours = rev(brewer.pal(10,"Spectral")))

