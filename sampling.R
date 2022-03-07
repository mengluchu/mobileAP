#### Usually, we have more ground stations close to traffic, does it affect modelling and how?
#### sample from known NO2 map 
## sample from close to road and far away from road. 
## 1. investigate the effects of sampling to prediction accuracy
## 2. evaluate the effects of sampling to prediction pattern
## 3. For different methods, machine learning and geostatistics, understand how the sampling affect method selection
## 4, Evaluate different priors in Bayesian Geostatistical methods.  
library(ranger)
library(sf)
library(dplyr)
library(terra)
library(APMtools)
library(xgboost)
library(ggplot2)

#library(raster)
a = read_sf("~/Documents/GitHub/mobileAP/AMS_MixedModels_25May2021.shp")
a = st_transform(a, 4326)
a = a["Mixed_NO2"]%>%filter(Mixed_NO2>1)
samp = st_sample(a, 10000)
samp = samp [!st_is_empty(samp)]
 
asp = as_Spatial(a)

plot(no2)

ams = read_sf("~/Documents/GitHub/AP_AMS/AMS_polygon/")
ams=st_transform(ams, 4326)

s = rast(list.files("~/Documents/GitHub/global-mapping/denl/ams", full.names = T))
s = crop(s, ams)

rs_template=  rast("~/Documents/GitHub/global-mapping/denl/ams/road_class_2_25.tif")

no2 = terra::rasterize(vect(a), rs_template,  field = "Mixed_NO2", fun = mean)
no2 =terra::crop(no2, s[[1]])


add(s)<-no2
 
r225 = s[["road_class_2_25"]]
## much faster than raster::rasterize 
no2r =as.raster(no2)
s= stack(no2,s)

nlayers(s)
samp = st_cast(samp, "POINT")
samp_sp =  as(samp, "Spatial")

# sample all traffic points
samp_r= extract(r225, vect(samp_sp), xy = T)
sf_r225 = st_as_sf(samp_r, coords = c("x","y"),crs = 4326)
quantile(sf_r225$road_class_2_25, c(0.25, 0.5, 0.92, 0.95)) # 9 percent of the data, from 100,000
traffic = sf_r225%>%filter(road_class_2_25>0.1) 
st_crs(traffic)
#notraffic 
samp_r= extract(s,vect(samp_sp), xy= T)
sf_r = st_as_sf(samp_r , coords = c("x","y"),crs = 4326)
no_traf = sf_r%>%filter(road_class_2_500<0.1&road_class_1_500<0.1)

plot(traffic)
nrow(traffic)
plot(no_traf)

create_dataset = function(traffic = traffic, no_traffic=no_traf,ras_stack = s2, n_traffic= 300, n_no_traffic= 50)
{
  sample_traffic = st_sample(traffic, n_traffic)%>%st_cast("POINT")
  st_crs(sample_traffic) = 4326
  
  sample_no_traffic = st_sample(no_traffic, n_no_traffic)%>%st_cast("POINT")
  plot(sample_traffic)
  plot(sample_no_traffic)
  st_crs(sample_no_traffic) = 4326
  tn = c(sample_no_traffic, sample_traffic)
  sptn = as(tn, "Spatial")
  extract(ras_stack, vect(sptn), xy= T)%>%data.frame()

}

xgbwrap = function (sr, df_var, y_var, xgbname = "xgb.tif", max_depth,
                    eta, gamma, nrounds, subsample, verbose = 0, xgb_lambda,
                    xgb_alpha, grepstring=NULL)
{ 
  if(!is.null(grepstring))
    dfpredict = subset_grep(df_var, grepstring)
  else dfpredict = df_var%>%dplyr::select(-all_of(y_var))
  sr = raster::subset(sr, names(dfpredict))
  re = names(sr)
  pre_mat3 = df_var %>% dplyr::select(re) # make sure df and rasterstack have the same sequence
  stopifnot(all.equal(names(sr), names(pre_mat3)))
  yvar = df_var %>% dplyr::select(y_var) %>% unlist()
  dfmatrix = as.matrix(pre_mat3)
  x = xgboost(data = dfmatrix, label = yvar, max_depth = max_depth,
              subsample = subsample, eta = eta, gamma = gamma, nrounds = nrounds,
              verbose = verbose, lambda = xgb_lambda, alpha = xgb_alpha)
  return(list(x, sr))
  
}


predfun <- function(model, data) {
  v <- predict(model, as.matrix(data))
}



total = 800
s1 = c(50, 150,300,500,700) # traffic 
s2 = total - s1
#s2 =rep(200,5)

set.seed(1)
for (i in 1:5){
data0  = create_dataset (traffic = traffic, no_traffic=no_traf,ras_stack = s, n_traffic= s1[i], n_no_traffic= s2[i])
data_ = data0%>%dplyr::select(-ID,-x,-y)%>%na.omit()
ra =ranger(Mixed_NO2~., data = data_)
print(ra$r.squared)

# form validation dataset
train = st_as_sf(data0,coords= c("x","y"),crs =4326)
sel = st_distance(train,traffic)

vali_traf= traffic[which(apply(sel, 2, min)>50),]%>%na.omit()%>%sample_n(100, replace = T)
 # validation points are at least 50m away from training points

sel = st_distance(train,no_traf)
vali_notraf= no_traf[which(apply(sel, 2, min)>50),]%>%na.omit()%>%sample_n(100, replace = T)

st_geometry(vali_traf) = NULL
st_geometry(vali_notraf) = NULL

pred <- predict(ra, data = vali_traf)
pred2 <- predict(ra, data = vali_notraf)

print(error_matrix(vali_traf$Mixed_NO2, pred$predictions))
print(error_matrix(vali_notraf$Mixed_NO2, pred2$predictions))

nrow(vali_traf)
nrow(vali_notraf)


#xgb



y_var = "Mixed_NO2"
inde_var = data_ 
xgbname = paste0("xgb", i,".png") 
bst = xgbwrap(s, df_var = inde_var, y_var = y_var , xgbname = xgbname,
              nrounds = 1000, eta = 0.007, gamma =5,max_depth = 6, xgb_alpha = 0, xgb_lambda = 2,
              subsample=0.7,grepstring = NULL)

b = bst[[1]]
rst = bst[[2]]

pred_name = paste("~/Downloads/xgb", i, ".tif", sep = "_")

ov = predict(rst, b, fun = predfun)

#beginCluster()
#ov <- clusterR(rst, predict, args=list(b, fun = predfun))
writeRaster(ov, pred_name, overwrite = TRUE)
#endCluster()

}
library(rasterVis)
library(RColorBrewer)

 myTheme<- rasterTheme(region = c(brewer.pal(7, "YlGnBu"), "orange","red"),strip.background = list(col = 'transparent'),
                      strip.border = list(col = 'transparent'),  axis.line = list(col = "transparent") )
 myTheme2<- rasterTheme(region = c(brewer.pal(11, "RdBu")),strip.background = list(col = 'transparent'),
                       strip.border = list(col = 'transparent'),  axis.line = list(col = "transparent") )

  name_attr = paste0("traf_", s1)
  
pdf("~/Downloads/result/xgb_fixedtotal.pdf")
levelplot(
  raster::stack(list.files("~/Downloads/", pattern = "xgb*" , full.names = TRUE)),
  par.settings = myTheme, names.attr = name_attr)
dev.off()


pdf("~/Downloads/result/xgbdif_fixedtotal.pdf")
rs = rast(list.files("~/Downloads/", pattern = "xgb*" , full.names = TRUE))
dif = rs - no2
dif = ifel(dif>-20, dif,-20) 
 
levelplot(dif,
  par.settings = myTheme2, names.attr = name_attr)
dev.off()
#predict(sr, bst, fun = predfun)

#################
# Inla
#####3

library(INLA)
source("~/Documents/mobileAP/INLA_functions.R")

covnames = c("b0", "nightlight_450", "population_1000", "population_3000",
             "road_class_1_5000", "road_class_2_100", "road_class_3_300",  
             "trop_mean_filt", "road_class_1_100")

#scale covariates 
d2 = data0
d2$b0 = 1 # intercept

d2 = d2%>%rename(coox = x, cooy=y, y = Mixed_NO2)
d2$real = d2$y
nrow(d2)
head(d2)
d2= d2%>%dplyr::select(covnames, real, coox, cooy, y)%>%data.frame

# Variables for stacked generalization
# d$lasso = d$lasso10f_pre
# d$rf = d$rf10f_pre
# d$xgb = d$xgb10f_pre
#d$Countrycode  = as.factor(d$Countrycode)

#spatial and non-spatial model
formula = as.formula(paste0('y ~ 0 + ', paste0(covnames, collapse = '+'), " + f(s, model = spde)"))
#be careful how the mesh is built, it canbe very slow. The mesh should have relatively large maximum age but small cutoff ratio.
#dp:test
lres <- fnFitModelINLA(d2, dp = d2, covnames, formula = formula, TFPOSTERIORSAMPLES = TRUE, family = "gaussian")
lres[[1]]$summary.fixed
nrow(d2)
#d3 = d2 [1:500, ]
#formula2 = as.formula(paste0('y ~ 0 + ', paste0(covnames, collapse = '+')))
#lres2 <- fnFitModelINLA(d3, dp = d3, covnames, formula = formula2, TFPOSTERIORSAMPLES = TRUE, family = "gaussian")
#lres2[[1]]$summary.fixed

# model comparison 
#lres[[1]]$summary.fixed/lres2[[1]]$summary.fixed
#(lres[[1]]$summary.fixed- lres2[[1]]$summary.fixed)/(lres[[1]]$summary.fixed+ lres2[[1]]$summary.fixed)

#===============================================
# model hyperparameter, figure 5 supplementary
#==================================================
#par(mfrow=c(1,1))
# Intercept
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

#
sdf = as.data.frame(s,xy = T)
sdf$b0=1
head(sdf)
sdf = sdf %>%rename(coox = x, cooy=y)

lres <- fnFitModelINLA(d2,  sdf, covnames = covnames, formula = formula, TFPOSTERIORSAMPLES = FALSE, family = "gaussian")

index.pred <- inla.stack.index(lres[[2]], "pred")$data

pred_mean = lres[[1]]$summary.fitted.values[index.pred, "mean"] 
predf = data.frame(pred_mean, x= lres[[5]][, 1],y= lres[[5]][, 2])
names(predf)
rast(predf, type = "xyz") 

preras = terra::rasterize(x = as.matrix(predf[,2:3]), y = s[[1]], values =predf[,1], fun=mean)
plot(preras)
pdf ("~/Downloads/INLApred.pdf")
levelplot(
 preras,  par.settings = myTheme, names.attr = "INLA")              
dev.off()



pdf("~/Downloads/result/INLAdif_fixedtotal.pdf")

dif = preras - no2
dif = ifel(dif>-20, dif,-20) 

levelplot(dif,
          par.settings = myTheme2, names.attr = "INLA - real")
dev.off()
#predict(sr, bst, fun = predfun)
