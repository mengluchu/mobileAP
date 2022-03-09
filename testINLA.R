covnames = c("b0", "nightlight_450", "population_1000", "population_3000",
             "road_class_1_5000", "road_class_2_100", "road_class_3_300",  
             "trop_mean_filt", "road_class_1_100")
formula = as.formula(paste0('y ~ 0 + ', paste0(covnames, collapse = '+'), " + f(s, model = spde)"))

#scale covariates 
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

