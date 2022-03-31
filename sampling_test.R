#### Usually, we have more ground stations close to traffic, does it affect modelling and how?
#### sample from known NO2 map 
## sample from close to road and far away from road. 
## 1. investigate the effects of sampling to prediction accuracy
## 2. evaluate the effects of sampling to prediction pattern
## 3. For different methods, machine learning and geostatistics, understand how the sampling affect method selection
## 4, Evaluate different priors in Bayesian Geostatistical methods.  
library(ranger)
library(rasterVis)
library(RColorBrewer)
library(INLA)
source("~/Documents/mobileAP/INLA_functions.R")
library("gstat")
library("automap")
library(tidyr)
library(sf)
library(dplyr)
library(terra)
library(APMtools)
library(xgboost)
library(ggplot2)
library(raster)

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
## much faster than raster::rasterize 

no2 = terra::rasterize(vect(a), rs_template,  field = "Mixed_NO2", fun = mean)
no2 =terra::crop(no2, s[[1]])
plot(no2)

add(s)<-no2
r225 = s[["road_class_2_25"]]

s= raster::stack(no2,s)

nlayers(s)
samp = st_cast(samp, "POINT")
samp_sp =  as(samp, "Spatial")

# sample all traffic points
samp_r= terra::extract(s, vect(samp_sp), xy = T)
sf_r225 = st_as_sf(samp_r, coords = c("x","y"),crs = 4326)
quantile(sf_r225$road_class_2_25, c(0.25, 0.5, 0.92, 0.95)) # 9 percent of the data, from 100,000
traffic2 = sf_r225%>%filter(road_class_2_25>0.1&road_class_1_25<0.1)
traffic1 = sf_r225%>%filter(road_class_1_25>0.1&road_class_2_25<0.1)
traffic3 = sf_r225%>%filter(road_class_1_25>0.1&road_class_2_25>0.1)

traffic = rbind(traffic1, traffic2, traffic3)
# same as traffic = sf_r%>%filter(road_class_2_25>0.1|road_class_1_25>0.1)
nrow(traffic) # 2059
#notraffic 
samp_r= terra::extract(s,vect(samp_sp), xy= T)
sf_r = st_as_sf(samp_r , coords = c("x","y"),crs = 4326)
no_traf = sf_r%>%filter(road_class_2_500<0.1&road_class_1_500<0.1)
mid_traf = sf_r%>%filter(road_class_2_100<0.1&road_class_1_100<0.1)
sum(mid_traf$elevation%in%no_traf$elevation)
midno_traf = mid_traf[!mid_traf$elevation%in%no_traf$elevation,] #differences between mid_traf and no_traf not close to any of road 1 and road 2within 100
nrow(mid_traf) # 7206
traffic_l = sf_r%>%filter(road_class_2_100>0.1|road_class_1_100>0.1) # the rest 
nrow(traffic_l) #2794
# 10000
ratio = 2794/7203

ggplot(traffic1)+geom_sf()
ggsave("~/Downloads/highway.png")
ggplot(traffic3)+geom_sf()
ggsave("~/Downloads/highpr.png")

ggplot(traffic2)+geom_sf()
ggsave("~/Downloads/primary.png")
ggplot(traffic)+geom_sf()
ggsave("~/Downloads/traffuc.png")

ggplot(no_traf)+geom_sf()
ggsave("~/Downloads/notraffic.png")


nrow(traffic)
plot(no_traf)

#calculation starts here
total = 100
#s1 = c(50, 150,300,500,700) # traffic 7000, 
ratio #traffic / rural ratio
s1 = c(floor(total *ratio*0.2),floor(total *ratio*0.5),floor(total *ratio*1),floor(total *ratio*1.5),floor(total *ratio*2)) # traffic 7000, 
s2 = total - s1
#s2 =rep(200,5)
s1

i =3#traf = 38
osrates = c(1.5,2, 2.5)
set.seed(2)
for (osrate in osrates)
{data0  = create_dataset (traffic = traffic_l, no_traffic=mid_traf,ras_stack = s, n_traffic= s1[i], n_no_traffic= s2[i], oversample_traf = osrate)
print(nrow(data0))
#inlacalc(data0 = data0,  i = i, s = s )
xgbcalc( data0 = data0, i = osrate, s = s)
#lmcalc( data0 = data0, i = osrate, s = s)
}
data0sf = st_as_sf(data0, coords = c("x","y"))
data0sf = st_buffer(data0sf, 0.001)
valiras =  mask(no2, vect(data0sf), inverse= TRUE) 
#valiras = no2
x1 = rast(paste0("Downloads/xgb_",osrates[1],"_.tif"))
x2 = rast(paste0("Downloads/xgb_",osrates[2],"_.tif"))
x3 = rast(paste0("Downloads/xgb_",osrates[3],"_.tif"))

#x0 = rast("Downloads/normallm/LM_3_.tif")
x0 = rast("Downloads/normalxgb/xgb_3_.tif")
add(x0) <-c(x1,x2,x3) 
error0 = x0-valiras
names(error0) = c("balanced", paste0("oversample traffic rate:", osrates)) 
#png("Downloads/xgboversample_error_100_3.png")
levelplot(error0 , par.setting =myTheme2)
#dev.off()
d = data.frame(balanced = errormat(x0[[1]],valiras), oversample2 = errormat(x0[[2]],valiras), oversample3 = errormat(x0[[3]],valiras),oversample4 =errormat(x0[[4]],valiras))
d= d[4:6,]
#d$var = c("rmse","rmse_traf","rmse_ru","r2", "r2_traf", "r2_ru")
d$var = c("r2", "r2_traf", "r2_ru")
d =gather(d,key, value, -var)
ggplot(data=d, aes(x=var, y=value, fill =key)) +geom_bar(stat="identity", position=position_dodge())

#ggsave("Downloads/xgbrmser2_700_5.png")
 

 
set.seed(1)
par(mfrow = c(1,1))
for (i in 1:5){
  data0  = create_dataset (traffic = traffic_l, no_traffic=mid_traf,ras_stack = s, n_traffic= s1[i], n_no_traffic= s2[i])
  #inlacalc(data0 = data0,  i = i, s = s )
  xgbcalc( data0 = data0, i = i, s = s)
  
  lmcalc(data0 = data0,  i = i, s = s )
  ukcalc(data0 = data0,  i = i, s = s )
}


i = 5
r1 = rast(stestuk[[i]])
r1 = rast(stestlm[[i]])
r1 = rast(stestxgb[[i]])

r1= r1-no2+no2
add(r1)<-no2

#xgb 19->r2 0.8, 60 r2 0.9
n = 60
  
res = sapply(2:n, testres, r1)

lm_ = plotr2(60, res = res)
xgb_ =plotr2(60, res = res)
uk_ =plotr2(60, res = res)

r2p= grid.arrange(lm_,uk_,xgb_, nrow = 1, ncol = 3)
ggsave("~/Downloads/result100/scale5.pdf", r2p)
?grid.arrange
r2 = aggregate(r1, fact=19, fun="mean",na.rm = T)
r2
20*25
summary(res)
levelplot( r2[[1]]- r2[[2]], par.setting =myTheme2)
plot(r2)
mean
layerCor(r2,fun = "pearson", na.rm= T)
 
#stestinla  =raster::stack(list.files("~/Downloads/", pattern = "INLA*" , full.names = TRUE))
stestxgb  =raster::stack(list.files("~/Downloads/", pattern = "xgb*" , full.names = TRUE))
stestlm  =raster::stack(list.files("~/Downloads/", pattern = "LM_*" , full.names = TRUE))
stestuk =raster::stack(list.files("~/Downloads/", pattern = "UK_*" , full.names = TRUE))
for(i in 1:5)
{  xgbdif= rast(stestxgb[[i]])-no2
 # inladif =rast(stestinla [[i]])-no2 
  lmdif =rast(stestlm [[i]])-no2 
  ukdif =rast(stestuk [[i]])-no2 
  
  xgbdif = ifel(xgbdif>-20, xgbdif,-20) 
  #inladif = ifel(inladif>-20, inladif,-20) 
  lmdif = ifel(lmdif>-20, lmdif,-20) 
  ukdif = ifel(ukdif>-20, ukdif,-20) 
  ukdif = ifel(ukdif<20, ukdif,-20) 
  xgbdif = ifel(xgbdif<20, xgbdif,20) 
  #inladif = ifel(inladif<20, inladif,20) 
  lmdif = ifel(lmdif<20, lmdif,-20) 
  #add(xgbdif)<-inladif
  add(xgbdif)<-lmdif
  add(xgbdif)<-ukdif
  
  pdf(paste0("~/Downloads/result50/comparetoreal_",i,".pdf"))
  print(levelplot(
    xgbdif,par.settings = myTheme2, names.attr = c("xgb","lm", "uk")))
  dev.off()
}
###
 
  
    ss0 = data.frame()
  
for (i in 1:5)
  {
  #  inl = errormat(rast(stestinla[[i]]), no2)
    u=errormat(rast(stestuk[[i]]), no2)
    l = errormat(rast(stestlm[[i]]), no2)
    x = errormat(rast(stestxgb[[i]]), no2)
    ss= rbind(u,l,x)
    ss0 = rbind(ss0, ss)
  }
  names(ss0 ) = c("rmse_all", "rmse_traf", "rmse_far", "r2_all", "r2_traf", "r2_ru")
  #plot(raster(nrows = dim(ss0)[1], ncols =dim(ss0)[2]/2, vals = as.matrix(ss0[, 1:3])))
  ss0
  plotem = function(ss0, em = 1){
    x1 = ss0 [ seq (1,nrow(ss0), by= 3), em]
    x2 = ss0 [ seq (2,nrow(ss0), by= 3), em]
    x3 = ss0 [ seq (3,nrow(ss0), by= 3), em]
    #x4 = ss0 [ seq (4,nrow(ss0), by= 4), em]
    nameem = names(ss0)
    #rmse1 = data.frame(x2,x3,x4, 1:5)
    rmse1 = data.frame(x1,x2,x3, 1:5)
    #names(rmse1) = c("INLA_SPDE", "UK", "LM", "XGB", "ID")
    names(rmse1) = c( "UK", "LM", "XGB", "ID")
    
    rmse1 = gather(rmse1, key, value, -ID)
    a = ggplot(rmse1, aes(ID, value))+geom_line(aes(color = factor(key)))+geom_point(aes(color = factor(key)))+
      xlab ( "sampling scheme")+ ylab (nameem[em]) + scale_fill_discrete(name=nameem[em])
  }
  #ggsave ("~/Downloads/rmse.pdf")
  rmse1 = plotem(ss0 = ss0, em = 1)
  rmsetr= plotem(ss0 = ss0, em = 2)
  rmseru = plotem(ss0 = ss0, em = 3)
  library(gridExtra)
  rmse_ = grid.arrange(rmse1, rmsetr, rmseru)
  ggsave("~/Downloads/result50/a100.png", rmse_)
  
  
  r2  = plotem(ss0 = ss0, em = 4)
  r21 = plotem(ss0 = ss0, em = 5)
  r22 = plotem(ss0 = ss0, em = 6)
  a800 = grid.arrange(r2, r21, r22)
  ggsave("~/Downloads/result50/a100_r2.png", a800)
  
  
  
#create validation dataset, not needed if just mapping
for (i in 1){
data0  = create_dataset (traffic = traffic, no_traffic=no_traf,ras_stack = s, n_traffic= s1[i], n_no_traffic= s2[i])
data_ = data0%>%dplyr::select(-ID,-x,-y)%>%na.omit()

#ra =ranger(Mixed_NO2~., data = data_)
#print(ra$r.squared)
 
train = st_as_sf(data0,coords= c("x","y"),crs =4326)
sel = st_distance(train,traffic)

vali_traf= traffic[which(apply(sel, 2, min)>50),]%>%na.omit()%>%sample_n(100, replace = T)
 # validation points are at least 50m away from training points

sel = st_distance(train,no_traf)
vali_notraf= no_traf[which(apply(sel, 2, min)>50),]%>%na.omit()%>%sample_n(100, replace = T)

st_geometry(vali_traf) = NULL
st_geometry(vali_notraf) = NULL

#pred <- predict(ra, data = vali_traf)
#pred2 <- predict(ra, data = vali_notraf)

#print(error_matrix(vali_traf$Mixed_NO2, pred$predictions))
#print(error_matrix(vali_notraf$Mixed_NO2, pred2$predictions))
  
}

#test INLA 
#lres <- fnFitModelINLA(d2,  d2, covnames = covnames, formula = formula, TFPOSTERIORSAMPLES = FALSE, family = "gaussian")
#for (i in 1:5){
#  data0  = create_dataset (traffic = traffic, no_traffic=no_traf,ras_stack = s, n_traffic= s1[i], n_no_traffic= s2[i])
#  data_ = data0%>%dplyr::select(-ID,-x,-y)%>%na.omit()
#  inlacalc(data0, data_)
#}
  
## INLA and XGB calculation
 

#########
# UK
#########
    



v_emp_uKsamp <- gstat::variogram(formula,train)
plot(v_emp_uKsamp)
fvs = fit.variogram(v_emp_uKsamp, vgm(30, "Ste", 1.3, 10))
1   Nug 23.52582 0.000000   0.0
2   Ste 12.61439 1.304142   0.5
plot(v_emp_uKsamp,fvs)
fvs
 
end_time <- Sys.time()
end_time - start_time
OK$var1.pred

predf = data.frame(OK$var1.pred, x= sdf$x ,y=sdf$y)
preras = terra::rasterize(x = as.matrix(predf[,2:3]), y = s[[1]], values =predf[,1], fun=mean)

terra::writeRaster(preras,paste("~/Downloads/UK", i, ".tif", sep = "_"),
                   overwrite =T)  


aa =rast("~/Downloads/UK_1_.tif")
aa-no2
levelplot(aa-no2)

i
##########
# try restoring variogram 
#############
 
#v_emp_OK <- gstat::variogram(Mixed_NO2~1,train)
#plot(v_emp_OK)

y_var = "Mixed_NO2"
covnames = c( "nightlight_450", "population_1000", "population_3000",
              "road_class_1_5000", "road_class_2_100", "road_class_3_300",  
              "trop_mean_filt", "road_class_1_100")
formula = as.formula(paste0(y_var,"~", paste0(covnames, collapse = '+')))

no2all = as.data.frame(s[["Mixed_NO2"]], xy = T)
no2sp=st_as_sf(no2all, coords = c("x","y"),crs = 4326)

#v_emp_uK <- gstat::variogram(Mixed_NO2~1,no2sp) # takes very long so I saved it. 
plot(v_emp_uK)
fv = fit.variogram(v_emp_uK, vgm(40, "Ste", 2, 30))
plot(v_emp_uK,fv)

#v_mod_uK <- automap::autofitVariogram(Mixed_NO2~1, as(train, "Spatial"))$var_model
#save(v_emp_uK,file = "~/Downloads/variogramsampleall.rdr")
#model    psill    range kappa
#1   Nug 29.38473 0.000000   0.0
#2   Ste 53.29643 1.988308   0.5
# lower sill is similar to using residuals of a eg linear regression model 
total = 800
s1 = c(50, 150,300,500,700) # traffic 
s2 = total - s1
i = 2
nrow(mid_traf)
data0  = create_dataset (traffic = traffic, no_traffic=mid_traf,ras_stack = s, n_traffic= s1[i], n_no_traffic= s2[i])
train = st_as_sf(data0,coords= c("x","y"),crs =4326)%>%na.omit()
nrow(train)

v_emp_uKsamp <- gstat::variogram(Mixed_NO2~1,train)
plot(v_emp_uKsamp)
fvs = fit.variogram(v_emp_uKsamp, vgm(30, "Ste", 1.3, 120))
plot(v_emp_uKsamp,fvs)
fvs
#best model
#800, 150
#model    psill    range kappa
#1   Nug 30.83981 0.000000   0.0
#2   Ste 62.85634 1.974983   0.5

model    psill    range kappa
1   Nug 34.53604 0.000000   0.0
2   Ste 72.20830 1.607947   0.5
4000
1   Nug 21.72725 0.000000   0.0
2   Ste 63.57828 1.586773   0.5


v_emp_uKsamp <- gstat::variogram(formula,train)
plot(v_emp_uKsamp)
fvs = fit.variogram(v_emp_uKsamp, vgm(30, "Ste", 1.3, 10))
1   Nug 23.52582 0.000000   0.0
2   Ste 12.61439 1.304142   0.5
plot(v_emp_uKsamp,fvs)
fvs
sdf = as.data.frame(s, xy = T)
sdfsp =st_as_sf(sdf, coords = c("x","y"),crs = 4326)
OK <- krige(
  formula, train, 
sdfsp,                 # locations to interpolate at
  model = fvs           # the variogram model fitted above
)    

plot(OK)
plot(no2)
#sdf = as.data.frame(no2, xy = T)
predf = data.frame(OK$var1.pred, x= sdf$x ,y=sdf$y)
preras = terra::rasterize(x = as.matrix(predf[,2:3]), y = s[[1]], values =predf[,1], fun=mean)

terra::writeRaster(preras,paste("~/Downloads/reUK", i, ".tif", sep = "_"),
                   overwrite =T)  


aa =rast("~/Downloads/newUK_2_.tif")
aa-no2
levelplot(aa-no2,par.settings = myTheme)
pdf("~/Downloads/OKdif.pdf")
levelplot(ukdif,par.settings = myTheme2)
dev.off()

