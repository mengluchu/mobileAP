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
 
lmcalc = function(data0, i ,s){
  y_var = "Mixed_NO2"
  
  covnames = c( "nightlight_450", "population_1000", "population_3000",
               "road_class_1_5000", "road_class_2_100", "road_class_3_300",  
               "trop_mean_filt", "road_class_1_100")
  formula = as.formula(paste0(y_var,"~", paste0(covnames, collapse = '+')))
  data1 = data0%>%dplyr::select(covnames,y_var)%>%na.omit()  
  sdf = as.data.frame(s, xy = T)
  sdf$b0=1
   
  fit1 = lm(formula, data = data1) 
  pr = predict(fit1, newdata = sdf)
   #beginCluster()
  
  predf = data.frame(pr, x=sdf$x,y= sdf$y)
  preras = terra::rasterize(x = as.matrix(predf[,2:3]), y = s[[1]], values =predf[,1], fun=mean)
  
  terra::writeRaster(preras,paste("~/Downloads/LM", i, ".tif", sep = "_"),
                     overwrite =T)  
  
  #ov <- clusterR(rst, predict, args=list(b, fun = predfun))
 # writeRaster(ov, pred_name, overwrite = TRUE)
  #endCluster()
}



xgbcalc = function(data0, i ,s){
  y_var = "Mixed_NO2"
  data0 = data0%>%dplyr::select(-ID,-x,-y)%>%na.omit()
  
  xgbname = paste0("xgb", i,".png") 
  bst = xgbwrap(s, df_var = data0, y_var = y_var , xgbname = xgbname,
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

inlacalc = function(data0, i , s ){
  d2 = data0
  d2$b0 = 1 # intercept
  d2 = d2%>%rename(coox = x, cooy=y, y = Mixed_NO2)
  d2$real = d2$y
  
  covnames = c("b0", "nightlight_450", "population_1000", "population_3000",
               "road_class_1_5000", "road_class_2_100", "road_class_3_300",  
               "trop_mean_filt", "road_class_1_100")
  formula = as.formula(paste0('y ~ 0 + ', paste0(covnames, collapse = '+'), " + f(s, model = spde)"))
  
  d2= d2%>%dplyr::select(covnames, real, coox, cooy, y)%>%data.frame
  # form validation dataset
  train = st_as_sf(data0,coords= c("x","y"),crs =4326)
  sel = st_distance(train,traffic)
  vali_traf= traffic[which(apply(sel, 2, min)>50),]%>%na.omit()%>%sample_n(100, replace = T)
  # validation points are at least 50m away from training points
  sel = st_distance(train,no_traf)
  vali_notraf= no_traf[which(apply(sel, 2, min)>50),]%>%na.omit()%>%sample_n(100, replace = T)
  
  st_geometry(vali_traf) = NULL
  st_geometry(vali_notraf) = NULL
  
  sdf = as.data.frame(s, xy = T)
  sdf$b0=1
  sdf = sdf %>%rename(coox = x, cooy=y)
  
  lres <- fnFitModelINLA(d2,  sdf, i = i , covnames = covnames, formula = formula, TFPOSTERIORSAMPLES = FALSE, family = "gaussian")
  
  index.pred <- inla.stack.index(lres[[2]], "pred")$data
  
  pred_mean = lres[[1]]$summary.fitted.values[index.pred, "mean"] 
  predf = data.frame(pred_mean, x= sdf$x,y= sdf$y)
  preras = terra::rasterize(x = as.matrix(predf[,2:3]), y = s[[1]], values =predf[,1], fun=mean)
  
  terra::writeRaster(preras,paste("~/Downloads/INLA", i, ".tif", sep = "_"),
                     overwrite =T)  
}

ukcalc = function(data0, i, s)
{
  train = st_as_sf(data0,coords= c("x","y"),crs =4326)%>%na.omit()
  y_var = "Mixed_NO2"
  covnames = c( "nightlight_450", "population_1000", "population_3000",
                "road_class_1_5000", "road_class_2_100", "road_class_3_300",  
                "trop_mean_filt", "road_class_1_100")
  formula = as.formula(paste0(y_var,"~", paste0(covnames, collapse = '+')))
  v_emp_uK1 <- gstat::variogram(formula,train)
  
  v_mod_uK <- automap::autofitVariogram(formula, as(train, "Spatial"))$var_model
  png(paste0("~/Downloads/vgm", i, ".png"))
  plot(v_emp_uK1,v_mod_uK)
  dev.off()
  write.table(v_mod_uK,file = paste0("~/Downloads/vgm", i, ".txt"),sep = "\t")
  sdf = as.data.frame(s, xy = T)
  sdfsp =st_as_sf(sdf, coords = c("x","y"),crs = 4326)
  uK <- krige(
    formula, train, 
    sdfsp,                 # locations to interpolate at
    model = v_mod_uK           # the variogram model fitted above
  ) 
  
  predf = data.frame(uK$var1.pred, x= sdf$x ,y=sdf$y)
  preras = terra::rasterize(x = as.matrix(predf[,2:3]), y = s[[1]], values =predf[,1], fun=mean)
  terra::writeRaster(preras,paste("~/Downloads/UK", i, ".tif", sep = "_"),
                     overwrite =T)  
  
}
}

myTheme2<- rasterTheme(region = c(brewer.pal(11, "RdBu")),strip.background = list(col = 'transparent'),
                       strip.border = list(col = 'transparent'),  axis.line = list(col = "transparent") )

#just wrap things together
createmaps = function(){
  myTheme<- rasterTheme(region = c(brewer.pal(7, "YlGnBu"), "orange","red"),strip.background = list(col = 'transparent'),
                        strip.border = list(col = 'transparent'),  axis.line = list(col = "transparent") )
  myTheme2<- rasterTheme(region = c(brewer.pal(11, "RdBu")),strip.background = list(col = 'transparent'),
                         strip.border = list(col = 'transparent'),  axis.line = list(col = "transparent") )
  
  name_attr = paste0("traf_", s1)
  
  pdf("~/Downloads/result/xgb_5model.pdf")
  levelplot(
    raster::stack(list.files("~/Downloads/", pattern = "xgb*" , full.names = TRUE)),
    par.settings = myTheme, names.attr = name_attr)
  dev.off()
  
  pdf("~/Downloads/result/lm_5model.pdf")
  levelplot(
    raster::stack(list.files("~/Downloads/", pattern = "LM_*" , full.names = TRUE)),
    par.settings = myTheme, names.attr = name_attr)
  dev.off()
  
  
  pdf("~/Downloads/result/lmdif_5model.pdf")
  rs = rast(list.files("~/Downloads/", pattern = "LM_*" , full.names = TRUE))
  dif = rs - no2
  dif = ifel(dif>-20, dif,-20) 
  levelplot(dif,
            par.settings = myTheme2, names.attr = name_attr)
  dev.off()
  
  pdf("~/Downloads/result/xgbdif_5model.pdf")
  rs = rast(list.files("~/Downloads/", pattern = "xgb*" , full.names = TRUE))
  dif = rs - no2
  dif = ifel(dif>-20, dif,-20) 
  levelplot(dif,
            par.settings = myTheme2, names.attr = name_attr)
  dev.off()
  
  pdf("~/Downloads/result/INLA_5model.pdf")
  levelplot(
    raster::stack(list.files("~/Downloads/", pattern = "INLA_*" , full.names = TRUE)),
    par.settings = myTheme, names.attr = name_attr)
  dev.off()
  
  inlaras= raster::stack(list.files("~/Downloads/", pattern = "INLA_*" , full.names = TRUE))
  dif = rast(inlaras) - no2
  dif = ifel(dif>-20, dif,-20) 
  
  pdf("~/Downloads/result/INLA_difmodel.pdf")
  
  levelplot(
    dif, par.settings = myTheme2, names.attr = name_attr)
  dev.off()
  
  stest  =raster::stack(list.files("~/Downloads/", pattern = "INLA*" , full.names = TRUE))
  
  stest2  =raster::stack(list.files("~/Downloads/", pattern = "xgb*" , full.names = TRUE))
  stest3  =raster::stack(list.files("~/Downloads/", pattern = "LM_*" , full.names = TRUE))
  
  stest2 = stest2-stest+stest
  
  xgbdif15= stest2[[5]]-stest2[[1]]
  inladif15 =stest [[5]]-stest[[1]]
  lmdif15= stest3[[5]]-stest3[[1]]
  
  sta = raster::stack(xgbdif15,inladif15,lmdif15)
  pdf("~/Downloads/result/INLA_xgb_lm15.pdf")
  levelplot(
    sta,par.settings = myTheme2, names.attr = c("xgb","inla","lm"))
  dev.off()
  
  
  xgbdif= rast(stest2[[1]])-no2
  inladif =rast(stest [[1]])-no2 
  lmdif =rast(stest3 [[1]])-no2 
  xgbdif = ifel(xgbdif>-20, xgbdif,-20) 
  inladif = ifel(inladif>-20, inladif,-20) 
  lmdif = ifel(lmdif>-20, lmdif,-20) 
  
  xgbdif = ifel(xgbdif<20, xgbdif,20) 
  inladif = ifel(inladif<20, inladif,20) 
  lmdif = ifel(lmdif<20, lmdif,-20) 
  
  add(xgbdif)<-inladif
  add(xgbdif)<-lmdif
  pdf("~/Downloads/result/comparetoreal_1.pdf")
  levelplot(
    xgbdif,par.settings = myTheme2, names.attr = c("xgb","inla","lm"))
  dev.off()
  
  
  pdf("~/Downloads/result/lm_inla1.pdf")
  levelplot(stest3[[1]]-stest[[1]],par.settings = myTheme2)
    dev.off()
    
    pdf("~/Downloads/result/lminla5.pdf")
    levelplot(stest3[[5]]-stest[[5]],par.settings = myTheme2)
    dev.off()
    
    xgbdif= rast(stest2[[5]])-no2
    inladif =rast(stest [[5]])-no2 
    lmdif =rast(stest3 [[5]])-no2 
    xgbdif = ifel(xgbdif>-20, xgbdif,-20) 
    inladif = ifel(inladif>-20, inladif,-20) 
    lmdif = ifel(lmdif>-20, lmdif,-20) 
    
    xgbdif = ifel(xgbdif<20, xgbdif,20) 
    inladif = ifel(inladif<20, inladif,20) 
    lmdif = ifel(lmdif<20, lmdif,-20) 
    
    add(xgbdif)<-inladif
    add(xgbdif)<-lmdif
    pdf("~/Downloads/result/comparetoreal_5.pdf")
    levelplot(
      xgbdif,par.settings = myTheme2, names.attr = c("xgb","inla","lm"))
    dev.off()
    
  
  
  
  xgbdif= rast(stest2[[1]])-no2
  inladif =rast(stest [[1]])-no2 
  dif = xgbdif -inladif 
  dif = ifel(dif>-10, dif,-10) 
  dif = ifel(dif<10, dif,10) 
  
  
  # both overestimate, who is better?
  # it is obvious that INLA model overestimate more at local roads, and even more with more traffic samples
  dif2 = ifel (xgbdif > 0, dif,NA )
  dif3 = ifel (inladif >0, dif2,NA )
  pdf("~/Downloads/result/xgbinla1_dif_overrestimate.pdf")
  levelplot(dif3,par.settings = myTheme2, names.attr = c("overestimate"), margin =F)
  dev.off() 
  
  # both underestimate, who is better?
  dif2 = ifel (xgbdif < 0, dif,NA )
  dif3 = ifel (inladif <0, dif2,NA )
  pdf("~/Downloads/result/xgbinla1_dif_under.pdf")
  levelplot(dif3,par.settings = myTheme2, names.attr = c("underestimate"), margin =F)
  dev.off() 
  
  xgbdif= rast(stest2[[5]])-no2
  inladif =rast(stest [[5]])-no2 
  dif = xgbdif -inladif 
  dif = ifel(dif>-10, dif,-10) 
  dif = ifel(dif<10, dif,10) 
  
  
  # both overestimate, who is better?
  # it is obvious that INLA model overestimate more at local roads, and even more with more traffic samples
  dif2 = ifel (xgbdif > 0, dif,NA )
  dif3 = ifel (inladif >0, dif2,NA )
  pdf("~/Downloads/result/xgbinla5_dif_overrestimate.pdf")
  levelplot(dif3,par.settings = myTheme2, names.attr = c("overestimate"), margin =F)
  dev.off() 
  
  # both underestimate, who is better?
  dif2 = ifel (xgbdif < 0, dif,NA )
  dif3 = ifel (inladif <0, dif2,NA )
  pdf("~/Downloads/result/xgbinla5_dif_under.pdf")
  levelplot(dif3,par.settings = myTheme2, names.attr = c("underestimate"), margin =F)
  dev.off() 
  
  xgbdif= rast(stest2[[5]])-no2
  inladif =rast(stest [[5]])-no2 
  pdf("~/Downloads/result/INLA_xgb_dif.pdf")
  difxgbinla = stest[[5]]-stest2[[5]]
  levelplot(difxgbinla,par.settings = myTheme2)
  dev.off()
  
  
  
  
 lmdif= rast(stest3[[5]])-no2
  inladif =rast(stest [[5]])-no2 
  dif = lmdif -inladif 
  dif = ifel(dif>-10, dif,-10) 
  dif = ifel(dif<10, dif,10) 
  
  
  # both overestimate, who is better?
  # it is obvious that INLA model overestimate more at local roads, and even more with more traffic samples
  dif2 = ifel (lmdif > 0, dif,NA )
  dif3 = ifel (inladif >0, dif2,NA )
  pdf("~/Downloads/result/lminla5_dif_overrestimate.pdf")
  levelplot(dif3,par.settings = myTheme2, names.attr = c("overestimate"), margin =F)
  dev.off() 
  
  # both underestimate, who is better?
  dif2 = ifel (lmdif < 0, dif,NA )
  dif3 = ifel (inladif <0, dif2,NA )
  pdf("~/Downloads/result/lminla5_dif_under.pdf")
  levelplot(dif3,par.settings = myTheme2, names.attr = c("underestimate"), margin =F)
  dev.off() 
  
 #first
  lmdif= rast(stest3[[1]])-no2
  inladif =rast(stest [[1]])-no2 
  dif = lmdif -inladif 
  dif = ifel(dif>-10, dif,-10) 
  dif = ifel(dif<10, dif,10) 
  
  
  # both overestimate, who is better?
  # it is obvious that INLA model overestimate more at local roads, and even more with more traffic samples
  dif2 = ifel (lmdif > 0, dif,NA )
  dif3 = ifel (inladif >0, dif2,NA )
  pdf("~/Downloads/result/lminla1_dif_overrestimate.pdf")
  levelplot(dif3,par.settings = myTheme2, names.attr = c("overestimate"), margin =F)
  dev.off() 
  
  # both underestimate, who is better?
  dif2 = ifel (lmdif < 0, dif,NA )
  dif3 = ifel (inladif <0, dif2,NA )
  pdf("~/Downloads/result/lminla1_dif_under.pdf")
  levelplot(dif3,par.settings = myTheme2, names.attr = c("underestimate"), margin =F)
  dev.off() 
  
  xgbdif= rast(stest2[[5]])-no2
  inladif =rast(stest [[5]])-no2 
  pdf("~/Downloads/result/INLA_xgb_dif.pdf")
  difxgbinla = stest[[5]]-stest2[[5]]
  levelplot(difxgbinla,par.settings = myTheme2)
  dev.off()
}


fnConstructMesh = function(coo, cutoff_ratio=0.0001, i){
  # meshbuilder()
  # offset: size of the inner and outer extensions around the data locations
  (offset1 = 1/10*max(dist(coo)))  #1/8
  (offset2 = 1/2*max(dist(coo)))
  # max.edge: maximum allowed triangle edge lengths in the region and in the extension
  (maxedge1 = quantile(dist(coo), 0.9))/20 #1/30
  (maxedge2 = quantile(dist(coo), 0.9))/3  #1/5
  # cutoff: minimum allowed distance between points used to avoid building many small triangles around clustered locations
  (cutoff = cutoff_ratio*quantile(dist(coo), 0.5)) #WAS max, [5] 500 -> 1501->1800->2231
  mesh = inla.mesh.2d(loc = coo, offset = c(offset1, offset2), cutoff = cutoff, max.edge = c(maxedge1, maxedge2))
  #png("~/Downloads/mesh", i, ".png")
  plot(mesh)
  points(coo, col = "red")
  write.csv(mesh$n,  file = paste0("~/Downloads/meshnodes", i,".csv"))
  
  #dev.off()
  print(mesh$n)
  return(mesh)
}

#calculation starts here
total = 100
#s1 = c(50, 150,300,500,700) # traffic 7000, 
s1 = c(floor(total *ratio*0.2),floor(total *ratio*0.5),floor(total *ratio*1),floor(total *ratio*1.5),floor(total *ratio*2)) # traffic 7000, 
s2 = total - s1
#s2 =rep(200,5)

set.seed(1)
par(mfrow = c(1,1))
for (i in 1:5){
  data0  = create_dataset (traffic = traffic_l, no_traffic=mid_traf,ras_stack = s, n_traffic= s1[i], n_no_traffic= s2[i])
  #inlacalc(data0 = data0,  i = i, s = s )
  xgbcalc( data0 = data0, i = i, s = s)
  
  lmcalc(data0 = data0,  i = i, s = s )
  ukcalc(data0 = data0,  i = i, s = s )
}

testres = function(i,r1)
{
  r2 = aggregate(r1, fact=i+0.1, fun="mean",na.rm = T)
  a = summary(lm(values(r2[[1]])~ values(r2[[2]])))$r.squared
  dif = r2[[1]]- r2[[2]]
  b = sqrt(mean(values((r2[[1]]- r2[[2]])^2), na.rm=T))
  if (a>0.8)
  { print(i)
    print(res(r2))}
  return(list(a,b))
}

plotr2 = function(n ,res)
{
  res = data.frame(unlist(res[1,]), unlist(res[2,]))
  res = cbind(res, ID= 1:(n-1))
  print(summary(res))
  names(res) = c( "r2","RMSE", "ID")
  resg = gather(res,key, value, -ID)
  ggplot(resg, aes(ID, value))+geom_line(aes(color = factor(key)))+geom_point(aes(color = factor(key)))+xlab ( "aggregation factor")
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
  errormat = function(rasters, no2)
{
    dif= rasters-no2
    totalrmse= sqrt(mean(values(dif^2), na.rm=T))
    r2 = summary(lm(values(rasters)~ values(no2)))$r.squared
    
    traf = ifel(s[["road_class_2_100"]]>0.1|s[["road_class_1_100"]]>0.1, rasters, NA)
    ru = ifel(s[["road_class_2_100"]]<0.1&s[["road_class_1_100"]]<0.1, rasters, NA)
    diftraf = traf-no2
    difru = ru-no2
    trafrmse= sqrt(mean(values(diftraf^2), na.rm=T))
    rurmse= sqrt(mean(values(difru^2), na.rm=T))
    r2_tr =summary(lm(values(traf)~ values(no2)))$r.squared
    r2_ru =summary(lm(values(ru)~ values(no2)))$r.squared
    
    c(totalrmse, trafrmse, rurmse, r2, r2_tr, r2_ru)
  }
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

