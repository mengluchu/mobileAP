
create_dataset = function(traffic = traffic, no_traffic=no_traf,ras_stack = s2, n_traffic= 300, n_no_traffic= 50, oversample_traf=1)
{
  sample_traffic = st_sample(traffic, n_traffic)%>%st_cast("POINT")%>% st_sample(floor( n_traffic*oversample_traf), replace = T)%>%st_cast("POINT") # oversample rate, 1 = no oversample
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