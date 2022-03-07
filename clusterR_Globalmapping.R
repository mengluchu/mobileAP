resolution = 25
predictor_dir = "~/Documents/GitHub/global-mapping/"
#for denl or a folder with multiple cities
#citynames = list.files(predictor_dir)
#city = citynames[1]
#city
#city = "daresalam" #though it's a country
city = "london"
lf = list.files(paste0(predictor_dir, city, "/"), full.names = T)
sr = stack(lf)
#install.packages("snow")

ipak <- function(pkg){

  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos='http://cran.muenster.r-project.org')
  sapply(pkg, require, character.only = TRUE)
}
packages <- c( "raster", "dplyr", "devtools", "rgdal","Matrix","xgboost", "data.table" , "randomForest", "glmnet" ,"rasterVis", "sf", "tmap"   )
ipak(packages)
install_github("mengluchu/APMtools")
library(APMtools)

# prepare data for mapping

# prepare rasters for mapping

#use this to merge roads if needed: sr[[names(sr)[grepl("road_class",names(sr))]]]


xgbwrap = function (sr, df_var, y_var, xgbname = "xgb.tif", max_depth,
          eta, gamma, nrounds, subsample, verbose = 0, xgb_lambda,
          xgb_alpha, grepstring)
{
  dfpredict = subset_grep(df_var, grepstring)
  sr = raster::subset(sr, names(dfpredict))
  re = names(sr)
  pre_mat3 = df_var %>% dplyr::select(re)
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

#predict(sr, bst, fun = predfun)

glo = read.csv("~/Documents/GitHub/Global mapping/predictorF2021_data/glo4var_2021.csv")
timeall= c("wkd_day_value", "wkd_night_value","wnd_day_value", "wnd_night_value", "mean_value")
prestring = "road|nightlight|population|temp|wind|trop|indu|elev|radi"

beginCluster()

#3 mins per map, 7 cores
for (i in 1:4)
{ print(i)

  y_var = timeall[i]
  varstring = paste(prestring,y_var,sep="|")
  inde_var = glo %>%dplyr::select(matches(varstring))


  if (resolution == 100){
    inde_var = inde_var%>%select(-c(industry_25,industry_50,road_class_1_25,road_class_1_50,road_class_2_25,road_class_2_50, road_class_3_25, road_class_3_50))}

  bst = xgbwrap(sr, df_var = inde_var, y_var = y_var ,xgbname = xgbname,
                nrounds = 1000, eta = 0.007, gamma =5,max_depth = 6, xgb_alpha = 0, xgb_lambda = 2,
                subsample=0.7,grepstring = prestring)

  b = bst[[1]]
  rst = bst[[2]]

  pred_name = paste("~/Downloads/xgb",city, y_var, ".tif", sep = "_")
  ov <- clusterR(rst, predict, args=list(b, fun = predfun))
  writeRaster(ov, pred_name, overwrite = TRUE)
  }
#ov <- clusterR(sr, calc, args=list(fun = b))


endCluster()

# show on openStreetMap
library(tmap)
webshot::install_phantomjs()
library(mapview)
for (i in 1:4)
{
  y_var = timeall[i]
  pred_name = paste("~/Downloads/adis_prediction/xgb",city, y_var, ".tif", sep = "_")
  ov = raster(pred_name)

locations_sf = st_as_sf(glo, coords = c("Longitude","Latitude"))

# kind of dirty way for visualisation, so that the data has the same legend.

locations_sf$forvis = ifelse(locations_sf[[y_var]]>maxValue(ov), maxValue(ov), locations_sf[[y_var]])
locations_sf$forvis = ifelse(locations_sf$forvis<minValue(ov), minValue(ov), locations_sf$forvis)

#+tm_shape(lnd)+tm_lines()

tmap_options(basemaps = "OpenStreetMap")

osm_valuemean = tm_shape(ov)+
  tm_raster(names(ov), palette = "BrBG", n = 9, alpha = 0.8, title = "Predictions")+
  tm_shape(locations_sf) +
  tm_dots( "forvis",  size = 0.05, title = "NO2 value",
           popup.vars = y_var, palette = "BrBG", n = 9,alpha = 0.8)+
  tm_layout(legend.outside = TRUE)

#+tm_view(basemaps = c('OpenStreetMap'))

tmap_save(osm_valuemean, paste0("~/Downloads/xgb",city, y_var,".html"))

lff <- tmap_leaflet(osm_valuemean)
mapshot(lff,file = paste0("~/Downloads/xgb",city, y_var,".pdf"))
}

library(RColorBrewer)
myTheme<- rasterTheme(region = c(brewer.pal(7, "YlGnBu"), "orange","red"),strip.background = list(col = 'transparent'),
                        strip.border = list(col = 'transparent'),  axis.line = list(col = "transparent") )

pdf("~/Downloads/Daressalam.pdf")
levelplot(
#stack(list.files("~/Downloads/AMS_prediction/", full.names = TRUE)),

stack(list.files("~/Downloads/daresalam_prediction/", full.names = TRUE)),
#
par.settings = myTheme, names.attr = c("wkd_day","wkd_night","wnd_day","wnd_night"))
dev.off()


#validate with AMS data
library(dplyr)
no2ams =read.csv("~/Documents/GitHub/AP_AMS/NO2_4WEKEN.csv", sep = ";")
no2ams%>%summary()
#&Lopend_gemiddelde<45
no2ams = no2ams%>%filter(Jaar==2017)%>%filter(Lopend_gemiddelde>0.2)%>%select(-WKT_LNG_LAT, -WKT_LAT_LNG, -LAT, -LNG)
loc = read.csv("~/Documents/GitHub/AP_AMS/NO2_LOCATIES.csv", sep = ";")
no2ams= merge(no2ams,loc,by.x = "CodeW", by.y = "Code") # get coordinates
no2ams=no2ams%>%st_as_sf(coords = c("LNG", "LAT"), crs = 4326)
plot(no2ams)

srams = stack(list.files("~/Downloads/AMS_prediction/", full.names = TRUE))
ams = extract(srams,no2ams,  sp = T)

ave1 = (ams$xgb_ams_wkd_day_value_* 15+ ams$xgb_ams_wkd_night_value_*9)/24

ave2 = (ams$xgb_ams_wnd_day_value_* 15+ ams$xgb_ams_wnd_night_value_*9)/24
ave = (ave1 *5  + ave2*2 )/7

summary(lm(ams$xgb_ams_wkd_day_value_~ams$Lopend_gemiddelde)) #night 0.36, day = 0.27

summary(lm(ave~ams$Lopend_gemiddelde)) ##0.36
summary(lm(ams$xgb_ams_mean_value_~ams$Lopend_gemiddelde)) #0.35

library(ggplot2)
ggplot()+geom_point(aes(y= ave, x = ams$Lopend_gemiddelde))+geom_abline()+xlab("AMS station measurements")+ ylab("predicted NO2")
ggsave("~/Downloads/AMSvalidate.png")
abline(a=0, b =1, col = "red")

# no big difference compared to using mean_value to fit the model
summary(lm(ams$mean~ams$Lopend_gemiddelde)) #0.35

ggplot()+geom_point(aes(y= ams$xgb_ams_mean_value_, x = ams$Lopend_gemiddelde))+geom_abline()+xlab("AMS station measurements")+ ylab("predicted NO2 using mean")
ggsave("~/Downloads/AMSvalidate.png")
abline(a=0, b =1, col = "red")
