library(rgdal)
library(raster)
library(usdm)
library(devtools)
library(glmmTMB)
library(boot)
library(MuMIn)
library(caret)
library(DHARMa)
library(sp)
library(SearchTrees)
library(randomForest)
library(rpart)
library(parallel)
library(lme4)
library(visreg)
library(ggeffects)
library(ggplot2)
library(sjPlot)
library(QCA)

rm(list=ls())
#load("dat.RData")

#setwd("/Users/maina/Documents/CCVA/Mangrove/Results_FINAL/Dataset_Final/Unnormalised/")
setwd('/Volumes/Data/Projects/CCVA/Mangrove/Results_FINAL/Dataset_Final/Unnormalised')

all.data<-read.csv("all.data.csv")
data.std<-read.csv("data.std.csv")
load("dat.RData")

#ndvi<-raster("ndvi.nc")
#sample.pts<-sampleRandom(ndvi, size=30000, exp=10, na.rm=TRUE, xy=TRUE, ext=NULL, sp=TRUE)  
#dim(dat)
#ndvi.dat<-extract(ndvi,sample.pts )
# ndvi.tree <- createTree(coordinates(all.data[,5:6]))
# inds.ndvi <- knnLookup(ndvi.tree, newdat=coordinates(sample.pts[,2:3]), k=1)
# data.pts<-all.data[inds.ndvi,]
# PredictVar<-all.data[,c(5,6,8:17,25,27)]
#data.std<-apply(X = PredictVar, MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))})
# 
# respVar<-as.data.frame(all.data[,25])
# scaling <- function(y){ 
#   a<-min(y)
#   b<-max(y)
#   N<-length(y)
#   y_a<-y-a
#   b_a<-(b-a)
#   yscaled<-((y_a*(N-1))/(b_a*N))+(1/(2*N))
#   return(yscaled) 
# }
# 
# resp.std<-data.frame(apply(X = respVar, MARGIN = 2,
#                            FUN = function(x) scaling(x)) )
# 
# data.std.fin<-cbind(all.data[,c(7,26)],resp.std,all.data[,19:20],data.std)
# 
# colnames(data.std.fin)[3]<-"ndvi.std"
# 
# source("functions_analyses_interaction_glmmADMB.R")
#
# ##all data is in 'all.data.csv'
# #standardized data is in 'data.std.csv'
# 
# devtools::install_github("glmmTMB/glmmTMB/glmmTMB")
# 
# #setwd('~/OneDrive - Macquarie University/Projects/comred/CCVA/Results_FINAL/Dataset_Final/Unnormalised/')
# 
# 
# setwd("/Users/maina/Documents/CCVA/Mangrove/Results_FINAL/Dataset_Final/Unnormalised/")
# list.files()
# 
#       files<-list.files()
#       files
#      
#       elevation<-raster(files[6])
#       erosion<-raster(files[7])
#       gravity<-raster(files[8])
#       landcover<-raster(files[17])
#       ldi<-raster(files[12])
#       sla<-raster(files[17])
#       slope<-raster(files[18])
#       tidecm<-raster(files[20])
#       ccd<-raster(files[18])
        cdd85<-raster('cddrcp85.nc')
#       tx90p<-raster(files[19])
#       detided<-raster(files[8])
#       
#       
#       slr<-raster("MSL_Map_MERGED_Global_AVISO_NoGIA_Adjust.nc")
#       crop.extent<-extent(27.44199,51,-33.28764,-0.22713)
#       slr.wio<-crop(slr,crop.extent)

##future slr and aviso mean msl
a.26<-stack('/Volumes/Data/Projects/CCVA/Mangrove/Results_FINAL/Dataset_Final/otherdata/SLR_Data/msl_2050/sl.wio.rcp26.nc')
a.45<-stack('/Volumes/Data/Projects/CCVA/Mangrove/Results_FINAL/Dataset_Final/otherdata/SLR_Data/msl_2050/sl.wio.rcp45.nc')
a.60<-stack('/Volumes/Data/Projects/CCVA/Mangrove/Results_FINAL/Dataset_Final/otherdata/SLR_Data/msl_2050/sl.wio.rcp60.nc')
a.85<-stack('/Volumes/Data/Projects/CCVA/Mangrove/Results_FINAL/Dataset_Final/otherdata/SLR_Data/msl_2050/sl.wio.rcp85.nc')

aviso.msl<-raster('/Volumes/Data/Projects/CCVA/Mangrove/Results_FINAL/Dataset_Final/otherdata/SLR_Data/altimeter_msl.nc')
pr1<-crs(aviso.msl)

a.26pr<-projectRaster(a.26, crs=pr1, method='bilinear')
a.45pr<-projectRaster(a.45, crs=pr1, method='bilinear')
a.60pr<-projectRaster(a.60, crs=pr1, method='bilinear')
a.85pr<-projectRaster(a.85, crs=pr1, method='bilinear')

futmsl<-stack(a.26pr,a.45pr,a.60pr,a.85pr)

wcs.msl<-calc(futmsl, max)
bcs.msl<-calc(futmsl, min)

#futmsl.msl<-stack(futmsl,bcs.msl,wcs.msl)
futmsl.msl<-stack(bcs.msl,wcs.msl)

#names(futmsl.msl)<-c("msl26.2050","msl45.2050","msl60.2050","msl85.2050",'bcs.msl','wcs.msl')
names(futmsl.msl)<-c('bcs.msl','wcs.msl')



fmsl.df<-as.data.frame(futmsl.msl, na.rm=T, xy=T)
fmsl.tree <- createTree(coordinates(fmsl.df[,1:2]))
inds.fmsl <- knnLookup(fmsl.tree, newdat=coordinates(all.data[,2:3]), k=1)
fmsl.pts<-fmsl.df[inds.fmsl,]


dim(fmsl.pts)


##trend
msl.26<-stack('/Volumes/Data/Projects/CCVA/Mangrove/Results_FINAL/Dataset_Final/otherdata/SLR_Data/msl.timeseries/ssh_icdc/Sea_Level_Rise_EnsembleMean_20to50_RCP26.nc')
msl.45<-stack('/Volumes/Data/Projects/CCVA/Mangrove/Results_FINAL/Dataset_Final/otherdata/SLR_Data/msl.timeseries/ssh_icdc/Sea_Level_Rise_EnsembleMean_20to50_RCP45.nc')
#msl.60<-stack('/Volumes/Data/Projects/CCVA/Mangrove/Results_FINAL/Dataset_Final/otherdata/SLR_Data/msl.timeseries/rcp60.2020.2050.nc')
msl.85<-stack('/Volumes/Data/Projects/CCVA/Mangrove/Results_FINAL/Dataset_Final/otherdata/SLR_Data/msl.timeseries/ssh_icdc/Sea_Level_Rise_EnsembleMean_20to50_RCP85.nc')
#indices<-rep(1:31,each=12)
#msl.26.annual<-stackApply(msl.26, indices, fun = mean)
#msl.45.annual<-stackApply(msl.45, indices, fun = mean)
#msl.60.annual<-stackApply(msl.60, indices, fun = mean)
#msl.85.annual<-stackApply(msl.85, indices, fun = mean)                         

#time <- 1:nlayers(msl.26.annual)
time <- 1:nlayers(msl.85)

t.trend=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[2] }}

msl.26.trend=calc(msl.26, t.trend)
msl.45.trend=calc(msl.45, t.trend)
#msl.60.trend=calc(msl.60.annual, t.trend)
msl.85.trend=calc(msl.85, t.trend)



#futmsl.trend<-stack(msl.26.trend,msl.45.trend,msl.60.trend,msl.85.trend)
futmsl.trend<-stack(msl.26.trend,msl.45.trend,msl.85.trend)

wcs.trend<-calc(futmsl.trend, max)
bcs.trend<-calc(futmsl.trend, min)

futmsl.trend<-stack(futmsl.trend,wcs.trend,bcs.trend)
names(futmsl.trend)<-c('msl.26.trend','msl.45.trend','msl.85.trend','wcs.trend','bcs.trend')


fmsltrend.df<-as.data.frame(futmsl.trend, na.rm=T, xy=T)
fmsltrend.tree <- createTree(coordinates(fmsltrend.df[,1:2]))
inds.fmsltrend <- knnLookup(fmsltrend.tree, newdat=coordinates(all.data[,2:3]), k=1)
fmsltrend.pts<-fmsltrend.df[inds.fmsltrend,]

#convert to mm
fmsltrend.pts<-fmsltrend.pts[,3:7]*1000
headz)
colnames(all.data)
#all.data1<-cbind(all.data[,c(2:25,30)],fmsl.pts[,3:8],fmsltrend.pts1)
all.data1<-cbind(all.data[,c(2:36)],fmsl.pts[,3:4])

write.csv(all.data1,'all.data.csv')


##standardized msl
fmsl.std<-apply(X = fmsl.pts[,3:8], MARGIN = 2,FUN = function(x){(x - mean(all.data$avmsl,na.rm=T)) / (2*sd(all.data$avmsl,na.rm=T))})
##standardized trend
fmsltrend.std<-apply(X = fmsltrend.pts1, MARGIN = 2,FUN = function(x){(x - mean(all.data$slr,na.rm=T)) / (2*sd(all.data$slr,na.rm=T))})

data.std1<-cbind(data.std[,2:21],fmsl.std,fmsltrend.std)


write.csv(data.std1,'data.std.csv')


#slr26.2015<-raster(a.26,9)
#slr85.2015<-raster(a.85,9)
#slr26.2050<-raster(a.26,44)
#slr85.2050<-raster(a.85,44)
##crop.extent<-extent(27.44199,51,-33.28764,-0.22713)
#fslr.stack<-stack(slr26.2050,slr85.2050)
#fslr.stack.wio<-crop(fslr.stack,crop.extent)
##names(fslr.stack.wio)<-c('slr26.2050','slr85.2050')
#aviso.msl.wio<-crop(aviso.msl,crop.extent)

#fslr<-as.data.frame(fslr.stack.wio, na.rm=T, xy=T)
#fslr.tree <- createTree(coordinates(fslr[,1:2]))
#inds.fslr <- knnLookup(fslr.tree, newdat=coordinates(all.data[,1:2]), k=1)
#fslr.pts<-fslr[inds.fslr,]

#av.msl<-as.data.frame(aviso.msl.wio, na.rm=T, xy=T)
#avmsl.tree <- createTree(coordinates(av.msl[,1:2]))
#inds.avmsl <- knnLookup(avmsl.tree, newdat=coordinates(all.data[,1:2]), k=1)
#avmsl.pts<-av.msl[inds.avmsl,]
#colnames(avmsl.pts)[3]<-"avmsl"
#all.data1<-data.frame(cbind(all.data,avmsl.pts[,3],fslr.pts[,3:4]))
#colnames(all.data1)[24]<-"avmsl"
#write.csv(all.data1,'all.data.csv')
#all.data<-read.csv('all.data')
#PredictVar1<-all.data[,25:27]
#data.std1<-apply(X = PredictVar1, MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))})
#data.std2<-cbind(data.std,data.std1)
#colnames(data.std2)
#write.csv(data.std2,'data.std.csv')
---------------------------------------------------------------       
# tt<-raster("/Users/maina/Documents/CCVA/Mangrove/Results_FINAL/Dataset_Final/Unnormalised/TravelTime_7k.tif")
#       
# 
#       
#       detided.r<-resample(detided, ldi, method='ngb')
#       elev.r<-resample(elevation, ldi, method='ngb')
#       erosion.r<-resample(erosion, ldi, method='ngb')
#       grav.r<-resample(gravity, ldi, method='ngb')
#       sla.r<-resample(sla, ldi, method='ngb')
#       slope.r<-resample(slope, ldi, method='ngb')
#       tidecm.r<-resample(tidecm, ldi, method='ngb')
#       cdd.r<-resample(ccd, ldi, method='ngb')
#       tx90p.r<-resample(tx90p, ldi, method='ngb')
#       
#       predR<-stack(c(ldi,elev.r,erosion.r, grav.r,sla.r,slope.r,tidecm.r,cdd.r,tx90p.r))
#       
#       if (require(ncdf4)) {	
#         rnc <- writeRaster(predR, filename='all.pred.resample.nc', format="CDF", overwrite=TRUE)   
#       }
#       ####
#       
#       predR<-stack('all.pred.resample.nc')
#       names(predR)<-c('ldi','elev.r','erosion.r', 'grav.r','sla.r','slope.r','tidecm.r','cdd.r','tx90p.r')
#       #predR<-dropLayer(predR, 1)
# 
#       ##vif
#       predR.v <- vifstep(predR, th=10)
#       
#       > predR.v
#       1 variables from the 9 input variables have collinearity problem: 
#         
#         tidecm.r 
#       
#       After excluding the collinear variables, the linear correlation coefficients ranges between: 
#         min correlation ( grav.r ~ ldi ):  0.008634646 
#       max correlation ( tx90p.r ~ cdd.r ):  0.7199057 
#       
#       ---------- VIFs of the remained variables -------- 
#         Variables      VIF
#       1       ldi 1.304522
#       2    elev.r 1.337553
#       3 erosion.r 1.045668
#       4    grav.r 1.400207
#       5     sla.r 2.917363
#       6   slope.r 1.057425
#       7     cdd.r 3.278799
#       8   tx90p.r 5.304706
#       resp<-raster('vci.tif')
#       plot(resp
#            ldi<-raster(predR, 1)
#       respr<-resample(resp, ldi, method='ngb')
#       pred.response<-addLayer(predR, respr)
#       
#       predR<-scale(predR, center = TRUE, scale = TRUE)
#       pred.response<-addLayer(predR, respr)
#   
#       #######
#       #sample 100000 data points randomnly
#       #fit full model 
#       #scale to the extracted data
#       # model average the estimates
#       
#       
#       library(caret)
#       
#       
#     
#       models<-list()
#       r<-list()
#       for (i in seq(1:1000)){
#        r[[i]]<-data.frame(na.omit(extract(pred.response, (randomPoints(pred.response, 10000)))))
#        models[[i]]<-glm(vci  ~ erosion+Gravity+ sla+ slope+   cdd_hist+tx90p, data=r[[i]],family = Gamma(link = "log"))
#       }
#       
#       
#           for i in
#         v[[i]] <- data.frame(na.omit(sampleRegular(pred.response, 10000000)))
#         models[[i]]<-glm(vci  ~ erosion+Gravity+ sla+ slope+   cdd_hist+tx90p   , data =[[i]] , family = Gamma(link = "log"))
#       }
#       
#       scale(predR, center = TRUE, scale = TRUE)
#       
#       dat.random<-sampleRandom(dat, 10000000, na.rm=TRUE)
#       
#       
#       
#       
#       ##all.var.df <- as.data.frame( pred.response, xy=TRUE)
#        
#        
#      # std.fun<-function(x,y){(x - mls(ean(y)) / (2*sd(y))}
#      # all.var.df1<-all.var.df[,c(3,4,7:11)]
#       
#       #all.var.df.std<-scale(all.var.df1, center = TRUE, scale = TRUE)
#       
#       m_glm <- glm(vci  ~ erosion+Gravity+ sla+ slope+   cdd_hist+tx90p   , data = all.var.df, family = Gamma(link = "log"))
#       
#       _________________________________________________________
#       
#       dat<-stack("all.dat.resampled.nc")
#       names(dat)<-c('ldi','elev.r','erosion.r', 'grav.r','sla.r','slope.r','tidecm.r','cdd.r','tx90p.r','vci')
#       vci<-raster('vci.tif')
#       pts<-as.data.frame(vci, na.rm=T, xy=T)
#       elev<-as.data.frame(elevation, na.rm=T, xy=T)
#       eros<-as.data.frame(erosion, na.rm=T, xy=T)
#       grav<-as.data.frame(gravity, na.rm=T, xy=T)
#       ld<-as.data.frame(ldi, na.rm=T, xy=T)
#       sl<-as.data.frame(sla, na.rm=T, xy=T)
#       slp<-as.data.frame(slope, na.rm=T, xy=T)
#       tx90<-as.data.frame(tx90, na.rm=T, xy=T)
#       ccdp<-as.data.frame(ccd, na.rm=T, xy=T)
#       tide.cm<-as.data.frame(tidecm, na.rm=T, xy=T)
#       inund<-as.data.frame(detided, na.rm=T, xy=T)
#       slr.df<-as.data.frame(slr.wio, na.rm=T, xy=T)
#       ndvi.df<-as.data.frame(ndvi, na.rm=T, xy=T)
#       tt.df<-as.data.frame(tt, na.rm=T, xy=T)
        cdd85.df<-as.data.frame(cdd85, na.rm=T, xy=T)
#       elev.tree <- createTree(coordinates(elev[,1:2]))
#       inds.elev <- knnLookup(elev.tree, newdat=coordinates(pts[,1:2]), k=1)
#       elev.pts<-elev[inds.elev,]
#       
#       
#       eros.tree <- createTree(coordinates(eros[,1:2]))
#       inds.eros <- knnLookup(eros.tree, newdat=coordinates(pts[,1:2]), k=1)
#       eros.pts<-eros[inds.eros,]
#       
#       grav.tree <- createTree(coordinates(grav[,1:2]))
#       inds.grav <- knnLookup(grav.tree, newdat=coordinates(pts[,1:2]), k=1)
#       grav.pts<-grav[inds.grav,]
#       
#       ld.tree <- createTree(coordinates(ld[,1:2]))
#       inds.ld <- knnLookup(ld.tree, newdat=coordinates(pts[,1:2]), k=1)
#       ld.pts<-ld[inds.ld,]
#       
#       sl.tree <- createTree(coordinates(sl[,1:2]))
#       inds.sl <- knnLookup(sl.tree, newdat=coordinates(pts[,1:2]), k=1)
#       sl.pts<-sl[inds.sl,]
#       
#       slp.tree <- createTree(coordinates(slp[,1:2]))
#       inds.slp <- knnLookup(slp.tree, newdat=coordinates(pts[,1:2]), k=1)
#       slp.pts<-slp[inds.slp,]
#       
#       tx90.tree <- createTree(coordinates(tx90[,1:2]))
#       inds.tx90 <- knnLookup(tx90.tree, newdat=coordinates(pts[,1:2]), k=1)
#       tx90.pts<-tx90[inds.tx90,]
#       
#       ccdp.tree <- createTree(coordinates(ccdp[,1:2]))
#       inds.ccdp <- knnLookup(ccdp.tree, newdat=coordinates(pts[,1:2]), k=1)
#       ccdp.pts<-ccdp[inds.ccdp,]
#       
#       tide.cm.tree <- createTree(coordinates(tide.cm[,1:2]))
#       inds.tide.cm <- knnLookup(tide.cm.tree, newdat=coordinates(pts[,1:2]), k=1)
#       tide.cm.pts<-tide.cm[inds.tide.cm,]
#       
#       pts.tree <- createTree(coordinates(sector))
#       inds.pts <- knnLookup(pts.tree, newdat=coordinates(pts[,1:2]), k=1)
#       sector.pts<-sector@data[inds.pts,]
#       
#       inund.tree <- createTree(coordinates(inund[,1:2]))
#       inds.inund <- knnLookup(inund.tree, newdat=coordinates(pts[,1:2]), k=1)
#       inund.pts<-inund[inds.inund,]
#       
#       
#       slr.tree <- createTree(coordinates(slr.df[,1:2]))
#       inds.slr <- knnLookup(slr.tree, newdat=coordinates(pts[,1:2]), k=1)
#       slr.pts<-slr.df[inds.slr,]
#       colnames(slr.pts)[3]<-"slr"
#       
#       tt.tree <- createTree(coordinates(tt.df[,1:2]))
#       inds.tt <- knnLookup(tt.tree, newdat=coordinates(pts[,1:2]), k=1)
#       tt.pts<-tt.df[inds.tt,]
#       colnames(tt.pts)[3]<-"travel.time"
        cdd85.tree <- createTree(coordinates(cdd85.df[,1:2]))
        inds.cdd85 <- knnLookup(cdd85.tree, newdat=coordinates(all.data[,2:3]), k=1)
        cdd85.pts<-cdd85.df[inds.cdd85,]
        colnames(cdd85.pts)[3]<-"cdd.rcp85"
      
        all.data1<-cbind(all.data,cdd85.pts[,3])
        colnames(all.data1)[30]<-"cdd.rcp85"
        

#future slr




#       
#       
#     
#       dat<-as.data.frame(cbind(pts,elev.pts[,3],eros.pts[,3],grav.pts[,3],ld.pts[,3],sl.pts[,3],slp.pts[,3],tide.cm.pts[,3],tx90.pts[,3],ccdp.pts[,3],inund.pts[,3],sector.pts))
#       colnames(dat)<-c('x','y','vci','elevation','erosion','gravity','ldi','sla','slope','tidecm','tx90','ccd','inundation')
#       write.csv(dat,"raw.data.csv")

#sector<-readOGR(dsn='/Users/maina/Downloads/Archive','wio_sectors')
#       
#       coordinates(dat)<-~x+y
#       sect.data <- data.frame(sector)
#       
#       proj4string(dat) <- proj4string(sector)
#       #function over from package sp
#       pts.sect <- data.frame(xx=over(dat, as(sector,"SpatialPolygons")))
#       #dat1<-read.csv("raw.data.csv")
#       
#       
#       #dat1<-read.csv("all.data.csv")
#       #all.data<-data.frame(cbind(dat1,tt.pts[,3]))
#       #colnames(all.data)[26]<-"travel.time"
#       #colnames(all.data)[24]<-"ndvi"
#       # all.data$travel.time<-tt.pts$travel.time
#       
#       write.csv(all.data,"all.data.csv")
#       all.data<-read.csv("all.data.csv")
#       
#       #transformationa ccording to Daniel Zimprich
#       respVar<-as.data.frame(all.data[,26])
#       scaling <- function(y){ 
#         a<-min(y)
#         b<-max(y)
#         N<-length(y)
#         y_a<-y-a
#         b_a<-(b-a)
#         yscaled<-((y_a*(N-1))/(b_a*N))+(1/(2*N))
#         return(yscaled) 
#       }
#       
#       resp.std<-data.frame(apply(X = respVar, MARGIN = 2,
#                                  FUN = function(x) scaling(x)) )
#      
#      
#   
#       source("functions_analyses_interaction_glmmADMB.R")
#       PredictVar<-all.data[,c(5,6,8:17,25,27)]
#       #standasrdize predictors
#       data.std<-apply(X = PredictVar, MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))})
#       data.std<-cbind(all.data[,c(7,26)],resp.std,all.data[,19:20],data.std)
#       colnames(data.std)[3]<-"ndvi.std"
#       write.csv(data.std,"data.std.csv")
#       data.std<- read.csv("data.std.csv")
      

      #sample the data for model fitting
     


source("functions_analyses_interaction_glmmADMB.R")

all.data$vci<-round(all.data$vci, 2)
        
PredictVar<-all.data[,c(5:14,22,24,25)]
#PredictVar1<-all.data[,c(25:29)]
    
      train.index <- createDataPartition(all.data$vci, times=1, p = .10, list = FALSE)
      train <- all.data[ train.index,]
      #balance  <- all.data[-train.index,]
#standasrdize data 
PredictVar1<-train[,c(5:14,22:38)]   
data.std1<-data.frame(apply(X = PredictVar1, MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}))
data.std2<-data.frame(cbind(train$vci,as.character(train$layer),as.character(train$County),data.std1 ))
colnames(data.std2)[1]<-"vci"
colnames(data.std2)[2]<-"sector"
colnames(data.std2)[3]<-"country"
data.std2$raneff<-factor(1:nrow(data.std2))
      #train$pos <- numFactor(train$x, train$y)
      #train$group <- factor(rep(1, nrow(train)))
    f1 <- glmmTMB(vci ~ elevation +erosion + gravity + (1|layer)+mat(pos + 0 | layer), data=train)
    f2 <- fitme(vci ~ elevation +erosion + gravity+ (1|layer)+Matern(1 | x + y), data=train, family =Gamma())
    f3 <- glmmTMB(vci ~ elevation +erosion + gravity + (1|layer)+mat(pos + 0 | layer), data=train)
    # #1. Create list with all possible combinations between predictors 
      vifPredCombinations  <-  list()
      varnames<-colnames(PredictVar)#
      
      maxCombs  <-  getMaximumNOfCombs(varnames)
      for(j in 1:maxCombs) {
        vifPredCombinations  <-  append(runPredCombinations(j, varnames), vifPredCombinations)
      }
      
      ##2. filter the combinations above with VIF<1.5
      vifPredCombinations_new<- c()
      for(con in vifPredCombinations){
        r <- subset(PredictVar, select = con)
        conClasses   <-  unique(sapply(r, class))
        numOfClasses  <-  length(conClasses)
        twoNCols <- ncol(r)==2
        numDf <- r[sapply(r,is.numeric)]
        zeroNumDf<-ncol(numDf)==0
        numeriCols<- ncol(numDf)>1
        
        if(length(con)<2){
          vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
        }else{
          if (zeroNumDf) { 
            vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
          }
          if (numOfClasses==2 && twoNCols) { 
            vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
          }
          if(numeriCols && max(vif(numDf)["VIF"])<=1.2){##vif cutoff
            vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
          }
          next
        }
      }
  #sector
      modelText.vci<-mclapply(vifPredCombinations_new, prepareModelText,train )
      modList.vci<-mclapply(modelText.vci, evalTextModel,mc.preschedule = TRUE, mc.cores=12)
      #modList.vci<-lapply(modelText.vci, evalTextModel)
      modList.vci<-modList.vci[!sapply(modList.vci, is.null)]
      
      #modList_rslope<-lapply(modelText, evalTextModel)
      #modList_rslope<-modList_rslope[!sapply(modList_rslope, is.null)] 
      modelSel.ndvi<-model.sel(modList.ndvi, rank.args = list(REML = FALSE), extra =c(AIC, BIC))#sector interaction
      write.csv(modelSel.ndvi, 'modelSel.ndvi.csv')
      top.model.ndvi<-get.models(modelSel.ndvi, subset=delta<2)
      top.model.ndvi.ave<-model.avg(top.model.ndvi) 
      #modelAvg.ndvi<-model.avg(modelSel.ndvi)  
      
      ##vci
      modelText.vci<-mclapply(vifPredCombinations_new, prepareModelText,train )
      modList.vci<-lapply(modelText.vci, evalTextModel)
      modList.vci<-modList.vci[!sapply(modList.vci, is.null)]
      modelSel.vci<-model.sel(modList.vci, rank.args = list(REML = FALSE), extra =c(AIC, BIC))#sector interaction
      top.model.vci<-get.models(modelSel.vci, subset=delta<2)
      top.model.vci.ave<-model.avg( top.model.vci)
      #modelAvg.vci<-model.avg(modelSel.vci)  
      write.csv(modelSel.vci, 'modelSel.vci.csv')
      
      save.image("dat.RData")
      load("dat.RData")
      
      #simulate confn intervals
      #ndvi.coef<-lme4::bootMer(top.model.ndvi$`1`,FUN=function(x) unlist(fixef(x)), nsim=2)
      >confint(top.model.ndvi$`1`)
      confint(top.model.vci$`16`)
      #calculate means and SD of the original data
      scale.var.mean<-apply(PredictVar, 2, mean)
      scale.var.sd<-apply(PredictVar, 2, 2*sd)
      #vci ~ elevation + erosion + gravity + slope + slr + avmsl + (1 |      layer)
      #ndvi.std ~ elevation + erosion + gravity + ldi + slope + ccd +      slr + (1 | layer)
      

      #                                 2.5 %      97.5 %    Estimate
      #cond.(Intercept)                1.04482006  1.28023675  1.16252841
      #cond.elevation                  0.46268987  0.49635973  0.47952480
      #cond.erosion                   -0.31722150 -0.29180494 -0.30451322
      #cond.gravity                   -0.18845035 -0.15115445 -0.16980240
      #cond.ldi                        0.03415236  0.08888497  0.06151866
      #cond.slope                     -0.03816264 -0.01181647 -0.02498956
      #cond.ccd                       -0.57908949 -0.45524738 -0.51716844
      #cond.slr                       -0.19299125 -0.15949387 -0.17624256
      
      
      
     
     ##calculate relative weights
      
      
      
      
      
       
    
      #modelSel_county<-model.sel(modList_county, rank.args = list(REML = FALSE), extra =c(AIC, BIC))
      #write.csv(modelSel_county, 'modList_county.csv')
      
      #modelSel_rslope<-model.sel(modList_rslope, rank.args = list(REML = FALSE), extra =c(AIC, BIC))
      
      #write.csv(modelSel_rslope, 'modList_rslope.csv')
      

      #top.model<-get.models(modelSel, subset=delta<2)
      save.image("dat.RData")
      #top.model
      #library(MuMIn)
     # modelAvg.vci<-model.avg(modelSel.vci)
      
    #extract averaged coefficients
      Coef.Conditional<-data.frame(modelAvg.vci$coefficients[2,])
      coefTable<-cbind(data.frame(confint(modelAvg.vci, level = 0.975),confint(modelAvg.vci,level = 0.95), Coef.Conditional))
      colnames(coefTable)<-c("conf2.5","conf97.5","conf5","conf95","StdCoef")
      coefTable$var<-rownames(coefTable)
      
      write.csv(coefTable, "coefTable.vci.csv")
      coefTable<-read.csv("coefTable.vci.csv")
      
      ggplot(coefTable, aes(var,StdCoef)) + geom_point(stat="identity") +  ggtitle("Model Averaged coefficients")+ coord_flip() + geom_hline(yintercept=0)+geom_errorbar(data = coefTable, aes(ymin = conf5, ymax = conf95),colour = 'black', width = 0.4) + xlab("standardized coefficients")
      
      mod.pplots.ndvi<-get.models(modelSel.ndvi, subset=delta<410)
      
    library(MuMIn)
      library(sjPlot)
      theme_set(theme_sjplot())
      library(sjmisc)
      plot_model(top.model)
      tab_model(top.model)
      
      library(ncf)
      library(ade4)
      library(ecodist)
      
      
      X <- getME(top.model.vci$`16`, "X")
      beta <- fixef(top.model.vci$`16`)$cond
      beta_X <- sweep(X,MARGIN=2,STATS=beta,FUN="*")
      p_resid <- sweep(beta_X,MARGIN=1,STATS=residuals(top.model.vci$`16`),FUN="+")
      par(mfrow=c(2,3))
      for (i in 2:7) {
        plot(X[,i],p_resid[,i],xlab=colnames(X)[i],ylab="partial residuals")
      }
   
  #plot ndvi
  top.model.vci$`12`
  write.csv(p_resid,"p_resid_top.vci.csv")  
  colnames(p_resid)
  p_resid<-as.data.frame(p_resid[,2:7])
  sub.data.std<-subset(train, select=c(elevation,erosion,gravity,slope,ccd,slr,layer,County))
    
  sub.data.std$type<-"predictors"
  p_resid$type<-"partial.residuals"
  library(reshape2)
  library(ggplot2)
  sub.data.std.1<-melt(sub.data.std, id=c('type','layer','County'),value.name = "predictors")
  p_resid.t<-melt(p_resid, id=c('type'),value.name = "partial.residuals")
  dataforgg.vci<-cbind(sub.data.std.1[,2:5],p_resid.t[,3])
  colnames(dataforgg.vci)<-c('layer', 'Country','variable', 'predictors','partial.residuals')
  write.csv(dataforgg.ndvi,"dataforgg.ndvi.csv")
  
   dataforgg.ndvi<-read.csv("dataforgg.ndvi.csv")
  
   ggplot(dataforgg.vci, aes(predictors,partial.residuals)) + geom_point(size=4, alpha = 0.5, aes(colour = Country)) + facet_wrap(~variable, scales = "free") + 
    theme_bw() + 
    scale_colour_manual(values = c("Kenya" = "green","Tanzania" = "blue", "Mozambique" = "pink","Madagascar" = "black")) +
    theme(panel.grid.minor.x = element_line(colour='grey94'),
          panel.grid.minor.y = element_line(colour='grey94'),
          panel.grid.major.x = element_line(colour='lightgrey'),
          panel.grid.major.y = element_line(colour='lightgrey'))

  #vci
  top.model.vci$`9`
  write.csv(p_resid,"p_resid_top.vci.csv")  
  colnames(p_resid)
  p_resid<-as.data.frame(p_resid[,2:8])
  sub.data.std<-subset(train, select=c(elevation,erosion,gravity,ldi,slope,tx90,slr,layer,County))
  
  sub.data.std$type<-"predictors"
  p_resid$type<-"partial.residuals"
  library(reshape2)
  library(ggplot2)
  sub.data.std.1<-melt(sub.data.std, id=c('type','layer','County'),value.name = "predictors")
  p_resid.t<-melt(p_resid, id=c('type'),value.name = "partial.residuals")
  dataforgg.vci<-cbind(sub.data.std.1[,2:5],p_resid.t[,3])
  colnames(dataforgg.vci)<-c('layer', 'Country','variable', 'predictors','partial.residuals')
  write.csv(dataforgg.vci,"dataforgg.vci.csv")
  
  dataforgg.vci<-read.csv("dataforgg.vci.csv")
  
  
  ggplot(dataforgg.vci, aes(predictors,partial.residuals)) + geom_point(size=4, alpha = 0.5, aes(colour = Country)) + facet_wrap(~variable, scales = "free") + 
    theme_bw() + geom_line(y = fitted(top.model.vci$`9`), alpha = 0.5, lwd = 2, color = "red")
    scale_colour_manual(values = c("Kenya" = "green","Tanzania" = "blue", "Mozambique" = "pink","Madagascar" = "black")) +
    theme(panel.grid.minor.x = element_line(colour='grey94'),
          panel.grid.minor.y = element_line(colour='grey94'),
          panel.grid.major.x = element_line(colour='lightgrey'),
          panel.grid.major.y = element_line(colour='lightgrey'))
  
  
  
  library(ggeffects)
  ggpredict(top.model, c("vci", "elevation")) %>% plot(rawdata = TRUE, jitter = .01)
  model2<-get.models(modelSel, subset=4)
  
  ###################

  

##### Mantel's Test to test for overal autocorrelation (numerical representation)
      library(ade4)
      library(ecodist)
      
      #Generate distance matrices
      Site.dists <- dist(cbind(all.data$x, all.data$y))
      vci.dists <- dist(all.data$vci) 
      vci.distsresiduals <- dist(resid(top.model$`1`)) 
      
      #Mantel Test
      #mantel.results=mantel.rtest(BP.dists, Site.dists, #nrepet=999)  
      
      #Plot Mantel Test
      mantel1=mgram(vci.dists, Site.dists)
      mantel2=mgram(vci.distsresiduals, Site.dists)
    
      
      
      pdf("MantelTest.pdf")
      par(mfrow=c(2,2), mai=c(0.5,1,0.2,0.2), mgp=c(0.5,0.6,0), mar=c(1.5,2.5,1,0.2), oma=c(2.5,1.5,0,0.5), bg="white", fg="black", cex.axis=1.05)
      plot(mantel1, ylim=c(-1,1), ann=F);abline(0,0, lty=2,col="black",lwd = 1.5)
      legend("topleft", "(A) Raw data", bty="n", cex=1.1)
      mtext("Mantel r", side=2, adj=0.5, line=2, cex=1.1) 
      plot(mantel2, ylim=c(-1,1), ann=F);abline(0,0, lty=2,col="black",lwd = 1.5)
      legend("topleft", "(B) Residuals.model.1", bty="n", cex=1.1)
      mtext("Distance", side=1, adj=0.5, line=2, cex=1.1)  
      plot(mantel3, ylim=c(-1,1), ann=F);abline(0,0, lty=2,col="black",lwd = 1.5)
      legend("topleft", "(C) Residuals.model.2", bty="n", cex=1.1)
      mtext("Distance", side=1, adj=0.5, line=2, cex=1.1) 
      mtext("Mantel r", side=2, adj=0.5, line=2, cex=1.1) 
      dev.off()
      
      
      
      
      
      
      
      
      
      
      
      CorRes.vci<- spline.correlog(all.data$x, all.data$y, all.data$vci, xmax = FALSE,resamp=1000,latlon=TRUE, quiet=TRUE, filter=TRUE, na.rm=T)
      CorRes.vci.res <- spline.correlog(all.data$x, all.data$y, as.vector(resid(m1)),xmax = FALSE, resamp=1000,latlon=TRUE, quiet=TRUE, filter=TRUE, na.rm=T)
      
      
      data<-data.std
      
      data <- data[complete.cases(data), ]
      # create a function suitable for boot that will return the goodness of fit
      # statistics testing models against the full original sample.
      compare <- function(orig_data, i){
        # create the resampled data
        train_data <- data[i, ]
        test_data <- data # ie the full original sample
        m1<-glmmTMB(vci ~elevation+erosion+gravity+ slope+slr+avmsl+ (1|layer), family=beta_family(link = "logit"),data=train_data)
       
          # predict the values on the original, unresampled data
          m1p  <- predict(m1, newdata = test_data)
        
          # return a vector of 6 summary results
          results <- c(m1r = R2(inv.logit(m1p), test_data$bleach_intensity))
          
          return(results)
        }
          
          
          # perform bootstrap
          Repeats <- 100
          res <- boot(data, statistic = compare, R = Repeats)
          tests <- as.data.frame(res)
          
          names(tests) <- c("m1r")
          tests
      
      
      
      
      
      #BRT
      library(gbm)
      library(rpart)
      library(tree)
      library(PerformanceAnalytics)
          
      
      vif.step(PredictVar)
      "elevation", "erosion" ,  "gravity"  , "ldi"  ,     "sla"   ,    "slope"  ,   "tidecm"  ,  "tx90",     
       "ccd"   ,    "detided" ,  "slr" 
      
      all.data.mc<-all.data[,c(5:14,23)]
      library(caret)
      library(randomForest)
      train.index <- createDataPartition(all.data$layer, p = .7, list = FALSE)
      train <- all.data.mc[ train.index,]
      test  <- all.data.mc[-train.index,]
      
      #feature selection
      control <- rfeControl(functions = rfFuncs,
                            method = "repeatedcv",
                            repeats = 3,
                            verbose = FALSE)
      
      outcomeName<-'vci'
      predictors<-names(train)[!names(train) %in% outcomeName]
      
      vci_Pred_Profile <- rfe(train[,predictors], train[,outcomeName],
                               rfeControl = control)
      
      
      
      ##transorm ndvi to beta regression
      
      
      
      
      ###PLOTS
      ggplot(all.data.agg, aes( elevation,gravity, label=Group.1, color=Group.2)) + geom_point() + geom_text(check_overlap = TRUE)
      all.data.agg <- aggregate(x=all.data, by = list(all.data$layer,all.data$County), FUN="mean")
      all.data.agg = all.data.agg[-1,]
      
      
      ##spca
      install.packages("adegenet")
      library(adegenet)
      data.std.spca<-cbind(all.data[,3:4],data.std[,6:16]
      write.csv()
      
      
     ##Normalizse data
      #If e < c < i, then the membership function is incr                   easing from e to i. If i < c < e, then the membership function is decreasing from i to e.
      library(QCA)
      #reinforcing<-data.frame(PredictVar[,c("elevation") ])
      increasing<- function(x) {
        minT<-quantile(x, probs = seq(0, 1, 0.01))[[2]]
        maxT<-quantile(x, probs = seq(0, 1, 0.01))[[100]]
        #minT<-min(x)
        #maxT<-(max(x))
        mtxt<-paste0('calibrate(x, thresholds = "e=', minT,',c=',(maxT+minT)/2,',i=', maxT,'",logistic = FALSE)',  collapse="")
        reinf<-eval(parse(text = mtxt))
        return(reinf)
      }
     #p.exposureReinf <- apply( reinforcing, 2,  increasing)
      #colnames(p.exposureReinf)[1]<-'elevation'
      ##"erosion"     "gravity",  "travel.time" "slr"
     #reducing<-PredictVar[,c("erosion","gravity", "slr","slope","ccd" ) ]
      decreasing <- function(x){
        minT<-quantile(x, probs = seq(0, 1, 0.01))[[2]]
        maxT<-quantile(x, probs = seq(0, 1, 0.01))[[100]]
        #minT<-min(x)
       #maxT<-max(x)
        mtxt<-paste0('calibrate(x, thresholds = "e=', maxT,',c=',(maxT+minT)/2,',i=',minT ,'",logistic = FALSE)',  collapse="")
        reduc<-eval(parse(text = mtxt))
        return( reduc)
      }
      #p.exposureReduce <- apply( reducing, 2,  decreasing)
      
      
##slr
#min.slr<-min(reducing$slr)
#max.slr<-max(reducing$slr)    
#slr.std<-calibrate(reducing$slr, thresholds = "e=5.041236,c=3.874998,i=2.708759",logistic = FALSE)
      
  #combine the data spreadsheets  
  p.exposure<-data.frame(cbind(p.exposureReinf,p.exposureReduce))
      
   
  #ndvi 
  #Conditional model:
  (conditional average) 
  Estimate Std. Error z value Pr(>|z|)    
  (Intercept)  1.162528   0.060056   19.36  < 2e-16 ***
    elevation    0.479525   0.008589   55.83  < 2e-16 ***
    erosion     -0.304513   0.006484  -46.96  < 2e-16 ***
    gravity     -0.169802   0.009514  -17.85  < 2e-16 ***
    ldi          0.061519   0.013963    4.41 1.05e-05 ***
    slope       -0.024990   0.006721   -3.72 0.000201 ***
    ccd         -0.517168   0.031593  -16.37  < 2e-16 ***
    slr         -0.176243   0.008545  -20.62  < 2e-16 ***
   #vci   
  #Conditional model: vci ~ elevation + erosion + gravity + slope + slr + avmsl + (1 |      layer)
              Estimate Std. Error z value Pr(>|z|)    
  (Intercept)  0.294588   0.042541   6.925 4.37e-12 ***
    elevation    0.074646   0.004423  16.875  < 2e-16 ***
    erosion     -0.115132   0.003926 -29.325  < 2e-16 ***
    gravity     -0.088406   0.005546 -15.942  < 2e-16 ***
    slope       -0.016215   0.003754  -4.320 1.56e-05 ***
    slr         -0.049257   0.005044  -9.765  < 2e-16 ***
    avmsl       -0.341969   0.019268 -17.748  < 2e-16 ***
      
 w.elev<- 
 w.erosion <-
 w.grav<-
  w.slr<-
  w.ldi<-
  w.slope<-
  w.ccd<-
  w.  avmsl <- 
    
# tot.weight<-sum(w.elev,w.erosion,w.grav,w.tx90,w.slr,w.tt,w.ldi,w.slope,w.ccd)
 
tot.weight1<-sum(w.elev,w.erosion,w.grav,w.slr,w.slope,w.ccd)
 
 wt.elev<-w.elev/tot.weight1*1
 wt.erosion <-w.erosion/tot.weight1*1
 wt.grav<- w.grav/tot.weight1*1
 #wt.tx90<-w.tx90/tot.weight*1
 wt.slr<-w.slr/tot.weight1*1
 #wt.tt<-w.tt/tot.weight*1
 #wt.ldi<-w.ldi/tot.weight*1
 wt.slope<-w.slope/tot.weight*1
 wt.ccd<-w.ccd/tot.weight1*1
 
 #wt<-c(0.2555277,0.260458,0.1315283,0.120148,0.1164822,0.0632932,0.02368855,0.02887404)
 wt<-c(0.3246059, 0.3149046,0.07791313,0.1257441,0.02146122,0.1313389)
 
 
 #names(wt)<-c('elevation','erosion','gravity','tx90','slr','travel.time','ldi','slope')
 names(wt)<-c('elevation','erosion','gravity','slr', 'slope','ccd')
 
 
p.exposure<-p.exposure[names(wt)]
p.exposure.w<- sweep(p.exposure,MARGIN=2,STATS=wt,FUN="*")
p.exposure$Exposure<-rowSums(p.exposure.w)
 #p.exposure.f<-cbind(all.data[,c("layer","County")], p.exposure.w)
  
p.exposure<-cbind(all.data[,c("layer")], exposure)
colnames(p.exposure)[1]<-"sector"
exposure.agg.mean = aggregate(p.exposure,by = list(all.data$layer),FUN = max)
exposure.agg.median = aggregate(exposure,by = list(all.data$layer),FUN = median)

rownames(exposure.agg.median)<-exposure.agg.median[,1]
rownames(exposure.agg.mean)<-exposure.agg.mean[,1]

exposure.agg.2<-exposure.agg.median[,-c(1)]
exposure.agg.3<-exposure.agg.mean[,-c(1)]

exposure.agg.2<-as.matrix(exposure.agg.2)
exposure.agg.3<-as.matrix(exposure.agg.3) 
 
 
 
 
 #names(exposure.agg.med)[1]<-'Sector'
 #names(exposure.agg.med)[2]<-'Country'
 #exposure.agg.med<-exposure.agg.med[,-c(3,4)]
 #write.csv(exposure.agg,"exposure.agg.csv
 #exposure.agg.1<-melt( exposure.agg, id=c('Sector','Country'),value.name = "exposures")
 
# library(ggplot2)
 #cols <- rev(rainbow(7)[-7])
# ggplot(exposure.agg.1, aes(variable,Sector)) +
 #  geom_raster(aes(fill = exposures), interpolate = FALSE) +
 #  scale_fill_gradientn(colours = cols)
 
 library(devtools)
 install_github("jokergoo/ComplexHeatmap")
 library(ComplexHeatmap)
 
 rownames(exposure.agg)<-exposure.agg[,1]
 #exposure.agg.2<-exposure.agg[,-c(1:2)]
 #exposure.agg.2<-as.matrix(exposure.agg.2)
 #Heatmap(exposure.agg.2)
 library(circlize)
 library(seriation)
 col_fun = colorRamp2(c(0, 0.5, 1), c("green", "white", "red"))
 col_fun(seq(-3, 3))
 
 o = seriate(max(exposure.agg) - exposure.agg, method = "BEA_TSP")
 
 Heatmap(max(exposure.agg) - exposure.agg, name = "exposure", 
         row_order = get_order(o, 1), column_order = get_order(o, 2))
 
 
 ht<-Heatmap(exposure.agg.2, row_dend_reorder = TRUE,name = "exposure",row_km = 4, border = TRUE)    
 
 ht = Heatmap(exposure.agg.2, name = "exp")
 dim(ht)
 
 
 ht = Heatmap(exposure.agg.3, name = "exp", row_km = 2, column_km = 3,
              col = colorRamp2(c(0, 0.5, 1), c("green", "white", "red")),
              right_annotation = rowAnnotation(expo = exposure.agg.2[,7]))
 
 
 
 
 ############post hoc
 drp.ndvi<-drop1(top.model.ndvi$`1`)
 
 drp.vci<-drop1(top.model.vci$`16`)
 
 

 ##########predictions
 library(ggeffects)
 library(ggplot2)
 library(sjPlot)
 
 
 #elev.ndv<-plot_model(top.model.vci$`16`, type = "pred", terms = c("elevation","layer")))
 
#data.std$cdd.rcp85<-(all.data$cdd.rcp85 - mean(all.data$cdd.rcp85,na.rm=T)) / (2*sd(all.data$cdd.rcp85,na.rm=T))
#data.std1<-data.std[,2:34]
#write.csv(data.std1,'data.std.csv')
 new.data1<-data.std[,c(5,9,10,11,14,32,27,12,34)]
 colnames(new.data1)<-c("layer","elevation" ,"erosion" ,  "gravity"  , "slope",  "slr", "avmsl", 'ldi','ccd')
 fut.vci.logit<- predict(top.model.vci$`16`, newdata = new.data1)
 fut.ndvi.logit<- predict(top.model.ndvi$`1`, newdata = new.data1)
 fut.vci<-inv.logit(fut.vci.logit)
 fut.ndvi<-inv.logit(fut.ndvi.logit)
 future.exposure.vci<-decreasing(fut.vci)
 future.exposure.ndvi<-decreasing(fut.ndvi)
 
 ##new data
 
 
 
 
 
 #marginal predictions using visreg
 library(visreg)
 elev.pred<-visreg(top.model.vci$`16`,xvar="elevation",type="conditional")
 elev.pred.df<-elev.pred$fit
#elev.pred <- ggpredict(top.model.vci$`16`, terms = "elevation",condition = c(avmsl = 0.10419))

 ##manual prediction
 new.data2<-new.data1[,2:9]
 med<-apply(new.data2, 2, median)
 
 
 
 
 
 ##BRT methods
hist(data.std$ndvi)'

source("brtfunctions.R")

library(here)
library(PerformanceAnalytics)
 library(MASS)
 library(tree)
 library(gbm)
 library(rpart)
 library(dismo)
 
#Random forest methods
library(randomForest)
library(party)
library(dplyr)
library(pdp)
library(vip)
library(caret)
library(mlr)
libray(rpart)
library(rpart.plot)

package.list = c("party", "mlr", "rpart", "rpart.plot")
tmp.install = which(lapply(package.list, require, character.only = TRUE)==FALSE)
if(length(tmp.install)>0) install.packages(package.list[tmp.install], repos = "http://cran.us.r-project.org")
lapply(package.list, require, character.only = TRUE)

#to use cforest::party given correlation in data
#stratify using response
#add sector as a predictor
#predict and extract by sector
#use pdp for partial plots

setwd('/Volumes/Data/Projects/CCVA/Mangrove/Results_FINAL/Dataset_Final/Unnormalised')

all.data<-read.csv("all.data.csv")


#sample data geographically

all.data$vci<-round(all.data$vci, 2)

#[1] "elevation"   "erosion"     "gravity"     "ldi"         "sla"         "slope"      
#[7] "tidecm"      "tx90"        "ccd"         "detided"     "slr"         "travel.time"
#[13] "avmsl"

##vci
train.index <- createDataPartition(c(all.data$layer), times=2, p = .30, list = FALSE)
train <- all.data[train.index[,1],]
test<-all.data[ train.index[,2],]

PredictVar<-all.data[,c(5:14,22,24,25)]
predictors<-PredictVar[,-c(10,12)]
all.data.vci<-cbind(all.data[,"vci"],predictors)
all.data.ndvi<-cbind(all.data[,"ndvi"],predictors)
colnames(all.data.ndvi)[1]<-"ndvi"
colnames(all.data.vci)[1]<-"vci"

#create random data for threshold
random1<-sample(1:10,nrow(train), replace=TRUE)
random2<-sample(1:100,nrow(train), replace=TRUE)

#implementation 1
vci.rf.1<-randomForest(vci ~elevation+erosion+gravity+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl,importance=TRUE,ntree=2000,nodesize=5,data=train)
#vci.rf.2<-randomForest(vci ~elevation+erosion+travel.time+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl+ random1+random2,importance=TRUE,ntree=2002,nodesize=5,data=train)
#vci.rf.3<-randomForest(vci ~detided+erosion+gravity+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl+ random1+random2,importance=TRUE,ntree=2002,nodesize=5,data=train)
#vci.rf.4<-randomForest(vci ~detided+erosion+travel.time+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl+ random1+random2,importance=TRUE,ntree=2002,nodesize=5,data=train)

ndvi.rf.1<-randomForest(ndvi ~elevation+erosion+gravity+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl,importance=TRUE,ntree=2002,nodesize=5,data=train)
#ndvi.rf.1<-randomForest(vci ~elevation+erosion+travel.time+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl+ random1+random2,importance=TRUE,ntree=2002,nodesize=5,data=train)
#ndvi.rf.1<-randomForest(vci ~detided+erosion+gravity+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl+ random1+random2,importance=TRUE,ntree=2002,nodesize=5,data=train)
#ndvi.rf.1<-randomForest(vci ~detided+erosion+travel.time+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl+ random1+random2,importance=TRUE,ntree=2002,nodesize=5,data=train)

varImpPlot(vci.rf.1)
varImpPlot(ndvi.rf.1)

##use cforest
vci.crf.1<-cforest(vci ~elevation+erosion+gravity+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl,data=train,controls=cforest_unbiased(ntree=500, mtry=4))
#vci.crf.2<-cforest(vci ~elevation+erosion+gravity+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl+ layer, random1+random2,data=train,controls=cforest_unbiased(ntree=2000, mtry=4))

ndvi.crf.1<- cforest(ndvi ~elevation+erosion+gravity+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl,data=train,controls=cforest_unbiased(ntree=2000, mtry=4))
#ndvi.crf.2<-cforest(ndvi ~elevation+erosion+gravity+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl+ layer, random1+random2,data=train,controls=cforest_unbiased(ntree=2000, mtry=4))

#var importance
library(edarf)

variable_importance(vci.crf.1, nperm = 2, data = train)
vip(ndvi.rf.1, bar = FALSE, horizontal = FALSE, size = 1.5) 

#partial plot
partialPlot(ndvi.rf.1, pred.data = train, x.var = "erosion") 

p1 <- partial(ndvi.rf.1, pred.var = c("elevation", "ndvi"), plot = TRUE, chull = TRUE)



library(mlbench)
library(forestFloor)
library(AUC)
library(pRF)

#statistical significance
#p.test<-pRF(response=as.numeric(train$ndvi),predictors<-predictors,nperms=20,mtry=3,type="regression,alpha=0.05)


regr.task = makeRegrTask(data = all.data.ndvi, target = "ndvi")
regr.lrn = makeLearner("regr.cforest")
regr.mod = train(regr.lrn, regr.task, subset = train.index[,1])
regr.mod

task.pred = predict(regr.mod, task = regr.task, subset = train.index[,2])
performance(task.pred, measures=list(mse, rsq))
plot(as.data.frame(task.pred)[,c("truth","response")])

listMeasures(regr.task)

getParamSet(regr.lrn)


library(parallel)
library(parallelMap)

# Initialize paralelllization
parallelStartSocket(cpus = 10)


ps = makeParamSet(
  makeIntegerParam("ntree", lower = 1, upper = 100),
makeIntegerParam("mtry", lower =  1, upper = 10)
      )

ctrl = makeTuneControlGrid(resolution = 10)
rdesc = makeResampleDesc("CV", iters = 5)

tune.cforest = tuneParams(regr.lrn, task = regr.task, resampling = rdesc, par.set = ps, control = ctrl)
parallelStop()


#-----------------------
vci.rf.1<-randomForest(vci ~elevation+erosion+gravity+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl,importance=TRUE,ntree=2000,nodesize=5,data=train)
vci.crf.1<-cforest(vci ~elevation+erosion+gravity+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl,data=train,controls=cforest_unbiased(ntree=500, mtry=4))
ndvi.rf.1<-randomForest(ndvi ~elevation+erosion+gravity+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl,importance=TRUE,ntree=2002,nodesize=5,data=train)
ndvi.crf.1<- cforest(ndvi ~elevation+erosion+gravity+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl,data=train,controls=cforest_unbiased(ntree=2000, mtry=4))

p1 <- partial(vci.rf.1, pred.var = c("elevation", "erosion"), plot = TRUE, chull = TRUE)
p2 <- partial(boston_rf, pred.var = c("elevation", "erosion"), plot = TRUE, chull = TRUE,
palette = "magma")
grid.arrange(p1, p2, nrow = 1)  # Figure 7

