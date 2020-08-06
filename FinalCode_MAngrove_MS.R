#Random forest methods
rm(list=ls())
library(randomForest)
library(party)
library(dplyr)
library(pdp)
library(vip)
library(caret)
library(mlr)
library(rpart.plot)
library(raster)
library(data.table)
library(QCA)
library(rgdal)
library(sf)
#library(partykit)
library(parallel)
library(usdm)


#setwd("D:/mainambui/mangrove_ccva/")
setwd('/Volumes/Data/Projects/CCVA/Mangrove/Results_FINAL/Dataset_Final/Unnormalised/') 
load("rforestTune.all_gcb_rev1.RData")
#rforestTune.all_gcb_rev1.RData
#coordinates(all.data)<-~x+y
#typol<-readOGR(dsn='/Volumes/Data/Projects/CCVA/Mangrove/Sectors/',"mangrovewithtypology")
#crs(all.data)<-crs(typol)
#polys<-st_as_sf(typol)
#pts<-st_as_sf(all.data)
#all.data_rev <- as.data.frame(st_join(pts, polys, join = st_nearest_feature))


#all.data_rev$FIN_TYP[all.data_rev$FIN_TYP==7] <-'Non_filter_Arheic'
#all.data_rev$FIN_TYP[all.data_rev$FIN_TYP==1] <-'small_deltas'
#all.data_rev$FIN_TYP[all.data_rev$FIN_TYP==2] <-'tidal_systems'
#all.data_rev$FIN_TYP[all.data_rev$FIN_TYP==3] <-'lagoons'
#all.data_rev$FIN_TYP[all.data_rev$FIN_TYP==51] <-'largeRiversTidalDeltas'

#all.data<-read.csv("all.data.csv")
#all.data$FIN_TYP<-as.factor(all.data_rev$FIN_TYP)

#colnames(all.data)
#colnames(all.data)[40]<-"formation"
#all.data<-all.data[,2:40]
#write.csv(all.data,"all.data.csv")
#save.img("rforestTune.all_gcb_rev1.RData")
up

#train.index <- createDataPartition(c(all.data$layer), times=2, p = .30, list = FALSE)
train <- all.data[train.index[,1],]
test<-all.data[ train.index[,2],]
#PredictVar<-all.data[,c(5:14,22,24,25)]
PredictVar<-all.data[,c(5:14,22,24,25,39)]#new
predictors<-PredictVar[,-c(7,8,10,12,14)]
all.data.vci<-cbind(all.data[,"vci"],predictors)
all.data.ndvi<-cbind(all.data[,"ndvi"],predictors)
colnames(all.data.ndvi)[1]<-"ndvi"
colnames(all.data.vci)[1]<-"vci"

regr.task = makeRegrTask(data = all.data.ndvi, target = "ndvi")
regr.task.v = makeRegrTask(data = all.data.vci, target = "vci")

regr.lrn = makeLearner("regr.cforest")
regr.lrn.v = makeLearner("regr.cforest")

regr.mod = train(regr.lrn, regr.task, subset = train.index[,1])
regr.mod.v = train(regr.lrn.v, regr.task.v, subset = train.index[,1])

regr.mod
regr.mod.v

task.pred = predict(regr.mod, task = regr.task, subset = train.index[,2], cores=10)
task.pred.v = predict(regr.mod.v, task = regr.task.v, subset = train.index[,2],cores=10)

performance(task.pred, measures=list(mse, rsq))
#mse         rsq 
0.006779552 0.606626456 

performance(task.pred.v, measures=list(mse, rsq))
#mse         rsq 
0.004794041 0.471933114 

plot(as.data.frame(task.pred)[,c("truth","response")])
plot(as.data.frame(task.pred.v)[,c("truth","response")])

listMeasures(regr.task)
listMeasures(regr.task.v)


##fine tune the model
getParamSet(regr.lrn)
getParamSet(regr.lrn.v)

library(parallel)
library(parallelMap)
 # Initialize paralelllization
parallelStartSocket(cpus = 18)

ps = makeParamSet(
makeIntegerParam("ntree", lower = 1, upper = 100),
makeIntegerParam("mtry", lower =  1, upper = 10 ))

ctrl = makeTuneControlGrid(resolution = 10)
rdesc = makeResampleDesc("CV", iters = 5)

tune.cforest = tuneParams(regr.lrn, task = regr.task, resampling = rdesc, par.set = ps, control = ctrl)
tune.cforest.v = tuneParams(regr.lrn.v, task = regr.task.v, resampling = rdesc, par.set = ps, control = ctrl)

plotHyperParsEffect(generateHyperParsEffectData(tune.cforest), x = "ntree", y = "mtry", z = "mse.test.mean",
plot.type = "heatmap")

tune.cforest$x
$ntree
[1] 100

$mtry
[1] 3

plotHyperParsEffect(generateHyperParsEffectData(tune.cforest.v), x = "ntree", y = "mtry", z = "mse.test.mean", plot.type = "heatmap")

tune.cforest.v$x
> tune.cforest.v$x
$ntree
[1] 78

$mtry
[1] 8

parallelStop()

save.image("rforestTune.all_gcb_rev1.RData")

regr.lrn.best = setHyperPars(makeLearner("regr.cforest"), ntree = tune.cforest$x$ntree, mtry = tune.cforest$x$mtry)
regr.lrn.best.v = setHyperPars(makeLearner("regr.cforest"), ntree = tune.cforest.v$x$ntree, mtry = tune.cforest.v$x$mtry)

regr.mod.best = train(regr.lrn.best, regr.task, subset = train.index[,1])
regr.mod.best.v = train(regr.lrn.best.v, regr.task.v, subset = train.index[,1])

save.image("rforestTune.all_gcb_rev1.RData")

vimp = unlist(getFeatureImportance(regr.mod.best)$res)
vimp.v = unlist(getFeatureImportance(regr.mod.best.v)$res)

#system.time(trial<-getFeatureImportance(regr.mod.best,conditional=TRUE,OOB=TRUE, cores=12))

#remotes::install_github("mlr3learners/mlr3learners.partykit")

barplot(vimp*100)
barplot(vimp.v)

#save.image('revision.anal.RData')
#load('revision.anal.RData')

##
pdp.all = generatePartialDependenceData(regr.mod.best, regr.task, individual = F, n = c(10, 25000))
pdp.all.vci = generatePartialDependenceData(regr.mod.best.v, regr.task.v, individual = F, n = c(10, 25000))
#plotPartialDependence(pdp.all.vci)

colnames(pdp.all$data)<-c('response',"elevation", "erosion","human pressure","land dev. intensity","sea level anomaly","slope","drought","sea level rate","mean sea level")      
pdp.ndvi<-melt(pdp.all$data,id.vars = "response")

pdp.ndvi$indicator<-"ndvi"

colnames(pdp.all.vci$data)<-c('response',"elevation", "erosion","human pressure","land dev. intensity","sea level anomaly","slope","tidal range","heat waves","drought","sea level rate","mean sea level")      
pdp.vci<-melt(pdp.all.vci$data,id.vars = "response")

pdp.vci$indicator<-"vci"

pdp.data<-rbind(pdp.ndvi,pdp.vci)


library(gridExtra)
library(ggthemes)
theme_ms <- function(base_size=12, base_family="Helvetica") {
  library(grid)
  (theme_bw(base_size = base_size, base_family = base_family)+
      theme(text=element_text(color="black"),
            axis.title=element_text(face="bold", size = rel(1.3)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text=element_text(size = rel(1), color = "black"),
            legend.title=element_text(face="bold"),
            legend.text=element_text(face="bold"),
            axis.ticks = element_line(),
            legend.background=element_rect(fill="transparent"),
            legend.key.size = unit(0.8, 'lines'),
            legend.position = "bottom",
            panel.border=element_rect(color="black",size=1),
            #panel.grid=element_blank(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
            
      ))
}

tiff("Figure3_mangroveMs.tiff", units="in", width=8, height=7, res=300)
fig<-ggplot(data=pdp.data,aes(x=value,y=response, col=indicator)) + geom_line() + facet_wrap(~variable, scales="free",ncol=3)+ xlab("Environmental driver") + ylab("Ecological condition index")
fig+theme_ms()
dev.off()

pdf("Figure3_mangroveMss.pdf") 
fig<-ggplot(data=pdp.data,aes(x=value,y=response, col=indicator)) + geom_line() + facet_wrap(~variable, scales="free",ncol=3)+ xlab("Environmental driver") + ylab("Ecological condition index")
fig+theme_ms()
dev.off()

##plot top variables
#make future predictions
#prepare data
#prediction data
predictors.future<-predictors
predictors.future$avmsl<-all.data$wcs.msl
predictors.future$slr<-all.data$slr
#predictors.future$slr<-all.data$wcs.trend
predictors.future$ccd<-all.data$cdd.rcp85
#predictors.future$tx90<-all.data$tx9085

summary(all.data[,c('avmsl','wcs.msl','slr','wcs.trend','ccd','cdd.rcp85','tx90','tx9085')])

##
future.ndvi<- predict(regr.mod.best, newdata = predictors.future, OOB=FALSE, type = "response", cores=24)
future.vci<- predict(regr.mod.best.v, newdata = predictors.future, OOB=FALSE, type = "response", cores=24)



##scale 
#If e < c < i, then the membership function is incr                   easing from e to i. If i < c < e, then the membership function is decreasing from i to e.
library(QCA)
#reinforcing<-data.frame(PredictVar[,c("elevation") ])
increasing<- function(x) {
  #minT<-quantile(x, probs = seq(0, 1, 0.01))[[2]]
  #maxT<-quantile(x, probs = seq(0, 1, 0.01))[[100]]
  minT<-min(x)
  maxT<-(max(x))
  mtxt<-paste0('calibrate(x, thresholds = "e=', minT,',c=',(maxT+minT)/2,',i=', maxT,'",logistic = FALSE)',  collapse="")
  reinf<-eval(parse(text = mtxt))
  return(reinf)
}

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

##scale
exposure<-data.frame(cbind(decreasing(future.ndvi$data$response),decreasing(future.vci$data$response)))
colnames(exposure)<-c("e.ndvi","e.vci")


##rasterrize
exposure.xy<-cbind(all.data[,2:3],exposure)
ndvi.dat<-raster("ndvi.nc")
ea.ext<-extent(ndvi.dat)
#nrow:13701
#ncol:9345
r <- raster(ea.ext, ncol=9345, nrow=13701)
ndvi.x <- rasterize(exposure.xy[, 1:2], r, exposure.xy[,3], fun=mean)
vci.x <- rasterize(exposure.xy[, 1:2], r, exposure.xy[,4], fun=mean)
exposure.stack<-stack(ndvi.x,vci.x)

writeRaster(exposure.stack, filename='exposure.stack.nc', format="CDF", overwrite=TRUE) 

rm(exposure.stack)

save.image("rforestTune.all_gcb_rev1.RData")#August5


##partial dependence and exposure
colnames(predictors.future)

#partial dependence
new.data1 <- within(predictors.future, elev.med <- ave(elevation, FUN = median) )
new.data2 <- within(new.data1, eros.med <- ave(erosion,  FUN = median) )
new.data3<-within(new.data2, grav.med <- ave(gravity,  FUN = median) )
new.data4 <- within(new.data3, ldi.med <- ave(ldi,  FUN = median) )
new.data5<-within(new.data4, sla.med <- ave(sla, FUN = median) )
new.data6 <- within(new.data5, slope.med <- ave(slope,  FUN = median) )
new.data7 <- within(new.data6, tide.med <- ave(tidecm,  FUN = median) )
new.data8<-within(new.data7, tx90.med <- ave(tx90,  FUN = median) )
new.data9 <- within(new.data8, ccd.med <- ave(ccd, FUN = median) )
new.data10 <- within(new.data9, slr.med <- ave(slr, FUN = median) )
new.data11 <- within(new.data10, avmsl.med <- ave(avmsl, FUN = median) )

colnames(predictors.future)
#[1] "elevation" "erosion"   "gravity"   "ldi"       "sla"       "slope"     "tidecm"   
#[8] "tx90"      "ccd"       "slr"       "avmsl"   

df1<- new.data11[,c("elevation","eros.med","grav.med","ldi.med","sla.med","slope.med","tide.med","tx90.med","ccd.med","slr.med", "avmsl.med","formation")]
colnames(df1)<-c("elevation","erosion","gravity","ldi","sla","slope","tidecm","tx90","ccd","slr","avmsl","formation" ) 

df2<- new.data11[,c("elev.med","erosion","grav.med","ldi.med","sla.med","slope.med","tide.med","tx90.med","ccd.med","slr.med", "avmsl.med","formation")]
colnames(df2)<-c("elevation","erosion","gravity","ldi","sla","slope","tidecm","tx90","ccd","slr","avmsl","formation" ) 

df3<- new.data11[,c("elev.med","eros.med","gravity","ldi.med","sla.med","slope.med","tide.med","tx90.med","ccd.med","slr.med", "avmsl.med","formation")]
colnames(df3)<-c("elevation","erosion","gravity","ldi","sla","slope","tidecm","tx90","ccd","slr","avmsl","formation" ) 

df4<- new.data11[,c("elev.med","eros.med","grav.med","ldi","sla.med","slope.med","tide.med","tx90.med","ccd.med","slr.med", "avmsl.med","formation")]
colnames(df4)<-c("elevation","erosion","gravity","ldi","sla","slope","tidecm","tx90","ccd","slr","avmsl","formation" ) 

df5<- new.data11[,c("elev.med","eros.med","grav.med","ldi.med","sla","slope.med","tide.med","tx90.med","ccd.med","slr.med", "avmsl.med","formation")]
colnames(df5)<-c("elevation","erosion","gravity","ldi","sla","slope","tidecm","tx90","ccd","slr","avmsl","formation" ) 

df6<- new.data11[,c("elev.med","eros.med","grav.med","ldi.med","sla.med","slope","tide.med","tx90.med","ccd.med","slr.med", "avmsl.med","formation")]
colnames(df6)<-c("elevation","erosion","gravity","ldi","sla","slope","tidecm","tx90","ccd","slr","avmsl","formation" ) 

df7<- new.data11[,c("elev.med","eros.med","grav.med","ldi.med","sla.med","slope.med","tidecm","tx90.med","ccd.med","slr.med", "avmsl.med","formation")]
colnames(df7)<-c("elevation","erosion","gravity","ldi","sla","slope","tidecm","tx90","ccd","slr","avmsl","formation" ) 

df8<- new.data11[,c("elev.med","eros.med","grav.med","ldi.med","sla.med","slope.med","tide.med","tx90","ccd.med","slr.med", "avmsl.med","formation")]
colnames(df8)<-c("elevation","erosion","gravity","ldi","sla","slope","tidecm","tx90","ccd","slr","avmsl","formation" ) 

df9<- new.data11[,c("elev.med","eros.med","grav.med","ldi.med","sla.med","slope.med","tide.med","tx90.med","ccd","slr.med", "avmsl.med","formation")]
colnames(df9)<-c("elevation","erosion","gravity","ldi","sla","slope","tidecm","tx90","ccd","slr","avmsl" ,"formation")

df10<- new.data11[,c("elev.med","eros.med","grav.med","ldi.med","sla.med","slope.med","tide.med","tx90.med","ccd.med","slr", "avmsl.med","formation")]
colnames(df10)<-c("elevation","erosion","gravity","ldi","sla","slope","tidecm","tx90","ccd","slr","avmsl","formation")

df11<- new.data11[,c("elev.med","eros.med","grav.med","ldi.med","sla.med","slope.med","tide.med","tx90.med","ccd.med","slr.med", "avmsl","formation")]
colnames(df11)<-c("elevation","erosion","gravity","ldi","sla","slope","tidecm","tx90","ccd","slr","avmsl","formation")

df12<- new.data11[,c("elev.med","eros.med","grav.med","ldi.med","sla.med","slope.med","tide.med","tx90.med","ccd.med","slr.med", "avmsl.med","formation")]
colnames(df12)<-c("elevation","erosion","gravity","ldi","sla","slope","tidecm","tx90","ccd","slr","avmsl","formation")

dfList<-list(df1,df2,df3,f4,df5,df6,df7,df8,df9,df10,df11,df12)

predict.ndvi<-mclapply(dfList, function(x) predict(regr.mod.best, newdata = x, OOB=FALSE, type = "response"))
ndvi.partial.exposure<- mclapply(predict.ndvi, function(x) decreasing(x$data$response))

predict.vci<-mclapply(dfList, function(x) predict(regr.mod.best.v, newdata = x, OOB=FALSE, type = "response"))
vci.partial.exposure<- mclapply(predict.vci, function(x) decreasing(x$data$response))


#save.image("rforestTune.all_gcb_rev1.RData")#aug62020


##convert lists to dataframe - ndvi
part.exp.df.ndvi <- as.data.frame(t(data.frame(matrix(unlist(ndvi.partial.exposure), nrow=length(ndvi.partial.exposure), byrow=T))))
colnames(part.exp.df.ndvi)<- c("elevation","erosion","gravity","ldi","sla","slope","tidecm","tx90","ccd","slr","avmsl","formation")
exposure.ndvi0<-cbind(all.data[,"layer"],part.exp.df.ndvi, exposure[,1])
colnames(exposure.ndvi0)[1]<-"sector"
colnames(exposure.ndvi0)[14]<-"e.ndvi"

drops <- c("formation","tidecm","tx90")
exposure.ndvi<-exposure.ndvi0[ , !(names(exposure.ndvi0) %in% drops)]

#all.exposure$exposure_index_fsum<-((all.exposure$e.ndvi+all.exposure$e.vci)-(all.exposure$e.vci*all.exposure$e.ndvi))
#all.exposure$exposure_index_ave<-(all.exposure$e.ndvi+all.exposure$e.vci)/2

##convert lists to dataframe - vci
part.exp.df.vci <- as.data.frame(t(data.frame(matrix(unlist(vci.partial.exposure), nrow=length(vci.partial.exposure), byrow=T))))
colnames(part.exp.df.vci)<- c("elevation","erosion","gravity","ldi","sla","slope","tidecm","tx90","ccd","slr","avmsl","formation")
exposure.vci0<-cbind(all.data[,"layer"],part.exp.df.vci, exposure[,2])
colnames(exposure.vci0)[1]<-"sector"
colnames(exposure.vci0)[14]<-"e.vci"
exposure.vci<-exposure.vci0[ , !(names(exposure.vci0) %in% drops)]

exposure.ndvi.1<-exposure.ndvi[,-1] #1
exposure.ndvi.1$slr <- 1-predict.ndvi[[10]]$data$response

exposure.vci.1<-exposure.vci[,-1] #2
exposure.vci.1$slr <- 1-predict.vci[[10]]$data$response

exposure_product_sum<-((exposure.vci.1+exposure.ndvi.1)-(exposure.vci.1*exposure.ndvi.1))#3
exposure_weighted_sum<-((0.47*exposure.vci.1)+(0.61*exposure.ndvi.1))/2


ndvi.vimp<-data.frame(vimp)*100
vci.vimp<-data.frame(vimp.v)*100
cumm.vimp<-data.frame(((vimp.v*100)+(vimp*100))-((vimp.v*100)*(vimp*100)))
colnames(cumm.vimp)[1]<-'vimp'
#ndvi.vimp1<-increasing(ndvi.vimp$vimp)
#vci.vimp1<-increasing(vci.vimp$vimp.v)

#wsum.vimp<-data.frame(((0.47*vci.vimp1)+(0.61*ndvi.vimp1))/2)
#colnames(wsum.vimp)<-"vimp"
##plot
#first create rownames vector for heatmap
country.mat = aggregate(exposure.ndvi.1,by = list(all.data$County,all.data$layer),FUN = median)
colnames(country.mat )[1]<-"Country"
tr<-unique(country.mat$Country)
val<-c('(Mz)','(Tz)','(Ken)','(Mdg)')
country.mat$country_1<-val[match(country.mat$Country,tr)]
rn<-paste(country.mat$country_1,country.mat$Group.2,sep="")

#overall
exposure.agg.median = aggregate(exposure_product_sum,by = list(all.data$layer),FUN = median)
#exposure.agg.median = aggregate(exposure_weighted_sum,by = list(all.data$layer),FUN = median)
rownames(exposure.agg.median)<-rn
exposure.agg.2<-exposure.agg.median[,-c(1)]
exposure.agg.2<-as.matrix(exposure.agg.2)
colnames(exposure.agg.2)<-c("elevation", "erosion","human pressure","land dev.intensity","sea level anomaly","slope","tidal range","heat waves","drought","sea level rate","mean sea level","formation","e.vci")      

#ndvi
exposure.median.ndvi = aggregate(exposure.ndvi.1,by = list(all.data$layer),FUN = median)
rownames(exposure.median.ndvi)<-rn
exposure.ndvi.2<-exposure.median.ndvi[,-c(1)]
exposure.ndvi.2<-as.matrix(exposure.ndvi.2)
colnames(exposure.ndvi.2)<-c("elevation", "erosion","human pressure","land dev.intensity","sea level anomaly","slope","tidal range","heat waves","drought","sea level rate","mean sea level","formation","e.ndvi")      

#vci
exposure.median.vci = aggregate(exposure.vci.1,by = list(all.data$layer),FUN = median)
rownames(exposure.median.vci)<-rn
exposure.vci.2<-exposure.median.vci[,-c(1)]
exposure.vci.2<-as.matrix(exposure.vci.2)
colnames(exposure.vci.2)<-c("elevation", "erosion","human pressure","land dev. intensity","sea level anomaly","slope","tidal range","heat waves","drought","sea level rate","mean sea level","formation","e.vci")      

#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(seriation)
library(dendextend)
col_fun = colorRamp2(c(0,0.25, 0.75, 1), c("blue","green", "yellow", "red"))
col_fun(seq(1, 10))
o = seriate(max(exposure.agg.2) - exposure.agg.2, method = "BEA_TSP")


##final plots
#Heatmap(exposure.agg.2[,c(1:11)], show_column_dend=FALSE,name = "mat", row_names_gp = gpar(fontsize = 9),border = FALSE,col = col_fun,row_gap = unit(2, "mm"),   row_km = 4, row_km_repeats = 1000, column_title_gp = gpar(font = 7),column_names_gp = gpar( fontsize = 11),cluster_column_slices = FALSE,heatmap_legend_param = list(col_fun = col_fun, title = "exposure index", legend_height = unit(8, "cm"),title_position = "lefttop-rot"),bottom_annotation = HeatmapAnnotation(vimp = anno_barplot(data.frame(ndvi.vimp$vimp*100)),annotation_name_side = "left"),right_annotation = rowAnnotation(Exposure = exposure.agg.2[,12],col = list(Exposure = col_fun),show_legend = c("Exposure" = FALSE)),left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:4),labels = c("group1", "group2", "group3", "group4"),  labels_gp = gpar(col = "white", fontsize = 10))))
#Final1
overall.plot<-Heatmap(exposure.agg.2[,c(1:11)], show_column_dend=FALSE,name = "mat",column_title = "(A) Cumulative exposure",column_title_gp = gpar(fontsize = 10,fontface="bold"), row_names_gp = gpar(fontsize = 8),border = FALSE,col = col_fun,row_gap = unit(2, "mm"),  row_km = 4, row_km_repeats = 1000, row_title=NULL,column_names_gp = gpar( fontsize = 9),cluster_column_slices = FALSE, heatmap_legend_param = list(col_fun = col_fun, title = "exposure index", direction="vertical",legend_height = unit(10, "cm"),title_position = "lefttop-rot"),bottom_annotation = HeatmapAnnotation(vimp = anno_barplot(data.frame(cumm.vimp$vimp*100)),annotation_name_side = "left",right_annotation = rowAnnotation(Exposure_overall = exposure.agg.2[,12],col = list(Exposure_overall = col_fun),annotation_name_gp=gpar(fontsize=8),show_legend = c("Exposure_overall" = FALSE)),left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = "grey"),labels = c("group1", "group2", "group3", "group4"),  labels_gp = gpar(col = "black", fontsize = 10))))
ndvi.plot<-Heatmap(exposure.ndvi.2[,c(1:11)], show_column_dend=FALSE,name = "mat",column_title = "(B) NDVI based exposure",column_title_gp = gpar(fontsize = 10,fontface="bold"), row_names_gp = gpar(fontsize = 8),border = FALSE,col = col_fun,row_gap = unit(2, "mm"),column_names_gp = gpar( fontsize = 9),cluster_column_slices = FALSE,heatmap_legend_param = list(col_fun = col_fun, title = "exposure index", legend_height = unit(8, "cm"),title_position = "lefttop-rot"),bottom_annotation = HeatmapAnnotation(vimp = anno_barplot(data.frame(ndvi.vimp$vimp*100)),show_annotation_name=c(bar=FALSE),annotation_name_side = "left"),right_annotation = rowAnnotation(Exposure_ndvi = exposure.ndvi.2[,12],col = list(Exposure_ndvi = col_fun),annotation_name_gp=gpar(fontsize=8),show_legend = c("Exposure_ndvi" = FALSE)))
vci.plot<-Heatmap(exposure.vci.2[,c(1:11)], show_column_dend=FALSE,name = "mat", column_title = "(C) VCI based exposure",column_title_gp = gpar(fontsize = 10,fontface="bold"),row_names_gp = gpar(fontsize =8),border = FALSE,col = col_fun,row_gap = unit(2, "mm"),column_names_gp = gpar( fontsize = 9),cluster_column_slices = FALSE,heatmap_legend_param = list(col_fun = col_fun, title = "exposure index", legend_height = unit(8, "cm"),title_position = "lefttop-rot"),bottom_annotation = HeatmapAnnotation(vimp = anno_barplot(data.frame(vci.vimp$vimp*100)),show_annotation_name=c(bar=FALSE),annotation_name_side = "left"),right_annotation = rowAnnotation(Exposure_vci = exposure.vci.2[,12],col = list(Exposure_vci = col_fun),annotation_name_gp=gpar(fontsize=8),show_legend = c("Exposure_vci" = FALSE)))

#Final2
overall.plot<-Heatmap(exposure.agg.2[,c(1:12)], show_column_dend=FALSE,name = "mat",column_title = "(A) Cumulative exposure",column_title_gp = gpar(fontsize = 8,fontface="bold"), row_names_gp = gpar(fontsize = 10),border = FALSE,col = col_fun,row_gap = unit(2, "mm"),  row_km = 4, row_km_repeats = 1000, row_title=NULL,column_names_gp = gpar( fontsize = 10),cluster_column_slices = FALSE, heatmap_legend_param = list(col_fun = col_fun, title = NULL, direction="horizontal",legend_width=unit(10, "cm"),title_position = "lefttop"),bottom_annotation = HeatmapAnnotation(vimp = anno_barplot(data.frame(cumm.vimp$vimp)),annotation_name_side = "left"),right_annotation = rowAnnotation(Exposure_overall = exposure.agg.2[,13],col = list(Exposure_overall = col_fun),annotation_name_gp=gpar(fontsize=10),show_legend = c("Exposure_overall" = FALSE)),left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = "grey"),labels = c("group1", "group2", "group3", "group4"),  labels_gp = gpar(col = "black", fontsize = 10))))
ndvi.plot<-Heatmap(exposure.ndvi.2[,c(1:12)], show_column_dend=FALSE,name = "mat",column_title = "(B) NDVI based exposure",column_title_gp = gpar(fontsize = 8,fontface="bold"), row_names_gp = gpar(fontsize = 10),border = FALSE,col = col_fun,row_gap = unit(2, "mm"),column_names_gp = gpar( fontsize = 10),cluster_column_slices = FALSE,heatmap_legend_param = list(col_fun = col_fun, title = NULL, direction="horizontal",legend_width=unit(10, "cm"),title_position = "lefttop"),bottom_annotation = HeatmapAnnotation(vimp = anno_barplot(data.frame(ndvi.vimp$vimp)),show_annotation_name=c(bar=FALSE),annotation_name_side = "left"),right_annotation = rowAnnotation(Exposure_ndvi = exposure.ndvi.2[,13],col = list(Exposure_ndvi = col_fun),annotation_name_gp=gpar(fontsize=10),show_legend = c("Exposure_ndvi" = FALSE)))
vci.plot<-Heatmap(exposure.vci.2[,c(1:12)], show_column_dend=FALSE,name = "mat",column_title = "(C) VCI based exposure",column_title_gp = gpar(fontsize = 8,fontface="bold"),row_names_gp = gpar(fontsize =10),border = FALSE,col = col_fun,row_gap = unit(2, "mm"),column_names_gp = gpar( fontsize = 10),cluster_column_slices = FALSE,heatmap_legend_param = list(col_fun = col_fun, title = NULL, direction="horizontal",legend_width=unit(10, "cm"),title_position = "lefttop"),bottom_annotation = HeatmapAnnotation(vimp = anno_barplot(data.frame(vci.vimp$vimp)),show_annotation_name=c(bar=FALSE),annotation_name_side = "left"),right_annotation = rowAnnotation(Exposure_vci = exposure.vci.2[,13],col = list(Exposure_vci = col_fun),annotation_name_gp=gpar(fontsize=10),show_legend = c("Exposure_vci" = FALSE)))


pdf("figure2_mangroveMs_revision.pdf") 
ht_list = overall.plot+ndvi.plot + vci.plot
draw(ht_list,heatmap_legend_side="bottom", ht_gap = unit(5, "mm"))
dev.off()

tiff("Figure2_mangroveMs_revision.tiff", units="in", width=10, height=8, res=300)
ht_list = overall.plot+ndvi.plot + vci.plot
draw(ht_list,heatmap_legend_side="bottom", ht_gap = unit(5, "mm"))
dev.off()





save.image("rforestTune.all.RData")
load("rforestTune.all.RData")

#parallelStartSocket(cpus = 12)
#pdp.vci = generatePartialDependenceData(regr.mod.best.v, regr.task.v, individual = F)
## Loading required package: mmpf
#plotPartialDependence(pdp.all)
      
#####Run above for VCI
num_cores=18
#vci.crf.1<-partykit::cforest(vci ~elevation+erosion+gravity+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl+formation,data=all.data.vci,cores=num_cores,ntree=89, mtry=3)
#ndvi.crf.1<-partykit::cforest(vci ~elevation+erosion+gravity+ldi+sla +slope + tidecm+tx90+ccd+slr+avmsl+formation,data=all.data.vci,cores=num_cores,ntree=89, mtry=3)

train1 <- all.data.ndvi[train.index[,1],]
ndvi.crf.1<-partykit::cforest(ndvi ~elevation+erosion+gravity+ldi+sla +slope +ccd+slr+avmsl,data=train1,cores=num_cores,ntree=100, mtry=3)
vci.crf.1<-partykit::cforest(vci ~elevation+erosion+gravity+ldi+sla +slope + ccd+slr+avmsl,data=all.data.vci,cores=num_cores,ntree=89, mtry=3)

#control = cforest_unbiased(ntree = 50))

system.time(varimpvci <- varimp(ndvi.crf.1,conditional=TRUE, ncores=18))



library(doParallel)
cl<-makeCluster(12)
registerDoParallel(cl)
pvci.1 <- partial(vci.crf.1, pred.va = c("elevation", "erosion","gravity"), chull = TRUE, parallel = TRUE)
pvci.2 <- partial(vci.crf.1, pred.var = c("ldi", "sla","slope"), chull = TRUE, parallel = TRUE)
pvci.3 <- partial(vci.crf.1, pred.var = c("tidecm", "tx90","ccd"), chull = TRUE, parallel = TRUE)
pvci.4 <- partial(vci.crf.1, pred.var = c("slr", "avmsl"), chull = TRUE, parallel = TRUE)
stopCluster()











