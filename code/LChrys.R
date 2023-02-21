library(ENMeval)
library(raster)
library(dplyr)
library(spocc)
library(maptools)
install.packages(raster)
library(raster)
library(rgdal)
library(sp)
library(maps)
library(mapproj)
library(ggplot2)
library(rasterVis)
library(blockCV)
library(sf)
library(devtools)
library(tibble)
library(ecospat)

alt_11 <- raster::getData('worldclim', var='alt', res=2.5, lon=-119, lat=30)
alt_11 <- crop(alt_11, bio1)
envs <- stack(bio_01)
envs <- crop(envs, bio1)
#envs <- mask(envs, bio1)
envs <- stack(envs, alt_11)
plot(envs)
memory.limit(size=10000000000)

Lchryschrys.occs <- select(RAW_Lchryschrys, 2:3)
Lchryschrys.occs.cells <- raster::extract(envs, Lchryschrys.occs, cellnumbers = TRUE)
Lchryschrys.occs.cellDups <- duplicated(Lchryschrys.occs.cells[,1])
Lchryschrys.occs <- Lchryschrys.occs[!Lchryschrys.occs.cellDups,]
Lchryschrys.occs.sp <- sp::SpatialPoints(Lchryschrys.occs)

Lchryschrys.occs.sf <- sf::st_as_sf(Lchryschrys.occs, coords=c("longitude","latitude")
                              , crs=raster::crs(envs))

eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
Lchryschrys.occs.sf <- sf::st_transform(Lchryschrys.occs.sf, crs = eckertIV)

Lchryschrys.occs.buf <- sf::st_buffer(Lchryschrys.occs.sf, dist = 100000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(envs))
plot(envs[[1]], main = names(envs)[1])
points(Lchryschrys.occs)
plot(Lchryschrys.occs.buf, border = "blue", lwd = 3, add = TRUE)
envs.bg <- raster::crop(envs, Lchryschrys.occs.buf)
envs.bg <- raster::mask(envs.bg, Lchryschrys.occs.buf)

plot(envs.bg, main = names(envs)[1])
points(Lchryschrys.occs)
plot(Lchryschrys.occs.buf, border = "blue", lwd = 3, add = TRUE)


corr <- layerStats(envs.bg, 'pearson', na.rm=T)
c <- corr$'pearson correlation coefficient'
write.csv(c, file='../Lchryschrysbioclimcorr.csv')
c <- data.matrix(read.csv("../Lchryschrysbioclimcorr.csv", header=T, row.names = 1, 
                          sep=","))
library(caret)
envtCor <- findCorrelation (c, cutoff=0.70, names=T, exact=T)
sort(envtCor)

envt.subset <- stack(envs$bio6, envs$bio7,
                     envs$bio8, envs$bio18)
envs.files <- c(envt.subset)
envs <- raster::stack(envs.files)
rasterVis::levelplot(envs[[-9]], main="Leptosiphon aureus Occurence Points", margin=F) +
  latticeExtra::layer(sp.points(Lchryschrys.occs.sp, col="pink"))

occs.z <- raster::extract(envs, Lchryschrys.occs)
Lchryschrys.occs.sim <- similarity(envs.bg, occs.z)
Lchryschrys.occs.mess <- Lchryschrys.occs.sim$similarity_min
Lchryschrys.occs.sp <- sp::SpatialPoints(Lchryschrys.occs)  
myScale <- seq(cellStats(Lchryschrys.occs.mess, min), cellStats(Lchryschrys.occs.mess, max), length.out = 100)
rasterVis::levelplot(Lchryschrys.occs.mess, main = "Environmental similarity", at = myScale, margin = FALSE) + 
  latticeExtra::layer(sp.points(Lchryschrys.occs.sp, col="black"))

cols <- RColorBrewer::brewer.pal(8, "Set1")
rasterVis::levelplot(Lchryschrys.occs.sim$mod, col.regions = cols, main = "Most different variable") + 
  latticeExtra::layer(sp.points(Lchryschrys.occs.sp, col="black"))
rasterVis::levelplot(Lchryschrys.occs.sim$mos, col.regions = cols, main = "Most similar variable") + 
  latticeExtra::layer(sp.points(Lchryschrys.occs.sp, col="black"))
Lchryschrys.occs.sf <- sf::st_as_sf(Lchryschrys.occs, coords=c("longitude","latitude")
                              , crs=raster::crs(envs))
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
Lchryschrys.occs.sf <- sf::st_transform(Lchryschrys.occs.sf, crs = eckertIV)
Lchryschrys.occs.buf <- sf::st_buffer(Lchryschrys.occs.sf, dist = 100000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(envs))
plot(envs[[1]], main = names(envs)[1])
points(Lchryschrys.occs)
plot(Lchryschrys.occs.buf, border = "blue", lwd = 3, add = TRUE)
envs.bg <- raster::crop(envs, Lchryschrys.occs.buf)
envs.bg <- raster::mask(envs.bg, Lchryschrys.occs.buf)
plot(envs.bg, main = names(envs)[1])
points(Lchryschrys.occs)
plot(Lchryschrys.occs.buf, border = "blue", lwd = 3, add = TRUE)
bg <- dismo::randomPoints(envs.bg[[1]], n = 10000) %>% as.data.frame()
colnames(bg) <- colnames(Lchryschrys.occs)
plot(envs.bg)
plot(envs.bg)
png("LchryschrysTrainingArea.png")
plot(envs.bg)
dev.off()
occs.z <- cbind(Lchryschrys.occs, raster::extract(envs, Lchryschrys.occs))
bg.z <- cbind(bg, raster::extract(envs, bg))
grp.n <- 6
kmeans <- kmeans(Lchryschrys.occs, grp.n)
occs.grp <- kmeans$cluster
evalplot.grps(pts = Lchryschrys.occs, pts.grp = occs.grp, envs = envs.bg)
bg.grp <- rep(0, nrow(bg))
evalplot.grps(pts = bg, pts.grp = bg.grp, envs = envs.bg)
centers <- kmeans$center
d <- raster::pointDistance(bg, centers, lonlat = TRUE)
bg.grp <- apply(d, 1, function(x) which(x == min(x)))
png("LchryschrysUserBGPartitions")
evalplot.grps(pts = bg, pts.grp = bg.grp, envs = envs.bg)
dev.off()
dev.new()
occsBg.sf <- sf::st_as_sf(rbind(Lchryschrys.occs, bg), coords = c("longitude","latitude"), 
                          crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
raster::crs(envs.bg) <- raster::crs(occsBg.sf)
sb <- blockCV::spatialBlock(speciesData = occsBg.sf, rasterLayer = envs.bg, 
                            theRange = 250000, k = 3, selection = "random")
foldExplorer(sb, envs.bg, Lchryschrys.occs.sp)
occs.grp <- sb$foldID[1:nrow(Lchryschrys.occs)]
bg.grp <- sb$foldID[(nrow(Lchryschrys.occs)+1):length(sb$foldID)]
png("Lchryschrys User Partitions.png")
evalplot.grps(pts = bg, pts.grp = bg.grp, envs = envs.bg)
dev.off()
dev.new()
dev.capture()
tune.args <- list(fc = c("L","H","Q","LH", "LQ", "LQH"), rm = 1:5)
user.grp <- list(occs.grp = round(runif(nrow(Lchryschrys.occs), 1, 2)), 
                 bg.grp = round(runif(nrow(bg), 1, 2)))
Lchryschrys.e.user <- ENMevaluate(Lchryschrys.occs, envs, bg, algorithm = "maxnet", overlap = T,
                            tune.args = tune.args, partitions = "user", user.grp = user.grp)

Lchryschrys.e.user@results
Lchryschrys.e.user@results.partitions
Lchryschrys.e.user
str(Lchryschrys.e.user, max.level=2)
eval.algorithm(Lchryschrys.e.user)
eval.tune.settings(Lchryschrys.e.user) %>% head()
eval.partition.method(Lchryschrys.e.user)
eval.results(Lchryschrys.e.user) %>% head()
eval.results.partitions(Lchryschrys.e.user) %>% head()
eval.models(Lchryschrys.e.user) %>% str(max.level = 1)
m1.mx <- eval.models(Lchryschrys.e.user)[["fc.LQH_rm.1"]]
m1.mx$betas
enframe(m1.mx$betas)
eval.predictions(Lchryschrys.e.user)
eval.occs(Lchryschrys.e.user) %>% head()
eval.bg(Lchryschrys.e.user) %>% head()
eval.occs.grp(Lchryschrys.e.user) %>% str()
eval.bg.grp(Lchryschrys.e.user) %>% str()
evalplot.stats(e = Lchryschrys.e.user, stats = "or.mtp", color = "fc", x.var = "rm")
evalplot.stats(e = Lchryschrys.e.user, stats = c("or.mtp", "auc.val"), color = "fc", x.var = "rm")
evalplot.stats(e = Lchryschrys.e.user, stats = c("or.mtp", "auc.val"), color = "fc", x.var = "rm", 
               error.bars = FALSE)
evalplot.stats(e = Lchryschrys.e.user, stats = c("or.mtp", "auc.val"), color = "fc", x.var = "rm", 
               dodge = 0.5)
evalplot.stats(e = Lchryschrys.e.user, stats = c("or.mtp", "auc.val"), color = "rm", x.var = "fc", 
               error.bars = FALSE, dodge = 0.6)
res <- eval.results(Lchryschrys.e.user)
opt.aicc <- res %>% filter(delta.AICc == 0)
opt.aicc
opttrain <- eval.models(Lchryschrys.e.user)[[opt.aicc$tune.args]]
dev.new()
png("Lchryschrys Optimal Performance Cloglog Plots.png")
plot(opttrain, type="cloglog", main= "L. chryschrysiculus Occurence Probability Model: Optimal Performance Metrics")
dev.off()
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq
mod.seq <- eval.models(Lchryschrys.e.user)[[opt.seq$tune.args]]
mod.seq$betas
png("Lchryschrys PV Occ Prob cloglog.png")
plot(mod.seq, type = "cloglog", main="L. chryschrysiculus Occurence Probability by PV")
dev.off()
dev.off()
pred.seq <- eval.predictions(Lchryschrys.e.user)[[opt.seq$tune.args]]
png("Lchryschrys sequential qualifier cloglog.ong")
plot(pred.seq)
dev.off()
points(eval.bg(Lchryschrys.e.user), pch = 3, col = eval.bg.grp(Lchryschrys.e.user), cex = 0.5)
points(eval.occs(Lchryschrys.e.user), pch = 21, bg = eval.occs.grp(Lchryschrys.e.user))
mod.simple <- eval.models(Lchryschrys.e.user)[['fc.L_rm.5']]
mod.complex <- eval.models(Lchryschrys.e.user)[['fc.LQH_rm.1']]
mod.simple$betas
length(mod.simple$betas)
mod.complex$betas
length(mod.complex$betas)
dev.n
png("Lchryschrys Simplest Model cloglog.png")
plot(mod.simple, type = "cloglog")
dev.off()
png("Lchryschrys Most Complex cloglog.png")
plot(mod.complex, type = "cloglog")
dev.off()
dev.new()
par(mfrow=c(2,1), mar=c(2,1,2,0))
plot(eval.predictions(Lchryschrys.e.user)[['fc.L_rm.5']], ylim = c(30,50), xlim = c(-90,-30), 
     legend = FALSE, main = 'Simplest Model Prediction')
plot(eval.predictions(Lchryschrys.e.user)[['fc.LQH_rm.1']],  ylim = c(30,50), xlim = c(-90,-30), 
     legend = FALSE, main = 'Most Complex Model Prediction')
png("L. chryschrysiculus H_RM4 (deltaAICc = 0) Prediction.png")
plot(eval.predictions(Lchryschrys.e.user)[['fc.H_rm.4']], ylim = c(30,50), xlim = c(-90,-30), 
     legend = FALSE, main = 'L. chryschrysiculus H_RM4 (deltaAICc = 0) Prediction')
dev.off()
mod.null <- ENMnulls(Lchryschrys.e.user , user.enm = Lchryschrys.e.user@algorithm, user.eval.type = "kspatial", mod.settings = list(fc = c("LQ"), rm = 3), no.iter = 100)
null.results(mod.null) %>% head()
null.results.partitions(mod.null) %>% head()
nullres <- null.results.partitions(mod.null)
null.emp.results(mod.null)
evalplot.nulls(mod.null, stats = c("or.10p", "auc.val"), plot.type = "histogram")
evalplot.nulls(mod.null, stats = c("or.10p", "auc.val"), plot.type = "violin")
rmm <- eval.rmm(Lchryschrys.e.user)
rmm$model$selectionRules <- "lowest 10 percentile omission rate, 
break ties with average validation AUC"
rmm$model$finalModelSettings <- "L1"
rmm$prediction$continuous$minVal <- cellStats(pred.seq, min)
rmm$prediction$continuous$maxVal <- cellStats(pred.seq, max)
rmm$prediction$continuous$units <- "suitability (cloglog transformation)"
rangeModelMetadata::rmmToCSV(rmm, "rmm_mx1.csv")

