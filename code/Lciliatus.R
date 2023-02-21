library(ENMeval)
library(raster)
library(dplyr)
library(spocc)
library(maptools)
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

Lbrev.occs <- select(RAW_Lbrev, 2:3)
Lbrev.occs.cells <- raster::extract(envs, Lbrev.occs, cellnumbers = TRUE)
Lbrev.occs.cellDups <- duplicated(Lbrev.occs.cells[,1])
Lbrev.occs <- Lbrev.occs[!Lbrev.occs.cellDups,]
Lbrev.occs.sp <- sp::SpatialPoints(Lbrev.occs)

Lbrev.occs.sf <- sf::st_as_sf(Lbrev.occs, coords=c("longitude","latitude")
                              , crs=raster::crs(envs))

eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
Lbrev.occs.sf <- sf::st_transform(Lbrev.occs.sf, crs = eckertIV)

Lbrev.occs.buf <- sf::st_buffer(Lbrev.occs.sf, dist = 100000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(envs))
plot(envs[[1]], main = names(envs)[1])
points(Lbrev.occs)
plot(Lbrev.occs.buf, border = "blue", lwd = 3, add = TRUE)
envs.bg <- raster::crop(envs, Lbrev.occs.buf)
envs.bg <- raster::mask(envs.bg, Lbrev.occs.buf)

plot(envs.bg, main = names(envs)[1])
points(Lbrev.occs)
plot(Lbrev.occs.buf, border = "blue", lwd = 3, add = TRUE)


corr <- layerStats(envs.bg, 'pearson', na.rm=T)
c <- corr$'pearson correlation coefficient'
write.csv(c, file='../Lbrevbioclimcorr.csv')
c <- data.matrix(read.csv("../Lbrevbioclimcorr.csv", header=T, row.names = 1, 
                          sep=","))
library(caret)
envtCor <- findCorrelation (c, cutoff=0.70, names=T, exact=T)
sort(envtCor)

envt.subset <- stack(envs$bio6, envs$bio7,
                     envs$bio8, envs$bio18)
envs.files <- c(envt.subset)
envs <- raster::stack(envs.files)
rasterVis::levelplot(envs[[-9]], main="Leptosiphon aureus Occurence Points", margin=F) +
  latticeExtra::layer(sp.points(Lbrev.occs.sp, col="pink"))

occs.z <- raster::extract(envs, Lbrev.occs)
Lbrev.occs.sim <- similarity(envs.bg, occs.z)
Lbrev.occs.mess <- Lbrev.occs.sim$similarity_min
Lbrev.occs.sp <- sp::SpatialPoints(Lbrev.occs)  
myScale <- seq(cellStats(Lbrev.occs.mess, min), cellStats(Lbrev.occs.mess, max), length.out = 100)
rasterVis::levelplot(Lbrev.occs.mess, main = "Environmental similarity", at = myScale, margin = FALSE) + 
  latticeExtra::layer(sp.points(Lbrev.occs.sp, col="black"))

cols <- RColorBrewer::brewer.pal(8, "Set1")
rasterVis::levelplot(Lbrev.occs.sim$mod, col.regions = cols, main = "Most different variable") + 
  latticeExtra::layer(sp.points(Lbrev.occs.sp, col="black"))
rasterVis::levelplot(Lbrev.occs.sim$mos, col.regions = cols, main = "Most similar variable") + 
  latticeExtra::layer(sp.points(Lbrev.occs.sp, col="black"))
Lbrev.occs.sf <- sf::st_as_sf(Lbrev.occs, coords=c("longitude","latitude")
                              , crs=raster::crs(envs))
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
Lbrev.occs.sf <- sf::st_transform(Lbrev.occs.sf, crs = eckertIV)
Lbrev.occs.buf <- sf::st_buffer(Lbrev.occs.sf, dist = 100000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(envs))
plot(envs[[1]], main = names(envs)[1])
points(Lbrev.occs)
plot(Lbrev.occs.buf, border = "blue", lwd = 3, add = TRUE)
envs.bg <- raster::crop(envs, Lbrev.occs.buf)
envs.bg <- raster::mask(envs.bg, Lbrev.occs.buf)
plot(envs.bg, main = names(envs)[1])
points(Lbrev.occs)
plot(Lbrev.occs.buf, border = "blue", lwd = 3, add = TRUE)
bg <- dismo::randomPoints(envs.bg[[1]], n = 10000) %>% as.data.frame()
colnames(bg) <- colnames(Lbrev.occs)
plot(envs.bg)
plot(envs.bg)
png("LbrevTrainingArea.png")
plot(envs.bg)
dev.off()
occs.z <- cbind(Lbrev.occs, raster::extract(envs, Lbrev.occs))
bg.z <- cbind(bg, raster::extract(envs, bg))
grp.n <- 6
kmeans <- kmeans(Lbrev.occs, grp.n)
occs.grp <- kmeans$cluster
evalplot.grps(pts = Lbrev.occs, pts.grp = occs.grp, envs = envs.bg)
bg.grp <- rep(0, nrow(bg))
evalplot.grps(pts = bg, pts.grp = bg.grp, envs = envs.bg)
centers <- kmeans$center
d <- raster::pointDistance(bg, centers, lonlat = TRUE)
bg.grp <- apply(d, 1, function(x) which(x == min(x)))
png("LbrevUserBGPartitions")
evalplot.grps(pts = bg, pts.grp = bg.grp, envs = envs.bg)
dev.off()
dev.new()
occsBg.sf <- sf::st_as_sf(rbind(Lbrev.occs, bg), coords = c("longitude","latitude"), 
                          crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
raster::crs(envs.bg) <- raster::crs(occsBg.sf)
sb <- blockCV::spatialBlock(speciesData = occsBg.sf, rasterLayer = envs.bg, 
                            theRange = 250000, k = 3, selection = "random")
foldExplorer(sb, envs.bg, Lbrev.occs.sp)
occs.grp <- sb$foldID[1:nrow(Lbrev.occs)]
bg.grp <- sb$foldID[(nrow(Lbrev.occs)+1):length(sb$foldID)]
png("Lbrev User Partitions.png")
evalplot.grps(pts = bg, pts.grp = bg.grp, envs = envs.bg)
dev.off()
dev.new()
dev.capture()
tune.args <- list(fc = c("L","H","Q","LH", "LQ", "LQH"), rm = 1:5)
user.grp <- list(occs.grp = round(runif(nrow(Lbrev.occs), 1, 2)), 
                 bg.grp = round(runif(nrow(bg), 1, 2)))
Lbrev.e.user <- ENMevaluate(Lbrev.occs, envs, bg, algorithm = "maxnet", overlap = T,
                            tune.args = tune.args, partitions = "user", user.grp = user.grp)

Lbrev.e.user@results
Lbrev.e.user@results.partitions
Lbrev.e.user
str(Lbrev.e.user, max.level=2)
eval.algorithm(Lbrev.e.user)
eval.tune.settings(Lbrev.e.user) %>% head()
eval.partition.method(Lbrev.e.user)
eval.results(Lbrev.e.user) %>% head()
eval.results.partitions(Lbrev.e.user) %>% head()
eval.models(Lbrev.e.user) %>% str(max.level = 1)
m1.mx <- eval.models(Lbrev.e.user)[["fc.LQH_rm.1"]]
m1.mx$betas
enframe(m1.mx$betas)
eval.predictions(Lbrev.e.user)
eval.occs(Lbrev.e.user) %>% head()
eval.bg(Lbrev.e.user) %>% head()
eval.occs.grp(Lbrev.e.user) %>% str()
eval.bg.grp(Lbrev.e.user) %>% str()
evalplot.stats(e = Lbrev.e.user, stats = "or.mtp", color = "fc", x.var = "rm")
evalplot.stats(e = Lbrev.e.user, stats = c("or.mtp", "auc.val"), color = "fc", x.var = "rm")
evalplot.stats(e = Lbrev.e.user, stats = c("or.mtp", "auc.val"), color = "fc", x.var = "rm", 
               error.bars = FALSE)
evalplot.stats(e = Lbrev.e.user, stats = c("or.mtp", "auc.val"), color = "fc", x.var = "rm", 
               dodge = 0.5)
evalplot.stats(e = Lbrev.e.user, stats = c("or.mtp", "auc.val"), color = "rm", x.var = "fc", 
               error.bars = FALSE, dodge = 0.6)
res <- eval.results(Lbrev.e.user)
opt.aicc <- res %>% filter(delta.AICc == 0)
opt.aicc
opttrain <- eval.models(Lbrev.e.user)[[opt.aicc$tune.args]]
dev.new()
png("Lbrev Optimal Performance Cloglog Plots.png")
plot(opttrain, type="cloglog", main= "L. breviculus Occurence Probability Model: Optimal Performance Metrics")
dev.off()
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq
mod.seq <- eval.models(Lbrev.e.user)[[opt.seq$tune.args]]
mod.seq$betas
png("Lbrev PV Occ Prob cloglog.png")
plot(mod.seq, type = "cloglog", main="L. breviculus Occurence Probability by PV")
dev.off()
dev.off()
pred.seq <- eval.predictions(Lbrev.e.user)[[opt.seq$tune.args]]
png("Lbrev sequential qualifier cloglog.ong")
plot(pred.seq)
dev.off()
points(eval.bg(Lbrev.e.user), pch = 3, col = eval.bg.grp(Lbrev.e.user), cex = 0.5)
points(eval.occs(Lbrev.e.user), pch = 21, bg = eval.occs.grp(Lbrev.e.user))
mod.simple <- eval.models(Lbrev.e.user)[['fc.L_rm.5']]
mod.complex <- eval.models(Lbrev.e.user)[['fc.LQH_rm.1']]
mod.simple$betas
length(mod.simple$betas)
mod.complex$betas
length(mod.complex$betas)
dev.n
png("Lbrev Simplest Model cloglog.png")
plot(mod.simple, type = "cloglog")
dev.off()
png("Lbrev Most Complex cloglog.png")
plot(mod.complex, type = "cloglog")
dev.off()
dev.new()
par(mfrow=c(2,1), mar=c(2,1,2,0))
plot(eval.predictions(Lbrev.e.user)[['fc.L_rm.5']], ylim = c(30,50), xlim = c(-90,-30), 
     legend = FALSE, main = 'Simplest Model Prediction')
plot(eval.predictions(Lbrev.e.user)[['fc.LQH_rm.1']],  ylim = c(30,50), xlim = c(-90,-30), 
     legend = FALSE, main = 'Most Complex Model Prediction')
png("L. breviculus H_RM4 (deltaAICc = 0) Prediction.png")
plot(eval.predictions(Lbrev.e.user)[['fc.H_rm.4']], ylim = c(30,50), xlim = c(-90,-30), 
     legend = FALSE, main = 'L. breviculus H_RM4 (deltaAICc = 0) Prediction')
dev.off()
mod.null <- ENMnulls(Lbrev.e.user , user.enm = Lbrev.e.user@algorithm, user.eval.type = "kspatial", mod.settings = list(fc = c("LQ"), rm = 3), no.iter = 100)
null.results(mod.null) %>% head()
null.results.partitions(mod.null) %>% head()
nullres <- null.results.partitions(mod.null)
null.emp.results(mod.null)
evalplot.nulls(mod.null, stats = c("or.10p", "auc.val"), plot.type = "histogram")
evalplot.nulls(mod.null, stats = c("or.10p", "auc.val"), plot.type = "violin")
rmm <- eval.rmm(Lbrev.e.user)
rmm$model$selectionRules <- "lowest 10 percentile omission rate, 
break ties with average validation AUC"
rmm$model$finalModelSettings <- "L1"
rmm$prediction$continuous$minVal <- cellStats(pred.seq, min)
rmm$prediction$continuous$maxVal <- cellStats(pred.seq, max)
rmm$prediction$continuous$units <- "suitability (cloglog transformation)"
rangeModelMetadata::rmmToCSV(rmm, "rmm_mx1.csv")

