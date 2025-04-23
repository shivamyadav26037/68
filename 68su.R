#1supervised ###############
#1
library(terra)
raslist <- paste0('data/rs/LC08_044034_20170614_B', 2:7, ".tif")
landsat <- rast(raslist)
names(landsat) <- c('blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2')
samp <- readRDS("data/rs/lcsamples.rds")
plot(samp)
text(samp, samp$class)
#2
set.seed(1)
ptsamp <- spatSample(samp, 200, method="random")
plot(ptsamp, "class")
#3
nlcd <- rast('data/rs/nlcd-L1.tif')
names(nlcd) <- c("nlcd2001", "nlcd2011")
nlcd2011 <- nlcd[[2]]
nlcdclass <- c("Water", "Developed", "Barren", "Forest", "Shrubland",
               "Herbaceous", "Cultivated", "Wetlands")
classdf <- data.frame(value = c(1,2,3,4,5,7,8,9), names = nlcdclass)
levels(nlcd2011) <- classdf
classcolor <- c("#5475A8", "#B50000", "#D2CDC0", "#38814E",
                "#AF963C", "#D1D182", "#FBF65D", "#C8E6F8")
plot(nlcd2011, col=classcolor)
ptlonlat <- project(ptsamp, crs(nlcd2011))
points(ptlonlat)
#4
samp2011 <- spatSample(nlcd2011, size = 200, method="regular")
table(samp2011[,1])
df <- extract(landsat, ptsamp, ID=FALSE)
head(df)
sampdata <- data.frame(class = ptsamp$class, df)

library(rpart)
cartmodel <- rpart(as.factor(class)~., data = sampdata, method =
                     'class', minsplit = 5)
print(cartmodel)
plot(cartmodel, uniform=TRUE, main="Classification Tree")
text(cartmodel, cex = 1)
#5
classified <- predict(landsat, cartmodel, na.rm = TRUE)
plot(classified)

#6
lulc <- which.max(classified)
lulc
cls <- c("built","cropland","fallow","open","water")
df <- data.frame(id = 1:5, class=cls)
levels(lulc) <- df
lulc
mycolor <- c("darkred", "yellow", "burlywood", "cyan", "blue")
plot(lulc, col=mycolor)

set.seed(99)
k <- 5
j <- sample(rep(1:k, each = round(nrow(sampdata))/k))
table(j)
x <- list()
for (k in 1:5) {
  train <- sampdata[j!= k, ]
  
  test <- sampdata[j == k, ]
  cart <- rpart(as.factor(class)~., data=train, method = 'class',
                minsplit = 5)
  pclass <- predict(cart, test, na.rm = TRUE)
  pc <- apply(pclass, 1, which.max)
  pc <- colnames(pclass)[pc]
  x[[k]] <- cbind(test$class, pc)
}
y <- do.call(rbind, x)
y <- data.frame(y)
colnames(y) <- c('observed', 'predicted')
conmat <- table(y)
print(conmat)
n <- sum(conmat)
n
diag <- diag(conmat)
OA <- sum(diag) / n
OA
rowsums <- apply(conmat, 1, sum)
p <- rowsums / n
colsums <- apply(conmat, 2, sum)
q <- colsums / n
expAccuracy <- sum(p*q)
kappa <- (OA - expAccuracy) / (1 - expAccuracy)
kappa
PA <- diag / colsums
UA <- diag / rowsums
outAcc <- data.frame(producerAccuracy = PA, userAccuracy = UA)
outAcc




#2unsupervised ###########################
library(terra)
landsat5 <- rast('data/rs/centralvalley-2011LT5.tif')
names(landsat5) <- c('blue','green','red','NIR','SWIR1','SWIR2')

ndvi <- (landsat5[['NIR']] - landsat5[['red']]) / (landsat5[['NIR']] + landsat5[['red']])


e <- ext(-121.807, -121.725, 38.004, 38.072)


ndvi <- crop(ndvi, e)
ndvi

nr <- as.data.frame(ndvi, cell=TRUE)
str(nr)

set.seed(99)

kmncluster <- kmeans(nr[,-1], centers=10, iter.max = 500, nstart = 5, algorithm="Lloyd")

str(kmncluster)

knr <- rast(ndvi, nlyr=1)
knr[nr$cell] <- kmncluster$cluster
knr
mycolor <- c("#fef65b","#ff0000", "#daa520","#0000ff","#0000ff","#00ff00","#cbbeb5",
             "#c3ff5b", "#ff7373", "#00ff00", "#808080")
dev.new(width = 12, height = 6)  # Open a larger plot window
par(mfrow = c(1, 2))
plot(ndvi, col = rev(terrain.colors(10)), main = "Landsat-NDVI")
plot(knr, main = 'Unsupervised classification', col = mycolor, type = "classes")


#practical8################################################

#1
install.packages("geostan")
library(geostan)
library(ggplot2)
library(gridExtra)
data("georgia")
sp_diag(georgia$college, georgia, name = "College (%)")
W <- shape2mat(georgia, style = "W")
#2
moran_plot(georgia$college, W)
mc(georgia$college, W)

A <- shape2mat(georgia, "B")
#3
moran_plot(georgia$college, A)
x <- georgia$college
W <- shape2mat(georgia, "W")
mc(x, W)
gr(x, W)
mc(x, W) + gr(x, W)
W <- shape2mat(georgia, "W")
x <- log(georgia$income)
Ii <- lisa(x, W)
head(Ii)
Ci <- lg(x, W)
head(Ci)
#4
Ci_map <- ggplot(georgia) + geom_sf(aes(fill=Ci)) +
  scale_fill_gradient(high = "navy",low = "white") + theme_void()
Li_map <- ggplot(georgia) +
  geom_sf(aes(fill=Ii$Li)) + 
  scale_fill_gradient2(name = "Ii") +
  theme_void()
gridExtra::grid.arrange(Ci_map, Li_map, nrow = 1)
x <- log(georgia$income)
rho <- aple(x, W)
n <- nrow(georgia)
ess <- n_eff(rho = rho, n = n)
c(nominal_n = n, rho = rho, MC = mc(x, W), ESS = ess)
C <- shape2mat(georgia)

fit <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)),
                data = georgia,
                re = ~ GEOID,
                family = poisson(),
                C = C,
                refresh = 0 )
print(fit) 
#5
sp_diag(fit, georgia)


