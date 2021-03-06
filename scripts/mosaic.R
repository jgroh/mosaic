library(geosphere)
library(MASS)
library(ggmap)
library(ggplot2)
library(visreg)

# ----- Read data files -----

herb <- read.csv("data/hrbm-pheno-mosaic.csv", stringsAsFactors = FALSE)
coll <- read.csv("data/ntrl-pops-pheno.csv", stringsAsFactors = FALSE)
porc <- read.csv("data/porcupine-vlys-2018-pheno.csv", stringsAsFactors = FALSE)
mini <- read.csv("data/miniana.csv", stringsAsFactors = FALSE)
color <- read.csv("data/sep-color-mosaic.csv", stringsAsFactors = FALSE)

# columns with traits used for analysis
cols <- c("species.label", "anther.exsertion", 
          "corolla.width", "spur.length", "blade.length", 
          "blade.width", "sepal.length", "sepal.width")
 
# objects for indexing herbarium data by species
fo <- which(herb[, "species.label"] == "formosa")
fl <- which(herb[, "species.label"] == "flavescens")


# ----- PCA on all data-----

pooled.dat <- rbind.data.frame( herb[, cols], coll[, cols], porc[, cols])
pca.pooled.input <- pooled.dat[,-1]


# run PCA
PCA.pooled <- prcomp(pca.pooled.input, scale. = T)
summary(PCA.pooled)

# plot
plot(PCA.pooled$x[,2] ~ PCA.pooled$x[,1], 
     col = as.numeric(as.factor(pooled.dat$species.label)), lwd = 3,
     xlab = "PC1 (43.5% of variance)", ylab = "PC2 (30.4 % of variance)")
legend('topleft', legend = c("flavescens", "formosa", "hybrid"), 
       col = 1:3, pch = 1, bty = "n", pt.lwd = 3)

# inspect trait loadings
PCA.pooled$rotation #PC2 has all negative loadings, captures size variation

# calculate overall size for use later on in size-correction
overall.size <- mean(PCA.pooled$x[, 2])


# ----- Define training data set for fitting discriminant function -----

# calculate minimum distances to a record of the other species
min.fo2fl <- vector()
min.fl2fo <- vector()

lon <- herb$Geo_LongDegree
lat <- herb$Geo_LatDegree

for(i in 1:length(fo)){ #calculates min distance to record of flavescens for each formosa
  fo2fl <- vector()
  for(j in 1:length(fl)) {
    fo2fl[j] <- distHaversine(c(lon[fo][i], lat[fo][i]), 
                              c(lon[fl][j], lat[fl][j]) )
  }
  min.fo2fl[i]  <- min(fo2fl)/1000
}

for(i in 1:length(fl)){ #calculates min distance to record of formosa for each flavescens
  fl2fo <- vector()
  for(j in 1:length(fo)) {
    fl2fo[j] <- distHaversine(c(lon[fl][i], lat[fl][i]), 
                              c(lon[fo][j], lat[fo][j]) )
  }
  min.fl2fo[i]  <- min(fl2fo)/1000
}

min.dist <- vector()
min.dist[fo] <- min.fo2fl
min.dist[fl] <- min.fl2fo


# omit 40% from each species data set with closest location to a record of the other species
# This value was chosen by visualizing the training set in the PCA (see below) and choosing
# a value which results in nearly non-overlapping species phenotype distributions
v <- min.dist

v[fl][ order( min.dist[fl])[1:51]] <- NA
v[fo][ order( min.dist[fo])[1:51]] <- NA

# next remove any that were a priori putative hybrids
v[which(herb[, "comment"] == "putative hybrid")] <- NA
length(which(is.na(v))) #additional 7 putative hybrids removed 

#finalize training and prediction sets
ind.trn <- which(!is.na(v))

lda.input <- herb[ind.trn, cols] 

dim(lda.input)
table(lda.input$species.label)


# Examine separation training data set
g <- predict(PCA.pooled, newdata = pooled.dat[ind.trn,])
plot(g[,2] ~ g[,1], col = as.factor(pooled.dat[ind.trn,]$species.label), lwd=2)

#check whether size scores are significantly different between species in the training set
gfo <- g[,2][pooled.dat[ind.trn,]$species.label == "formosa"]
gfl <- g[,2][pooled.dat[ind.trn,]$species.label == "flavescens"]
t.test(gfo, gfl)
# ----- fit LDA model -----
LDA.trn <- lda(species.label ~ ., data = lda.input)
plot(LDA.trn)

# trait loadings and group means
LDA.trn

# compute standard errors
se.fo <- apply(lda.input[lda.input$species.label == "formosa", -1], 
            MARGIN = 2, FUN = function(x) {sd(x)/sqrt(length(x))} )
se.fo
se.fl <- apply(lda.input[lda.input$species.label == "flavescens", -1], 
            MARGIN = 2, FUN = function(x) {sd(x)/sqrt(length(x))} )
se.fl

# for calculating classification accuracy later
lda.input$species.label <- as.factor(lda.input$species.label)



# ----- Define size-correction function for natural populations  -----
sigmoid <- function(x){1/(1+ exp(-x))} 

size.correct.pops <- function(x, scaling = "standard"){ 
  
  trait.data <- x[, cols[-1]]
  
  x$site <- as.factor(x$site)
  
  #get size scores for input data according to PC2 values
  size.scores <- predict(PCA.pooled, newdata = trait.data)[, 2] 
  
  #predict lda scores based on LDA of training set
  lda.raw.scores <- predict(LDA.trn, newdata = trait.data, prior = c(.5,.5))$x
  
  #adjust scores for size using linear regression
  lda.adjusted <- vector()
  
  size.models <- list(length = length(levels(x$site))) #regression of lda ~ size for each population
  
  for(i in 1:length(levels(x$site))){
    pop.ind <- which(x$site == levels(x$site)[i]) #index of all individuals in that pop
    
    pop.set <- cbind.data.frame(lda = lda.raw.scores[pop.ind], size = size.scores[pop.ind])
    
    pop.lm <- lm(lda ~ size, data = pop.set) #regression of predicted lda scores against size 
    size.models[[i]] <- pop.lm
  }
  
  for(i in 1:nrow(x)){ # loop projects individuals to a common size, then displaces lda score from the lsmean by the residual from population regression line at its original size
    
    pop <- x$site[i]
    pop.num <- which(levels(x$site) == pop)
    
    pred.dat <- cbind.data.frame(size = size.scores[i])
    
    lda.predicted <- predict( size.models[[pop.num]], newdata = pred.dat) #predicted lda score at given size according to regression line
    residual <- lda.raw.scores[i] - lda.predicted
    lda.adjusted[i] <- predict( size.models[[pop.num]], newdata = data.frame(size = overall.size)) + residual
    
  }
  
  #re-scale values to range between 0 and 1
  
  if(scaling == "sigmoid"){
    lda.adjusted <- sigmoid(lda.adjusted)
  } else{
    lda.adjusted <- (lda.adjusted - min(lda.adjusted))/(max(lda.adjusted) - min(lda.adjusted))
  }
  return(lda.adjusted)
}


# ----- Define Size-correction for herbarium specimens -----

size.correct.herbarium <- function(x, scaling = "standard"){  
  
  trait.data <- x[, cols[-1]]
  

  x$species.label <- as.factor(x$species.label)
  
  #get size scores for input data according to PC2 values
  size.scores <- predict(PCA.pooled, newdata = trait.data)[, 2] 
  
  #predict lda scores based on LDA of all specimens or of training set
  predictions <-  predict(LDA.trn, newdata = trait.data, prior = c(.5,.5))
  lda.raw.scores <- predictions$x
  
  size.models <- list(length = length(levels(x$species.label)))
  
  for(i in 1:length(levels(x$species.label))){
    spc.ind <- which(x$species.label == levels(x$species.label)[i]) #index of all individuals in that pop
    
    spc.set <- cbind.data.frame(lda = lda.raw.scores[spc.ind], size = size.scores[spc.ind])
    
    spc.lm <- lm(lda ~ size, data = spc.set) #regression of predicted lda scores against size 
    size.models[[i]] <- spc.lm
  }
  
  lda.adjusted.herb <- vector()
  for(i in 1:nrow(x)){ # loop projects individuals to a common size, then displaces lda score from the lsmean by the residual from population regression line at its original size
    
    spc <- x$species.label[i]
    spc.num <- which(levels(x$species.label) == spc)
    
    pred.dat <- cbind.data.frame(size = size.scores[i])
    # if( x[i, "species"] == "flavescens" ){
    
    lda.predicted <- predict( size.models[[spc.num]], newdata = pred.dat)
    residual <- lda.raw.scores[i] - lda.predicted
    lda.adjusted.herb[i] <- predict( size.models[[spc.num]], newdata = data.frame(size = overall.size)) + residual
    
  }
  
  if(scaling == "sigmoid"){
    sigmoid <- function(x){1/(1+ exp(-x))}
    hybrid.index.herb <- sigmoid(lda.adjusted.herb)
    
  } else{
    hybrid.index.herb <- (lda.adjusted.herb - min(lda.adjusted.herb))/(max(lda.adjusted.herb) - min(lda.adjusted.herb))
  }
  return(list(lda.raw.scores = lda.raw.scores, 
              predicted.group = predictions$class, 
              lda.adjusted.herb = lda.adjusted.herb, 
              hybrid.index.herb = hybrid.index.herb))
}


# ----- FIGURE 1: Hybrid Index Heat Map -----

# Prepare data for ggmap 

plot.data <- herb[, c("species.label", "Geo_LongDegree", "Geo_LatDegree")]
herb.sc <- size.correct.herbarium(herb, scaling = "sigmoid")
herb.scores <- herb.sc$hybrid.index.herb


# check predictive accuracy
mean(as.character(herb.sc$predicted.group) == herb$species.label)

plot.data$fit <- herb.scores

# combine collection data sets to add population means to heat map plot

cols3 <- intersect(colnames(coll), colnames(porc))
all.coll.dat <- rbind(coll[, cols3], porc[ grep("porcupine", porc$site), cols3])
all.coll.dat[ grep("porcupine", all.coll.dat$site), ]$site <- rep("Porcupine Ridge",  length(grep("porcupine", all.coll.dat$site)))
pop.scores <- size.correct.pops(all.coll.dat, scaling = "sigmoid")

# quick look at hybrid index values by population
stripchart(pop.scores ~ all.coll.dat$site, vertical = TRUE, method = "jitter", las = 2)


# write corresponding coordinates based on order of population means 
pop.lons <- c(-117.546, -120, -120.918, -120.3998, -119.666, -121.828245, -125.546)
pop.lats <- c(51.266, 51.9, 40.05, 47.2924, 49.111, 51.112431, 50.215)
pop.sp <- c("flavescens", "formosa", "formosa", "flavescens", "flavescens", "hybrid", "formosa")

pop.means <- tapply(pop.scores, all.coll.dat$site, mean)

pop.df <- cbind.data.frame(species.label = pop.sp, Geo_LongDegree = pop.lons, 
                           Geo_LatDegree = pop.lats, fit = pop.means)
all.plot.data <- rbind.data.frame(plot.data, pop.df)


# reorder data so that intermediate colored points are not covered
hy.ind <- all.plot.data$fit

i0 <- which(hy.ind > 0.9)
i1 <- which(hy.ind < 0.1)
i2 <- which(( hy.ind < 0.9 & hy.ind > 0.8) | (hy.ind > 0.1 & hy.ind < 0.2))
i3 <- which(( hy.ind < 0.8 & hy.ind > 0.7) | (hy.ind > 0.2 & hy.ind < 0.3))
i4 <- which( (hy.ind < 0.7 & hy.ind > 0.6) | ( hy.ind > 0.3 & hy.ind < 0.4))
i5 <- which( hy.ind < 0.6 & hy.ind > 0.4)
 
all.plot.data <- all.plot.data[c(i0, i1, i2, i3, i4, i5), ]

# plot background terrain map
map <- get_map(location = c(-154, 30, -108, 63), source = "stamen", maptype = "terrain-background")


# add points colored according to fitted values
heat.map <- ggmap(map) + geom_point(data = all.plot.data,
                                    aes(x = Geo_LongDegree, 
                                        y = Geo_LatDegree, shape = species.label, 
                                        fill = fit), cex = 2.5, stroke = .3) + 
  scale_shape_manual(values = c(22,21, 23)) + 
  scale_fill_gradient(low = "#F5FF0A", high = "#F50000", limits = c(0,1), breaks = c(0, .25, 0.5, 0.75, 1)) +
  labs(title = "", x = "Longitude", 
       y = "Latitude", fill = "Hybrid index", 
       shape = "Verbatim species \ndetermination") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15), 
        legend.text=element_text(size=12), legend.title=element_text(size=15))

print(heat.map)
dev.off()


# ----- FIGURE 2: Geospatial Morphology Diagram -----

# rescale discriminant scores between 0 and 1
herb.scr <- (herb.sc$lda.raw.scores - min(herb.sc$lda.raw.scores))/(max(herb.sc$lda.raw.scores - min(herb.sc$lda.raw.scores)))

# get size sores
size.scores <- predict(PCA.pooled, newdata = herb)[, 2] 

# distance scores
log.min.dist <- log(min.dist)

# data frame
d <- cbind.data.frame(LD1 = herb.scr, size= size.scores, dist = log.min.dist, species = herb$species.label)

# model
z <- lm(LD1 ~ species*dist + size*species, data = d)

# plot
par(oma = c(0,1,0,0))
visreg(z, xvar = "dist", by = "species", overlay = T,
       line = list(col=c("black", "black"), lwd = 3),
       fill = list(col=c(grey(.5, alpha = .4))),
       points = list(bg=c("darkgray", "black"), cex = 1, pch=c(22,21),col=c("black","black"), lwd = .5),
       xlab = "Log(min km to record of other species)", 
       ylab = "Hybrid index",
       legend = F,
       cex.lab = 1.2,
       cex.axis = 1.2,
       ylim = c(0,1)
       )



# ----- FIGURE 3A: Population Hybrid Index Stripchart -----

par(mfrow=c(1,2), mar = c(6,4,6,1))

# combine data sets
cols4 <- cols3
cols4[1] <- colnames(mini)[1]
cols4[3] <- colnames(mini)[5]
mini4 <- mini[, cols4]
colnames(mini4) <- cols3

all.coll.dat2 <- rbind(coll[, cols3], porc[, cols3], mini4)
all.coll.dat2 <- all.coll.dat2[all.coll.dat2$site != "Manning Park" & all.coll.dat2$site != "Clearwater",]


# index sites
por <- grep("porcupine", all.coll.dat2$site)
adj <- grep("adj", all.coll.dat2$site)
all.coll.dat2[por ,]$site <- rep("Porcupine Ridge",  length(por))
all.coll.dat2[adj, ]$site <- rep("Adjacent valley", length(adj))
levels(as.factor(all.coll.dat2$site))
all.coll.dat2 <- all.coll.dat2[-which(all.coll.dat2$site == "Adjacent valley"),]
#all.coll.dat2$site

# get size-corrected hybrid index values for natural populations
pop.scores <- size.correct.pops(all.coll.dat2, scaling = "sigmoid")
all.coll.dat2$site <- factor(all.coll.dat2$site, 
                       levels = levels(as.factor(all.coll.dat2$site) )[order( tapply( 
                         pop.scores, all.coll.dat2$site, mean) ) ] )

# make plot
par(mar = par()$mar + c(0,1,0,0), bty = "o")
stripchart(pop.scores ~ all.coll.dat2$site, vertical = TRUE, method = "jitter", 
          pch = 1, lwd = 1.5, axes = "FALSE", 
           cex.lab = 1.2, ylab = "")

# population means
y  <-  tapply(pop.scores, all.coll.dat2$site, mean)

for(i in 1:6){
  lines(x = c(i-.35, i+.35), y = c(y[i], y[i]),  lwd = 2)
}

mtext(2, text = "Hybrid index", line = 3.5, cex = 1.2)
axis(2, lwd = 0, lwd.ticks = 1, las = 2, cex.axis = 1.2)
axis(1, at = 1:6, labels = c(1:5,7))
box()

#Which specimens annotated by Whittemore?
cbind(pop.scores[all.coll.dat2$site == "var. miniana"], mini$annotated)


# ----- FIGURE 3B: Population Sepal Color Stripchart ----- 
color <- color[color$site != "Clearwater",]

color$log.rg <- with(color, log(red.mean / green.mean))

color$site <- factor(color$site, 
                    levels = levels(as.factor(color$site) )[order( tapply( 
                      color$log.rg, color$site, mean) ) ] )

stripchart(log.rg ~ site, data = color, vertical = TRUE, pch = 1, las = 2, 
           method = "jitter" , ylab = "Sepal red vs. green reflectance (log R/G)", 
          xaxt = "n", yaxt = "n", lwd = 1.5, cex.axis = 1.2, 
           cex.lab = 1.2, ylim = c(-.05,.8))
y2 <- tapply(color$log.rg, color$site, mean)
for(i in 1:6){
  lines(x = c(i-.35, i+.35), y = c(y2[i], y2[i]),  lwd = 2)
}

axis(2, las = 2, cex.axis = 1.2)
axis(1, at = 1:6, labels = c(1,2,4,5,6,7))

# confirm order of populations on x axis
levels(color$site)

dev.off()


# ----- FIGURE 4: Elevational Morphology Cline ----- 

s <- predict(PCA.pooled, newdata = porc)[,2]
l <- predict(LDA.trn, newdata = porc)$x
hybrid <- (l - min(herb.sc$lda.adjusted.herb))/(max(herb.sc$lda.adjusted.herb) - min(herb.sc$lda.adjusted.herb))
w <- cbind.data.frame(size = s, l = l, hybrid = hybrid, elev = porc$elevation)

# fit model with size as a covariate
z <- lm(hybrid ~ elev + size, data = w)
summary(z)
par(oma = c(0,1,0,0))
visreg(z, xvar = "elev", points=list(pch = 1, lwd = 2, cex = 1, col = "black", cex = 1),
       fill.par=list(col = "lightgray"), line.par = list(col= "black"),cex.lab= 1.2,
       xlab = "Elevation (m)", ylab = "Hybrid index", cex.axis = 1.2)
anova(z)



