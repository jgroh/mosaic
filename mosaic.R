library(RColorBrewer)
library(geosphere)
library(MASS)
library(ggmap)
library(ggplot2)
library(extrafont)

######## Read data files ########

herb <- read.csv("data/hrbm-pheno-mosaic.csv", stringsAsFactors = FALSE)
coll <- read.csv("data/ntrl-pops-pheno.csv", stringsAsFactors = FALSE)
porc <- read.csv("data/porcupine-vlys-2018-pheno.csv", stringsAsFactors = FALSE)
color <- read.csv("data/sep-color-mosaic.csv", stringsAsFactors = FALSE)

# columns with traits used for analysis
cols <- intersect(intersect(colnames(herb), colnames(coll)), colnames(porc))

# for indexing herbarium data by species
fo <- which(herb[, 1] == "formosa")
fl <- which(herb[, 1] == "flavescens")


######## Pool all data for PCA ########

pooled.dat <- rbind.data.frame( herb[, cols], coll[, cols], porc[, cols])

# leave out categorical trait for PCA, not indicative of size
PCA.pooled <- prcomp(pooled.dat[, -c(1, which(cols == "blade.cleft"))], scale = TRUE)

plot(PCA.pooled$x[,2] ~ PCA.pooled$x[,1], 
     col = as.numeric(as.factor(pooled.dat$species.label)), lwd = 3,
     xlab = "PC1 (43.5% of variance)", ylab = "PC2 (30.4 % of variance)")
PCA.pooled$rotation #PC2 has all negative loadings, captures size variation

legend('topright', legend = c("flavescens", "formosa", "hybrid"), 
       col = 1:3, pch = 1, bty = "n", pt.lwd = 3)

# calculate overall size for use later on in size-correction
overall.size <- mean(PCA.pooled$x[, 2])


######## Define training data set for fitting discriminant function ########


# calculate minimum distances to a record of the other species

min.fo2fl <- vector()
min.fl2fo <- vector()

lon <- herb$Geo_LongDegree
lat <- herb$Geo_LatDegree

fo <- which(herb$species.label== "formosa")
fl <- which(herb$species.label=="flavescens")

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
v <- min.dist

v[fl][ order( min.dist[fl])[1:51]] <- NA
v[fo][ order( min.dist[fo])[1:51]] <- NA

# next remove any that were a priori putative hybrids

v[which(herb[, "comment"] == "putative hybrid")] <- NA

length(which(is.na(v))) #additional 7 putative hybrids removed 

#finalize training and prediction sets
ind.trn <- which(!is.na(v))
trn.set <- herb[ind.trn, ] 


# PCA to examine species separation in training set of herbarium specimens

PCA.trn <- prcomp(trn.set[, c(12:15, 17:19)], scale. = TRUE)

col.trn <- vector()
col.trn[which(trn.set[, "species.label"] == "formosa")] <- "firebrick1"
col.trn[which(trn.set[,"species.label"] == "flavescens")] <- "yellow3"
plot(PCA.trn$x[,2] ~ PCA.trn$x[, 1], pch = 1, col = col.trn, lwd = 3, 
     xlab = "PCA1", ylab = "PCA2")

PCA.trn$rotation #PC2 again captures size variation


######## fit LDA model ########

LDA.trn <- lda(species.label ~ ., data = trn.set[, c(3,  12:19)])


######## Define Size-Correction Functions for Natural Populations and Herbarium Specimens ########

size.correct.pops <- function(x, scaling = "standard"){ 
  
  trait.data <- x[, c("corolla.width", "spur.length","blade.length", "blade.width", 
                      "sepal.length", "sepal.width", "anther.exsertion", "blade.cleft")]
  x$site <- as.factor(x$site)
  
  # if( !is.numeric(trait.data$blade.cleft)){ #reformat variable
  #   b <- toupper(trait.data$blade.cleft)
  #   for(i in 1:nrow(trait.data)){
  #     if( is.na(b[i])){
  #       b[i] <- NA
  #     } else if (!is.na(b[i]) & (b[i] == "Y" | b[i] == "1")){
  #       b[i] <- 1 
  #     } else if (!is.na(b[i]) & (b[i] == "N" | b[i] == "0")){
  #       b[i] <- 0
  #     }
  #   }
  #   trait.data$blade.cleft <- as.numeric(b)
  # }
  
  #get size scores for input data according to PC2 values
  size.scores <- predict(PCA.pooled, newdata = trait.data[, -which(colnames(trait.data) == "blade.cleft")])[, 2] 
  
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
  sigmoid <- function(x){1/(1+ exp(-x))} 
  
  if(scaling == "sigmoid"){
    lda.adjusted <- sigmoid(lda.adjusted)
  } else{
    lda.adjusted <- (lda.adjusted - min(lda.adjusted))/(max(lda.adjusted) - min(lda.adjusted))
  }
  
  return(lda.adjusted)
}

###### For Herbarium Specimens...

size.correct.herbarium <- function(x, prior = "equal", scaling = "standard"){
  
  trait.data <- x[, c("corolla.width", "spur.length","blade.length", "blade.width", 
                      "sepal.length", "sepal.width", "anther.exsertion", "blade.cleft")]

  x$species.label <- as.factor(x$species.label)
  

  #get size scores for input data according to PC2 values
  size.scores <- predict(PCA.pooled, newdata = trait.data[, -which(colnames(trait.data) == "blade.cleft")])[, 2] 
  
  #predict lda scores based on LDA of all specimens or of training set
 
  #set priors
   priors <- list()
   
   if(prior == "unequal"){
  
      for(i in 1:nrow(x)){
        
        if( x[i , "species.label"] == "formosa"){
          
          priors[[i]] <- c(0.75, 0.25)
          
        } else if( x[i, "species.label"] == "flavescens"){
          
          priors[[i]] <- c(0.25, 0.75)
        } else{ priors[[i]] <- c(0.5, 0.5)}
      }
     
   } else{
     for(i in 1:nrow(x)){
       priors[[i]] <- c(.5,.5)
     }
   }
  
  lda.raw.scores <- vector()
  for (i in 1:nrow(x)){
    lda.raw.scores[i] <- predict(LDA.trn, newdata = trait.data[i,], prior = priors[[i]])$x
  }
  
  size.models <- list(length = length(levels(x$species.label)))
  
  for(i in 1:length(levels(x$species.label))){
    spc.ind <- which(x$species.label == levels(x$species.label)[i]) #index of all individuals in that pop
    
    spc.set <- cbind.data.frame(lda = lda.raw.scores[spc.ind], size = size.scores[spc.ind])
    
    spc.lm <- lm(lda ~ size, data = spc.set) #regression of predicted lda scores against size 
    size.models[[i]] <- spc.lm
  }
  
  lda.adjusted <- vector()
  for(i in 1:nrow(x)){ # loop projects individuals to a common size, then displaces lda score from the lsmean by the residual from population regressionline at its original size
    
    spc <- x$species.label[i]
    spc.num <- which(levels(x$species.label) == spc)
    
    pred.dat <- cbind.data.frame(size = size.scores[i])
    # if( x[i, "species"] == "flavescens" ){
    
    lda.predicted <- predict( size.models[[spc.num]], newdata = pred.dat)
    residual <- lda.raw.scores[i] - lda.predicted
    lda.adjusted[i] <- predict( size.models[[spc.num]], newdata = data.frame(size = overall.size)) + residual
    
  }
  
  if(scaling == "sigmoid"){
    sigmoid <- function(x){1/(1+ exp(-x))}
    lda.adjusted <- sigmoid(lda.adjusted)
    
  } else{
    lda.adjusted <- (lda.adjusted - min(lda.adjusted))/(max(lda.adjusted) - min(lda.adjusted))
  }
  return(list(lda.raw.scores, lda.adjusted))
}


######### FIGURE 1: Hybrid Index Heat Map ########

lon <- herb[,"Geo_LongDegree"]
lat <- herb[,"Geo_LatDegree"]

# Prepare data for ggmap 

plot.data <- herb[, c("species.label", "Geo_LongDegree", "Geo_LatDegree")]
herb.scores <- size.correct.herbarium(herb, prior = "equal", scaling = "sigmoid")[[2]]
plot.data$fit <- herb.scores

# combine collection data sets to add population means to heat map plot

cols2 <- intersect(colnames(coll), colnames(porc))
all.coll.dat <- rbind(coll[, cols2], porc[ grep("porcupine", porc$site), cols2])
all.coll.dat[ grep("porcupine", all.coll.dat$site), ]$site <- rep("Porcupine Ridge",  length(grep("porcupine", all.coll.dat$site)))
#all.coll.dat[ grep("adj", all.coll.dat$site), ]$site <- rep("Adjacent Valley",  length(grep("adj", all.coll.dat$site)))

pop.scores <- size.correct.pops(all.coll.dat, scaling = "sigmoid")

# quick look at hybrid index values by population
stripchart(pop.scores ~ all.coll.dat$site, vertical = TRUE, method = "jitter", las = 2)

# calculate popualtion means
pop.means <- tapply(pop.scores, all.coll.dat$site, mean)
pop.means

# write corresponding coordinates based on order of population means 
pop.lons <- c(-117.546, -120, -120.918, -120.3998, -119.666, -121.828245, -125.546)
pop.lats <- c(51.266, 51.9, 40.05, 47.2924, 49.111, 51.112431, 50.215)
pop.sp <- c("flavescens", "formosa", "formosa", "flavescens", "flavescens", "hybrid", "formosa")

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
       shape = "Verbatim species \ndetermination") 

print(heat.map)


# # heat map of sepal color
# 
# herb2 <- read.csv("data/herbarium-phenotypes-mosaic.csv", stringsAsFactors = FALSE)
# herb2$sepal.color.0.4[ which(herb2$sepal.color.0.4 == 3)] <- 4
# 
# # add points colored according to fitted values
# heat.map <- ggmap(map) + geom_point(data = herb2,
#                                     aes(x = Geo_LongDegree, 
#                                         y = Geo_LatDegree, shape = species.label, 
#                                         fill = sepal.color.0.4), cex = 2.5, stroke = .3) + 
#   scale_shape_manual(values = c(22, 21, 23)) + 
#   scale_fill_gradient(low = "#F5FF0A", high = "#F50000", limits = c(0,4), breaks = c(0, 1, 2, 3, 4)) +
#   labs(title = "", x = "Longitude", 
#        y = "Latitude", fill = "Hybrid index", 
#        shape = "Verbatim species \ndetermination") 
# 
# print(heat.map)



######### FIGURE 2: Geospatial Morphology Diagram ########

herb.scores.st <- size.correct.herbarium(herb, prior = "equal", scaling = "standard")[[2]]
log.min.dist <- log(min.dist)
colvec <- vector()
colvec[fo] <- "#F50000"
colvec[fl] <- "#F5FF0A"

#set up plot
dev.off()
par(mar = c(5,5,4,2))
plot(herb.scores.st ~ log.min.dist, type = "n",
     xlab = "Log(min km to record of other species)", col = colvec,
     ylab = "Scaled discriminant score", yaxt = "n", cex.lab = 1.2, cex.axis = 1.2, 
     xlim = c(0, 8))
axis(2, las = 2, cex.axis = 1.2)

#make models
fl.dat <- cbind.data.frame(morph = herb.scores.st[fl], log.min.dist = log.min.dist[fl])
flz <- lm(morph ~ log.min.dist, dat = fl.dat)
fo.dat <- cbind.data.frame(morph = herb.scores.st[fo], log.min.dist = log.min.dist[fo])
foz <- lm(morph ~ log.min.dist, dat = fo.dat)

# add confidence bands
xpfl2 <- seq(min(fl.dat$log.min.dist), max(fl.dat$log.min.dist), length.out=100)
prfl <- predict(flz, newdata = data.frame(log.min.dist = xpfl2),  interval = 'confidence')
polygon(c(rev(xpfl2), xpfl2), c(rev(prfl[ ,3]), prfl[ ,2]), col = rgb(190/255,190/255,190/255, .6), border = NA)

xpfo2 <- seq(min(fo.dat$log.min.dist), max(fo.dat$log.min.dist), length.out=100)
prfo <- predict(foz, newdata = data.frame(log.min.dist = xpfo2),  interval = 'confidence')
polygon(c(rev(xpfo2), xpfo2), c(rev(prfo[ ,3]), prfo[ ,2]), col = rgb(190/255,190/255,190/255, .6), border = NA)

#add regression lines
xpfl <- range(fl.dat$log.min.dist)
ypfl <- predict(flz, newdata = data.frame(log.min.dist = xpfl))
lines(x = xpfl, y = ypfl)

xpfo <- range(fo.dat$log.min.dist)
ypfo <- predict(foz, newdata = data.frame(log.min.dist = xpfo))
lines(x = xpfo, y = ypfo)

# add points
points(herb.scores.st ~ log.min.dist, pch = 21, bg = colvec)


######## FIGURE 3A: Population Hybrid Index Stripchart #######

all.coll.dat2 <- rbind(coll[, cols2], porc[, cols2])

# index sites
por <- grep("porcupine", all.coll.dat2$site)
adj <- grep("adj", all.coll.dat2$site)
all.coll.dat2[por ,]$site <- rep("Porcupine Ridge",  length(por))
all.coll.dat2[adj, ]$site <- rep("Adjacent valley", length(adj))
levels(as.factor(all.coll.dat2$site))

pop.scores <- size.correct.pops(all.coll.dat2, scaling = "sigmoid")

all.coll.dat2$site <- factor(all.coll.dat2$site, 
                       levels = levels(as.factor(all.coll.dat2$site) )[order( tapply( 
                         pop.scores, all.coll.dat2$site, mean) ) ] )

colvec8 <- brewer.pal(8, "Dark2")
colvec8.1 <- colvec8[c(6, 2,3,4,5, 8, 7, 1)]

par(mar = par()$mar + c(0,1,0,0), bty = "o")
stripchart(pop.scores ~ all.coll.dat2$site, vertical = TRUE, method = "jitter", 
           col = colvec8.1, pch = 1, lwd = 2, axes = "FALSE", 
           cex.lab = 1.2, ylab = "")
mtext(2, text = "Hybrid index", line = 3.5, cex = 1.2)
axis(2, lwd = 0, lwd.ticks = 1, las = 2, cex.axis = 1.2)
legend("bottomright", levels(new.dat$site), col = colvec8.1, pch = 1, pt.lwd = 2, bty = "n")
box()


######## FIGURE 3B: Population Sepal Color Stripchart ########


color$log.rg <- with(color, log(red.mean / green.mean))

color$site <- factor(color$site, 
                    levels = levels(as.factor(color$site) )[order( tapply( 
                      color$log.rg, color$site, mean) ) ] )
colvec8.2 <- colvec8.1[c(1,2,3,4,5, 8, 6)]
stripchart(log.rg ~ site, data = color, vertical = TRUE, pch = 1, las = 2, 
           method = "jitter" , ylab = "Sepal red vs. green reflectance (log R/G)", 
           col = colvec8.2, xaxt = "n", yaxt = "n", lwd = 2, cex.axis = 1.2, 
           cex.lab = 1.2, ylim = c(-.05,.8))
axis(2, las = 2, cex.axis = 1.2)

# legend('topleft', legend = levels(as.factor(color$site)), 
#        col = colvec8.2, pch = 1, pt.lwd = 2, bty = "n")


######## Fig 4: A. flavescens var. miniana Type Specimen Analysis ########

mini <- read.csv("data/miniana.csv", stringsAsFactors = FALSE)

type.scores <- size.correct.herbarium(mini, prior = "equal", scaling = "sigmoid")[[2]]

dev.off()
par(mar = par()$mar + c(0,1,0,0), bty = "o")
plot.new()
plot.window(xlim = c(.15,.85), ylim =c(0,1))

points(y = type.scores[which(mini$type == "paratype")], 
       x = jitter(rep(.25, 4), 3), pch = 1, lwd = 2)

points(y = type.scores[which(mini$type == "type")], 
       x = .5, pch = 1, lwd = 2)

points(y = type.scores[which(mini$type == "isotype")], 
       x = jitter(rep(.75, 8), 2), pch = c(2,2, rep(1, 6)), lwd = 2)

axis(2, lwd = 0, lwd.ticks = 1, las = 2, cex.axis = 1.2)
axis(1, lwd = 0, lwd.ticks = 1, at = c(.25, .5, .75), labels = c("paratypes", "holotype", "isotypes"))
box()
mtext(2, text = "Hybrid index", line = 3.5, cex = 1.2)


#### END