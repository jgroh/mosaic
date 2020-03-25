library(geosphere)
library(MASS)

# ===== Read Data Files =====
coll.dat <- read.csv("data/transect-pheno.csv", stringsAsFactors = FALSE)
coll.dat <- subset(coll.dat, comments != "developmental mutant")
hrbm.dat <- read.csv("data/hrbm-pheno.csv", stringsAsFactors = FALSE)
prev.dat <- read.csv("data/prev-coll-pheno.csv", stringsAsFactors = FALSE)

# ===== PCA on Pooled Data Set =====

# columns with traits used for analysis
cols <- intersect(intersect(colnames(hrbm.dat), colnames(coll.dat)), colnames(prev.dat))

# pool data based on shared columns
pooled.dat <- rbind.data.frame(hrbm.dat[, cols], coll.dat[, cols], prev.dat[, cols])

# leave out categorical trait for PCA, not indicative of size
PCA.pooled <- prcomp(pooled.dat[, -c(1, which(cols == "blade.cleft"))], scale = TRUE)
summary(PCA.pooled)

plot(PCA.pooled$x[,2] ~ PCA.pooled$x[,1], 
     col = as.numeric(as.factor(pooled.dat$species.label)), lwd = 3,
     xlab = "PC1 (41.9% of variance)", ylab = "PC2 (31.4 % of variance)")
PCA.pooled$rotation #PC2 has all negative loadings, captures size variation

legend('topright', legend = c("flavescens", "formosa", "hybrid"), 
       col = 1:3, pch = 1, bty = "n", pt.lwd = 3)

# calculate overall size for use later on in size-correction
overall.size <- mean(PCA.pooled$x[, 2])


# ===== Define Training Set for Fitting Discriminant Function =====

# 1. Filter by distance
# calculate minimum distances to a record of the other species

min.fo2fl <- vector()
min.fl2fo <- vector()

lon <- hrbm.dat$Geo_LongDegree
lat <- hrbm.dat$Geo_LatDegree

fo <- which(hrbm.dat$species.label== "formosa")
fl <- which(hrbm.dat$species.label=="flavescens")

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

v[which(hrbm.dat[, "comment"] == "putative hybrid")] <- NA

length(which(is.na(v))) #additional 7 putative hybrids removed 

#finalize training and prediction sets
ind.trn <- which(!is.na(v))
trn.set <- hrbm.dat[ind.trn, ] 

#trn.set <- hrbm.dat
# PCA to examine species separation in training set of herbarium specimens

PCA.trn <- prcomp(trn.set[, c(12:15, 17:19)], scale. = TRUE)

col.trn <- vector()
col.trn[which(trn.set[, "species.label"] == "formosa")] <- "firebrick1"
col.trn[which(trn.set[,"species.label"] == "flavescens")] <- "yellow3"
plot(PCA.trn$x[,2] ~ PCA.trn$x[, 1], pch = 1, col = col.trn, lwd = 3, 
     xlab = "PCA1", ylab = "PCA2")

PCA.trn$rotation #PC2 again captures size variation


# ===== Fit LDA model =====

colnames(trn.set)

LDA.trn <- lda(species.label ~ ., data = trn.set[, c(3,  12:15,19)])


# ===== Define Function to Calculate Size-Adjusted Hybrid Index =====

HybridIndex <- function(x){ 
  # This function will calculate predicted 'hybrid index' 
  # (LDA score from model defined using herbarium specimens)
  # for transect phenotypes and use linear regression to determine scaling relationship with size 
  # (value along PC2 from pooled data set).
  # Each individual will be displaced by its residual from the population regression line 
  # from the predicted discriminant score (using the population regression line) 
  # at the mean size from pooled data.
  
  # Input a data frame of trait data. 
  # Scaling can be either 'standard' or 'sigmoid.'
  
  trait.data <- x[, c("corolla.width", "spur.length","blade.length", "blade.width", 
                      "sepal.length", "sepal.width", "anther.exsertion", "blade.cleft")]
  
  # get size scores for input data according to PC2 values
  size.scores <- predict(PCA.pooled, newdata = trait.data[, -which(colnames(trait.data) == "blade.cleft")])[, 2] 
  
  # predict lda scores based on LDA of training set
  lda.raw.scores <- as.vector(predict(LDA.trn, newdata = trait.data, prior = c(.5,.5))$x)
  # *note as.vector was used to prevent error in lm 
  
  # empty vector for adjusted hybrid index scores 
  lda.adjusted <- vector()
  
  # calculate adjusted scores
  pop.set <- cbind.data.frame(lda = lda.raw.scores, size = size.scores)
  
  pop.lm <- lm(lda ~ size, data = pop.set) #regression of predicted lda scores against size 
  
  for(i in 1:nrow(x)){ # loop effectively projects individuals to a common size
    
    pred.dat <- cbind.data.frame(size = size.scores[i])
    
    lda.predicted <- predict(pop.lm, newdata = pred.dat) #predicted lda score at given size according to regression line
    residual <- lda.raw.scores[i] - lda.predicted
    lda.adjusted[i] <- predict( pop.lm, newdata = data.frame(size = overall.size)) + residual
    
  }
  
  #re-scale values to range between 0 and 1
  lda.adjusted <- (lda.adjusted - min(lda.adjusted))/(max(lda.adjusted) - min(lda.adjusted))
  #sigmoid <- function(x){1/(1+ exp(-x))} 
  #lda.adjusted <- sigmoid(lda.adjusted)
  
  return(lda.adjusted)
}


# calculate size-adjusted hybrid index scores for transect phenotypes
hyb <- HybridIndex(coll.dat)
hyb2 <- HybridIndex(coll.dat[which(coll.dat$site == "porcupine valley transect"),])

# regression of hybrid index against elevation
zdat <- cbind.data.frame(hyb, elev = coll.dat$elevation)
z <- lm(hyb ~ elev, data = zdat)

# make plot
plot(hyb ~ coll.dat$elevation, lwd = 2,
     ylab = "Hybrid index", xlab = "Elevation (m)", 
     cex.axis = 1.2, cex.lab = 1.2)

# add confidence bands
xpt <- seq(min(zdat$elev), max(zdat$elev), length.out=100)
prd <- predict(z, newdata = data.frame(elev = xpt),  interval = 'confidence')
polygon(c(rev(xpt), xpt), c(rev(prd[ ,3]), prd[ ,2]), col = rgb(190/255,190/255,190/255, .6), border = NA)

#add regression lines
xpt <- range(zdat$elev)
ypt <- predict(z, newdata = data.frame(elev = xpt))
lines(x = xpt, y = ypt, lwd = 3, col = "darkgreen")

# overlay points again 
points(zdat$hyb ~ zdat$elev, lwd = 2)




zdat2 <- cbind.data.frame(hyb2, elev = coll.dat$elevation[which(coll.dat$site == "porcupine valley transect")])
z2 <- lm(hyb2 ~ elev, data = zdat2)

plot(hyb2 ~ coll.dat$elevation[which(coll.dat$site == "porcupine valley transect")], lwd = 2,
     ylab = "Size-adjusted hybrid index", xlab = "Elevation (m)", 
     cex.axis = 1.2, cex.lab = 1.2)

xpt <- seq(min(zdat2$elev), max(zdat2$elev), length.out=100)
prd <- predict(z2, newdata = data.frame(elev = xpt),  interval = 'confidence')
polygon(c(rev(xpt), xpt), c(rev(prd[ ,3]), prd[ ,2]), col = rgb(190/255,190/255,190/255, .6), border = NA)

xpt <- range(zdat2$elev)
ypt <- predict(z2, newdata = data.frame(elev = xpt))
lines(x = xpt, y = ypt, lwd = 3, col = "darkgreen")

points(zdat2$hyb ~ zdat2$elev, lwd = 2)
anova(z2)













#predict species probabilities according to herbarium specimen logistic regression model
#coll.dat$pred.pops <- predict(z.log.reg, newdata = coll.dat[,c(5,6,7,9,12)], type = "response")

# PCA of collected specimen phenotypes
pca <- prcomp(coll.dat[ , c(5:12)], scale. = TRUE)

pca1 <- pca$x[,1]
pca2 <- pca$x[,2]

#inspect trait loadings
s <- summary(pca)
s
s$rotation 
#PC1 separates species, all traits load neg onto PC2, so this axis may represent "size"

#make biplot
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]

library(RColorBrewer)

dev.off()
cpal <- brewer.pal(7,"Accent")
cpal <- cpal[c(4, 5, 3, 6, 2, 1, 7)]
palette(cpal)

plot(pca2 ~ pca1, pch = 21, bg = as.factor(coll.dat$site),
     xlab = "PC1 (42.8% variance explained)", 
     ylab = "PC2 (26.2% variance explained")


# compare density of logistic regression response values from herb specimens to 
# field-collected specimens

#dev.off()
#par(mar = c(5,6,4,2))
#plot.new()
#plot(density(fit), ylim = c(0,1.5), col = "dodgerblue", lwd = 4, 
#     xlab = "", xaxt = "n", main = "", las = 2, cex.axis = 1.2, ylab = "")

#axis(1, cex.axis = 1.2)
#mtext(1, text = "Logistic regression scores (0 = flav, 1 = form)", line = 2.5, cex = 1.2)
#mtext(2, text = "Density", line = 3.5, cex = 1.2)


#lines(density(pred.pops), col = "firebrick1", lwd = 4)

#legend('topleft', 
#       legend = c("Herbarium specimens", "Field collections 2017-2018"), 
#       col = c("dodgerblue", "firebrick1"), lwd = 4, 
#       bty = "n")


################################
################################
#Stripchart of size-corrected hybrid index values
#names(pred.pops) <- as.numeric(as.factor(coll.dat$site))

library(RColorBrewer)
colvec7 <- brewer.pal(7, "Dark2")
set.seed(1)
colvec7.1 <- sample(brewer.pal(7, "Set2")) #good color scheme, don't rerun

coll.dat$pred.pops <- size.correct(coll.dat, scaling = "sigmoid")

coll.dat$site <- factor(coll.dat$site, 
                    levels = levels(as.factor(coll.dat$site) )[order( tapply( 
                      coll.dat$pred.pops, coll.dat$site, mean) ) ] )


par(mar =c(8,5,2,1))
stripchart(pred.pops ~ site, data = coll.dat, vertical = TRUE, 
           col = colvec7.1[ c(4, 3, 1, 5, 2) ], pch = 16, 
           xaxt = "n", las = 2, ylab = "Logistic Regression fitted value")
axis(1, labels = levels(coll.dat$site), at = 1:6, las = 2)





####################### LDA stripchart

p18 <- read.csv("/Users/jeff/Projects/Aquilegia/porcupine_transect_phenotypes.csv", stringsAsFactors = FALSE)
p18.species <- vector()
adj <- grep("adj", p18$site)
p18.species[adj] <- "formosa"
p18.species[-adj] <- "hybrid"
p18$species <- p18.species

cols <- intersect(colnames(coll.dat), colnames(p18))
new.dat <- rbind.data.frame(coll.dat[, cols], p18[ ,cols] )

new.dat$site <- as.character(new.dat$site)
new.dat$site[grep("porcupine", new.dat$site)] <- "Porcupine Ridge"
new.dat$site[grep("adj", new.dat$site)] <- "Adjacent valley"
new.dat$site <- as.factor(new.dat$site)
levels(new.dat$site)

#new.dat <- new.dat[-grep("Adj", new.dat$site), ]
new.dat$site <- droplevels(new.dat$site)

new.dat$pred.pops <- size.correct(new.dat, scaling = "sigmoid")
#minv <- min(pred.pops)
#maxv <- max(pred.pops)
#new.dat$pred.pops <- (pred.pops - minv)/(maxv - minv)

names(new.dat$pred.pops) <- as.numeric(as.factor(new.dat$site))




library(RColorBrewer)
#colvec7 <- brewer.pal(7, "Dark2")
#set.seed(1)
#colvec7.1 <- sample(brewer.pal(7, "Set2")) #good color scheme, don't rerun
colvec8 <- brewer.pal(8, "Dark2")
colvec8.1 <- colvec8[c(6, 2,3,4,5, 8, 7, 1)]

new.dat$site <- factor(new.dat$site, 
                        levels = levels(as.factor(new.dat$site) )[order( tapply( 
                          pred.pops, new.dat$site, mean) ) ] )

par(mar = par()$mar + c(0,1,0,0), bty = "o")
stripchart(pred.pops ~ site, data = new.dat, vertical = TRUE, 
           col = colvec8.1, pch = 1, 
           xaxt = "n", las = 2, yaxt = "n", ylab = "Size-adjusted discriminant score",
           method = "jitter", lwd = 2, cex.lab = 1.2)

#axis(1, labels = levels(new.dat$site), at = 1:8, las = 2)
axis(2, lwd = 0, lwd.ticks = 1, las = 2, cex = 1.2)
#axis(1, at = 1:8, labels = c("MK", "MR", "CM", "PV", "AV", "MP", "CW", "RL"))
legend("bottomright", levels(new.dat$site), col = colvec8.1, pch = 1, pt.lwd = 2, bty = "n")


######################

# morphology- altitude cline investigation -----------------------------------------

p18 <- read.csv("/Users/jeff/Projects/Aquilegia/porcupine_transect_phenotypes.csv", stringsAsFactors = FALSE)

ind.val1 <- grep("porcupine", p18$site)
ind.val2 <- grep("adj", p18$site)

#prediction according to LDA from herbarium training set
trsc.pred1 <- predict(LDA.trn, p18[,c(4:11)], type = "response", prior = c(.5, .5))$x
    #both valleys
     plot(trsc.pred1 ~ p18$elevation)
     #single valley
     plot(trsc.pred1[ind.val1] ~ p18$elevation[ind.val1], ylim = c(-3,3))

#prediction according to LDA from all herbarium data
trsc.pred2 <- predict(LDA.all, p18[,c(4:11)], type = "response", prior = c(.5, .5))$x
    #both valleys
    plot(trsc.pred2 ~ p18$elevation)
    #single valley
    plot(trsc.pred2[ind.val1] ~ p18$elevation[ind.val1], ylim = c(-2,2))

#prediction according to logistic regression from all herbarium data, don't include presence/absence trait
trsc.pred3 <- predict(LOGREG.all, newdata = p18, type = "response")
plot(trsc.pred3 ~ p18$elevation)


 tr.pred2 <- predict(z2, p18[,c(4:6, 8)], type = "response")
tr.pred3 <- predict(z3, p18[,c(4:10)], type = "response")


pch18 <- p18$site
adj <- grep("adj", pch18)
pch18[adj] <- 15
pch18[-adj] <- 16
pch18 <- as.numeric(pch18)


plot(tr.pred1 ~ p18$elevation, ylim = c(0,1), pch = pch18, xlim = c(1400, 2200),
     ylab = "Logistic Regression fitted value", xlab = "elevation (m)")



plot(tr.pred2 ~ p18$elevation, ylim = c(0,1), pch = pch18, xlim = c(1400, 2200),
     ylab = "Logistic Regression fitted value", xlab = "elevation (m)")
cline <- lm(tr.pred2 ~ p18$elevation)
abline(cline)

plot(tr.pred3 ~ p18$elevation, ylim = c(0,1), pch = pch18, xlim = c(1400, 2200),
     ylab = "Logistic Regression fitted value", xlab = "elevation (m)")
cline.3 <- lm(tr.pred3 ~ p18$elevation)
abline(cline.3)


pheno <- predict(z.herb, p18[, -c(1:3,12)], prior = c(.5,.5))$x




pch18 <- p18$site
adj <- grep("adj", pch18)
pch18[adj] <- 15
pch18[-adj] <- 16
pch18 <- as.numeric(pch18)

plot(pheno ~ p18$elevation, ylim = range(sc), pch = pch18, xlim = c(1400, 2200),
     ylab = "flavescens - formosa discriminant axis", xlab = "elevation (m)")

plot(pheno ~ p18$elevation, ylim = range(-4,4), pch = pch18, xlim = c(1400, 2200),
     ylab = "flavescens - formosa discriminant axis", xlab = "elevation (m)")



#run model for specimens in both valleys
df18 <- cbind.data.frame(elevation = p18$elevation, pheno)
z18 <- lm(pheno ~ elevation, data = df18)

xs <- range(df18$elevation) #endpoints for regression line
ys <- predict(z18, newdata = data.frame(elevation = xs))

newx18 <- seq(min(df18$elevation, na.rm = TRUE), max(df18$elevation, na.rm = TRUE), length.out=100)
py18 <- predict(z18, newdata = data.frame(elevation = newx18),  interval = 'confidence')

#plot regression line and confidence bands
polygon(c(rev(newx18), newx18), c(rev(py18[ ,3]), py18[ ,2]), col = rgb(190/255,190/255,190/255, .6), border = NA)
lines(y = ys, x = xs)


#subset to include only Porcupine valley
pv18 <- subset(p18, site == "porcupine valley transect")
pv.index <- which(p18$site == "porcupine valley transect")
plot(tr.pred2[pv.index] ~ p18$elevation[pv.index], pch = 16, 
     xlim = c(1800, 2150), ylim = c(0,1), 
     xlab = "Elevation (m)", ylab = "Flavescens - formosa discriminant axis")

#run model for this valley subset
pv18 <- subset(p18, site == "porcupine valley transect")
pvdf18 <- cbind.data.frame(elevation = p18$elevation[pv.index], pheno = tr.pred2[pv.index])
pvz18 <- lm(pheno ~ elevation, data = pvdf18)

#regression line and confidence bands
newxpv <- range(pvdf18$elevation)
newypv <- predict(pvz18, newdata = data.frame(elevation = newxpv))

xconf <-  seq(min(pvdf18$elevation, na.rm = TRUE), max(pvdf18$elevation, na.rm = TRUE), length.out=100)
yconf <- predict(pvz18, newdata = data.frame(elevation = xconf),  interval = 'confidence')
polygon(c(rev(xconf), xconf), c(rev(yconf[ ,3]), yconf[ ,2]), col = rgb(190/255,190/255,190/255, .6), border = NA)
lines(x = newxpv, y = newypv)

anova(pvz18)









plot.new()
plot.window(xlim = c(0,1), ylim = c(0,20))
tapply(pred.pops, coll.dat$site, FUN = function(x){ 
  lines(density(x)$y ~ density(x)$x,
        col =  colvec7.1[ mean(as.numeric(names(x))) ], 
        lwd = 4) })












plot.new()
par(mar = c(5,6,4,2))
plot(density(herb.scores), ylim = c(0,0.3), col = "dodgerblue", lwd = 4, 
     xlab = "", xaxt = "n", main = "", las = 2, cex.axis = 1.2, ylab = "")


predict.pops <- predict(LDA.herb, newdata = coll.dat[,5:12])$x #phenotype ordination scores along formosa-flavescens axis

plot.new()
par(mar = c(5,6,4,2))
plot(density(herb.scores), ylim = c(0,0.3), col = "dodgerblue", lwd = 4, 
     xlab = "", xaxt = "n", main = "", las = 2, cex.axis = 1.2, ylab = "")

axis(1, cex.axis = 1.2)
mtext(1, text = "flavescens - formosa discriminant axis", line = 2.5, cex = 1.2)
mtext(2, text = "Density", line = 3.5, cex = 1.2)

lines(density(predict.pops), col = "firebrick1", lwd = 4, add = TRUE)

legend('topleft', 
       legend = c("Herbarium specimens", "Field collections 2017-2018"), 
       col = c("dodgerblue", "firebrick1"), lwd = 4, 
       bty = "n")




plot.new()
plot.window(xlim = range(LDA.all.scrs), ylim = c(0,1))

predict.pops <- size.correct(coll.dat)
names(predict.pops) <- as.numeric(as.factor(coll.dat$site))

library(RColorBrewer)
colvec7 <- brewer.pal(7, "Dark2")
set.seed(1)
colvec7.1 <- sample(brewer.pal(7, "Set2")) #good color scheme, don't rerun

tapply(predict.pops, coll.dat$site, FUN = function(x){ 
                                      lines(density(x)$y ~ density(x)$x,
                                            col =  colvec7.1[ mean(as.numeric(names(x))) ], 
                                            lwd = 4) })

axis(1, cex.axis = 1.2, lwd = 0, lwd.ticks = 1)
axis(2, las = 2, cex.axis = 1.2, lwd = 0, lwd.ticks = 1)
mtext(1, text = "flavescens - formosa discriminant axis", cex = 1.2, line = 2.5)
mtext(2, text = "Density", cex = 1.2 , line = 3)
box()


legend('topleft', legend = levels(as.factor(coll.dat$site)), 
       col = colvec, lwd = 4, bty = "n")



# RDA axis ----------------------------------------------------------------

herb.rda.scrs <- scores(RDA.herb, display = c("sites", "species"), scaling = 2)


rda.axis <- predict(RDA.herb, type="wa", model="CCA", newdata = coll.dat[, 5:12],   scaling = 2 )
pca.axis <- predict(RDA.herb, type="wa", model="CA", newdata = coll.dat[, 5:12],  scaling = 2 )[, 1]
pca2.axis <- predict(RDA.herb, type="wa", model="CA", newdata = coll.dat[, 5:12],  scaling = 2 )[, 2]

dev.off()
plot.new()
par(mar = c(5,5,2,1), oma = c(5,0,0,0), xpd = NA)
plot.window(xlim = range( herb.rda.scrs$sites[,1]) + c(0,2), ylim = range(herb.rda.scrs$sites[,2]) ) 

colvec7.1[as.factor(coll.dat$site)]


dataEllipse(x = as.numeric(rda.axis), y = pca.axis, levels = 0.95, group.labels = rep("",6),
            groups = as.factor(coll.dat$site), grid = FALSE,
            pch = rep(16,6), cex = 0, center.cex = 0,
            col = colvec7.1,
            xlim = range( herb.rda.scrs$sites[,1]), ylim = range(herb.rda.scrs$sites[,2]), 
            xlab = "", ylab = "", axes = FALSE)
points(pca.axis ~ rda.axis, pch = 21, bg = colvec7.1[as.factor(coll.dat$site)])

axis(1, cex.axis = 1.2, lwd = 0, lwd.ticks = 1)
axis(2, cex.axis = 1.2, lwd = 0, lwd.ticks = 1, las = 2)
box()

mtext(1, text = "flavescens - formosa canonical RDA axis", cex = 1.1, line = 2.5)
mtext(2, text = "PCA1 of residual variation", cex = 1.1, line = 3)

legend('bottomright', legend = levels(as.factor(coll.dat$site)) [1:3], 
       col = colvec7.1[1:3], lwd = 4, bty = "n", inset = c(.55,-.6))
legend('bottomright', legend = levels(as.factor(coll.dat$site)) [4:6], 
       col = colvec7.1[4:6], lwd = 4, bty = "n", inset = c(.1,-.6))

#try plot of 1st two non-canonical PCA axes
dev.off()
par(mar = c(5,5,2,1), oma = c(5,0,0,0), xpd = NA)
dataEllipse(x = as.numeric(pca.axis), y = pca2.axis, levels = 0.95, group.labels = rep("",6),
            groups = as.factor(coll.dat$site), grid = FALSE,
            pch = rep(16,6), cex = 0, center.cex = 0,
            col = colvec7.1,
            xlim = range(pca.axis) + c(-.5,.5), ylim = range(pca2.axis) + c(-.5,.5), 
            xlab = "", ylab = "", axes = FALSE)


points(pca2.axis ~ pca.axis, pch = 21, bg = colvec7.1[as.factor(coll.dat$site)])

axis(1, cex.axis = 1.2, lwd = 0, lwd.ticks = 1)
axis(2, cex.axis = 1.2, lwd = 0, lwd.ticks = 1, las = 2)
box()

mtext(1, text = "Unconstrained PCA1", cex = 1.1, line = 2.5)
mtext(2, text = "Unconstrained PCA2", cex = 1.1, line = 3)

legend('bottomright', legend = levels(as.factor(coll.dat$site)) [1:3], 
       col = colvec7.1[1:3], lwd = 4, bty = "n", inset = c(.55,-.6))
legend('bottomright', legend = levels(as.factor(coll.dat$site)) [4:6], 
       col = colvec7.1[4:6], lwd = 4, bty = "n", inset = c(.1,-.6))




# Color phenotypes --------------------------------------------------------

color <- read.csv("/Users/jeff/Projects/Aquilegia/raw_data/color.csv", stringsAsFactors = FALSE)
color.sub <- color[-which(color$)]
color$log.rg <-  log( color$red.mean / color$green.mean) 

plot.new()

plot(density(color$log.rg))


?#hy <- which(d$species == "hybrid") #species indices
#fo <- which(d$species == "A. formosa")
#fl <- which(d$species == "A. flavescens")

#prnt <- d[-hy,] #data subset without marble range phenotypes
#prnt$species <- droplevels(prnt$species) #get rid of extra factor level
#prnt <- prnt[-8,] #MK08 which had pink flower

#fo.p <- which(prnt$species == "A. formosa") #species indices
#fl.p <- which(prnt$species == "A. flavescens") 

#mk <- which(prnt$site == "Mt. Kobau") #population indices
#cl <- which(prnt$site == "Clearwater")
#rl <- which(prnt$site == "Roberts Lake")

predict.pops <- predict(LDA.herb, newdata = coll.dat[,5:12])

predict.pops

fo.traits <- d[fo,4:10]
fl.traits <- d[fl,4:10]
hy.traits <- d[hy,4:10]




# source LDA axis from herbarium specimen analysis ---------------------------------------------------------------------



mk8 <- predict(z, d[8, -c(1,3,12:19)], prior = c(.5,.5))

z$scaling #gives the coefficients of the linear combination of traits

#make figure
par(bty = "l", mar = c(5,5,4,2))
plot.new()
sc <- predict(LDA.h)$x

#Plot RL
stripchart(sc[rl], vertical = FALSE, at = .65, ylim = c(0,1), xlim = range(sc) + c(0.1,0.3), 
           pch = 21, bg = "#F50000", col = "black", method = "jitter", axes= FALSE)
#Plot WG
stripchart(sc[cl], vertical = FALSE, at = .65, ylim = c(0,1), xlim = range(sc) + c(0.1,0.3), 
           pch = 23, bg = "#F50000", col = "black", method = "jitter", add = TRUE)
#Plot MK
stripchart(sc[mk], vertical = FALSE, at = .65, ylim = c(0,1), xlim = range(sc) + c(0.1,0.3), 
           pch = 24, bg = "#F5FF0A", col = "black", method = "jitter", add = TRUE)

#add densities
x.fla <- seq(min(sc[fl.p]), to=max(sc[fl.p]), length.out = 512)
y.fla <- density(sc[fl.p])$y
polygon(x.fla, y.fla, col = rgb(245/255,1,10/255,0.6), lwd = 2)

x.fo <- seq(min(sc[fo.p]), to=max(sc[fo.p]), length.out = 512)
y.fo <- density(sc[fo.p])$y
polygon(x.fo, y.fo, col = rgb(1,0,0,0.6), lwd = 2)

#Predict hybrids and MK08
mk8 <- predict(z, d[8, -c(1,3,12:19)], prior = c(.5,.5))
hyb <- predict(z, d[hy, -c(1,3,12:19)], prior = c(.5,.5))
hyb.lda.fl <- hyb$x
names(hyb.lda.fl) <- d[hy,1]

#plot MK8

#points(x = rep(mk8$x,3), y = rep(.75, 3), pch = c(24,25,11),
#       bg = c("#F5FF0A","#F5FF0A","black"), col = c("#F5FF0A","#F5FF0A","black"))

points(x = mk8$x, y = .75, pch = 24,
       bg = "#F5FF0A")

#plot hybrids 
points(x = hyb$x, y = jitter(rep(0.65,32), amount = .1), pch = 4, 
       lwd = 2, col = "#ff7880")

hy.coords <- hyb$x
x.hy <- seq(min(hy.coords), to=max(hy.coords),
            length.out = 512)
y.hy <- density(hy.coords)$y
polygon(x.hy, y.hy, col = rgb(1,120/255,128/255,0.5), lwd = 2)

axis(2, at = c(0,.1,.2,.3,.4,.5), labels = TRUE, las = 2, cex.axis = 1.2)
mtext(text = "Density", side = 2, line = 3, adj =.25, cex = 1.2)
axis(1, line = -.5, cex.axis = 1.2)
#axis(1, line = -9)

mtext(text = "Floral morphology discriminant axis", side = 1, line =2, cex = 1.2)
#mtext(text = "LDA1", side = 1, line =-7)






######## analyze phenotypes of 2018 populations ############
d18 <- read.csv("/Users/jeff/Projects/Aquilegia/collections.2018.csv", stringsAsFactors = FALSE)

#x <- predict(z, d18[, -c(1,2)], prior = c(.5,.5))

#points(x = x$x, y = jitter(rep(.65,22), amount = .1), pch = 4, 
       lwd = 2)

#combine mission ridge with mt. kobau, then predict hybrids
#dat <- d[,-c(12:18)]
#all <- rbind(dat, mr)
#all$site[which(all$site == "mission.ridge" | all$site == "")] <- "Mission Ridge"



########try pca
d18$site
d$site

d.17.18 <- merge(d18, d, all = T)


pca <- prcomp(d.17.18[,c(4:11)], scale = TRUE)

pc1 <- pca$x[,1]
pc2 <- pca$x[,2]

pchz <- vector()
pchz[which(d.17.18$site == "Roberts Lake")] <- 16
pchz[which(d.17.18$site == "Clearwater")] <- 18
pchz[which(d.17.18$site == "Mt. Kobau")] <- 17
pchz[which(d.17.18$site == "Mission Ridge")] <- 15
pchz[which(d.17.18$site == "Porcupine Peak")] <- 4
pchz[which(d.17.18$site == "Rogers Pass")] <- 6


bgz <- pchz
bgz[which(pchz == 21 | pchz == 23)] <- "#F50000"
bgz[which(pchz == 24)] <- "#F5FF0A"
bgz[which(pchz == 4)] <- "#ff7880"
bgz[which(pchz == 15)] <- "black"

clz <- pchz
clz[which(clz == 21 | clz == 23 | clz == 24 | clz == 15)] <- "black"
clz[which(clz == "4")] <- "#ff7880"


plot(pc2 ~ pc1, pch = pchz,  xlab = "PC1 (47.4% of variance)", ylab = "PC2 (22.3% of variance)")


points(pc2[hy] ~ pc1[hy], pch = 4, col = "#ff7880", lwd = 2)
summary(pca)




mrb <- which(all$species == "hybrid")
prns <- all[-mrb,]
prns$species <- as.character(prns$species)
prns$species[which(is.na(prns$species))] <- "A. flavescens" #change mission ridge to A. flavescens
prns$species

zz <- lda(species ~ ., data = prns[,-c(1,3,12:19)], prior = c(.5,.5))


par(bty = "l", mar = c(5,5,4,2))
plot.new()
sco <- predict(zz)$x

#create indices for plotting
prns$site[which(prns$site == "" | prns$site == "mission.ridge")] <- "Mission Ridge"
msr <- which(prns$site == "Mission Ridge")
rbs <- which(prns$site == "Roberts Lake")
clr <- which(prns$site == "Clearwater")
mtk <- which(prns$site == "Mt. Kobau")

#Plot RL
stripchart(sco[rbs], vertical = FALSE, at = .65, ylim = c(0,1), xlim = range(sc) + c(0.1,0.3), 
           pch = 21, bg = "#F50000", col = "black", method = "jitter", axes= FALSE)
#Plot WG
stripchart(sco[clr], vertical = FALSE, at = .65, ylim = c(0,1), xlim = range(sc) + c(0.1,0.3), 
           pch = 23, bg = "#F50000", col = "black", method = "jitter", add = TRUE)
#Plot MK
stripchart(sco[mtk], vertical = FALSE, at = .65, ylim = c(0,1), xlim = range(sc) + c(0.1,0.3), 
           pch = 24, bg = "#F5FF0A", col = "black", method = "jitter", add = TRUE)
#Plot Mission ridge
stripchart(sco[msr], vertical = FALSE, at = .65, ylim = c(0,1), xlim = range(sc) + c(0.1,0.3), 
           pch = 16, method = "jitter", add = TRUE)

#now predict hybrids
hybz <- predict(zz, d[hy, -c(1,3,12:19)], prior = c(.5,.5))
points(x = hybz$x, y = jitter(rep(0.65,32), amount = .1), pch = 4, 
       lwd = 2, col = "#ff7880")


#add densities
x.fla <- seq(min(sc[fl.p]), to=max(sc[fl.p]), length.out = 512)
y.fla <- density(sc[fl.p])$y
polygon(x.fla, y.fla, col = rgb(245/255,1,10/255,0.6), lwd = 2)

x.fo <- seq(min(sc[fo.p]), to=max(sc[fo.p]), length.out = 512)
y.fo <- density(sc[fo.p])$y
polygon(x.fo, y.fo, col = rgb(1,0,0,0.6), lwd = 2)


####################
shits <- all
shits$species <- as.character(shits$species)
shits$site[which(shits$site == "" | shits$site == "mission.ridge")] <- "Mission Ridge"
shits$species[hy] <- "A. flavescens"
mrd <- which(shits$site == "Mission Ridge")
shits$species[mrd] <- "predict"

prz <- shits[-mrd, ]
zzz <- lda(species ~ ., data = prz[,-c(1,3,12:19)], prior = c(.5,.5))

scoz <- predict(zzz)$x

#create indices for plotting
msrz <- which(prz$site == "Mission Ridge")
rbsz <- which(prz$site == "Roberts Lake")
clrz <- which(prz$site == "Clearwater")
mtkz <- which(prz$site == "Mt. Kobau")
prcz <- which(prz$site == "Porcupine Peak")

#Plot RL
stripchart(scoz[rbsz], vertical = FALSE, at = .65, ylim = c(0,1), xlim = range(sc) + c(0.1,0.3), 
           pch = 21, bg = "#F50000", col = "black", method = "jitter", axes= FALSE)
#Plot WG
stripchart(scoz[clrz], vertical = FALSE, at = .65, ylim = c(0,1), xlim = range(sc) + c(0.1,0.3), 
           pch = 23, bg = "#F50000", col = "black", method = "jitter", add = TRUE)
#Plot MK
stripchart(scoz[mtkz], vertical = FALSE, at = .65, ylim = c(0,1), xlim = range(sc) + c(0.1,0.3), 
           pch = 24, bg = "#F5FF0A", col = "black", method = "jitter", add = TRUE)
#Plot porcupine peak
stripchart(scoz[prcz], vertical = FALSE, at = .65, ylim = c(0,1), xlim = range(sc) + c(0.1,0.3), 
           pch = 16, method = "jitter", add = TRUE)


#Plot Mission ridge
stripchart(scoz[msrz], vertical = FALSE, at = .65, ylim = c(0,1), xlim = range(sc) + c(0.1,0.3), 
           pch = 16, method = "jitter", add = TRUE)

#now predict mission ridge
mrdz <- predict(zzz, shits[mrd, -c(1,3,12:19)], prior = c(.5,.5))
points(x = mrdz$x, y = jitter(rep(0.65,22), amount = .1), pch = 1)


########try pca

pca <- prcomp(all[,-c(1:3, 11:19)], scale = TRUE)

pc1 <- pca$x[,1]
pc2 <- pca$x[,2]

pchz <- vector()
pchz[which(all$site == "Roberts Lake")] <- 21
pchz[which(all$site == "Clearwater")] <- 23
pchz[which(all$site == "Mt. Kobau")] <- 24
pchz[which(all$site == "Mission Ridge")] <- 15
pchz[which(all$site == "Porcupine Peak")] <- 4

     
bgz <- pchz
bgz[which(pchz == 21 | pchz == 23)] <- "#F50000"
bgz[which(pchz == 24)] <- "#F5FF0A"
bgz[which(pchz == 4)] <- "#ff7880"
bgz[which(pchz == 15)] <- "black"

clz <- pchz
clz[which(clz == 21 | clz == 23 | clz == 24 | clz == 15)] <- "black"
clz[which(clz == "4")] <- "#ff7880"


plot(pc2 ~ pc1, pch = pchz, bg = bgz, col = clz, xlab = "PC1 (57.7% of variance)", ylab = "PC2 (23.9% of variance)")


points(pc2[hy] ~ pc1[hy], pch = 4, col = "#ff7880", lwd = 2)
summary(pca)



###########################################################################
###########################################################################
###########################################################################
################## Predict Mt Tatlow
t <- read.csv("/Users/jeff/Projects/Aquilegia/raw_data/mt.tatlow.csv")
head(t)
tatlow <- predict(z, t)

tatlow.pts <- tatlow$x

points(x = tatlow.pts, y = jitter(rep(0.65,6), amount = .2), pch = 4, 
       lwd = 2, col = "#FF3399")






###########################################################################
###########################################################################
# V-FOLD CROSS VALIDATION -------------------------------------------------
###########################################################################
###########################################################################
vlda = function(v,formula,data,cl){
  require(MASS)
  grps = cut(1:nrow(data),v,labels=FALSE)[sample(1:nrow(data))]
  pred = lapply(1:v,function(i,formula,data){
    omit = which(grps == i)
    z = lda(formula,data=data[-omit,])
    predict(z,data[omit,])
  },formula,data)
  
  wh = unlist(lapply(pred,function(pp)pp$class))
  table(wh,cl[order(grps)]) # I don't understand how this last part works but it does!
}

############################ error rate ###################################
mis <- vector()
for(i in 1:1000){
  xval <- vlda(5,species~., data=prnt[,-c(1,3,12:19)], cl=prnt$species)
  mis[i] <- sum(xval[row(xval) != col(xval)]) / sum(xval)
}
mean(mis) # 0.0001688312 error rate
###########################################################################

###########################################################################
############Rerun with smaller sample size (= to that in genetic analysis)#
###########################################################################

prnt.sub <- rbind(prnt[sample(1:19,15),], prnt[sample(20:77, 18),])
z.sub <- lda(species ~ ., data = prnt.sub[,-c(1,3,12:19)], prior = c(.5,.5))
clust.sub <- prnt.sub$species

mis.sub <- vector()
for(i in 1:1000){
  xval <- vlda(5,species~., data=prnt.sub[,-c(1,3,12:19)], cl=prnt.sub$species)
  mis.sub[i] <- sum(xval[row(xval) != col(xval)]) / sum(xval)
}
mean(mis.sub) # 0.007575758












###########################################################################
###########################################################################
# correlation between morph and genetic LDAs --------------------------------------------------------
###########################################################################
names(hyb$x) <- d[hy, "id"] #hyb is from morphometric analysis

gen.lda <- pred$ind.scores[-1] #genetic lda scores of hybrids
morph.lda <- hyb$x[c("PP_01", "PP_02", "PP_03", "PP_05", "PP_07", "PP_08", 
                   "PP_09", "PP_10")]
plot(morph.lda ~ gen.lda)
###########################################################################
###########################################################################





##########################################################################
#calculate distance from spurs to anthers
#spur length + blade length + anther exsertion
##########################################################################

d$dist <- d$spur.length + d$blade.length + d$anther.exsertion

boxplot(d$dist[fo], at = .5, xlim =c(0,3),ylim=c(1.5,4.5), 
        horizontal = T,col = "orangered2",pch = 16)

boxplot(d$dist[fl], at = 2.5,horizontal=T,
        col = "gold1", pch = 16, add = TRUE)

points(x = 1.6, y = 2)

text(c(1.6,2), "Calliope")

##########################################################################
##########################################################################
#estimate proportion of indiv. with celft laminae#########################
##########################################################################
##########################################################################

bl <- d$blade.cleft
fo.cleft <- sum(bl[fo]=="y")
fo.n <- length(bl[fo])
p.fo <- fo.cleft/fo.n
p.fo # 0.3793103
fo.se <- sqrt(p.fo*(1-p.fo)/fo.n)
fo.se #0.06371191
binom.agresti.coull(fo.cleft,fo.n)

fl.cleft <- sum(bl[fl]=="y")
fl.n <- length(bl[fl])

#for Marble Range flowers
hy.cleft <- sum(bl[hy]=="y")
hy.n <- length(bl[hy])
binom.agresti.coull(hy.cleft,hy.n)
# p= 0.21875, 95% CI 0.10734431 < x < 0.3904452

#z test for difference in proportions
#read in herbarium data to use proportion pooled from 
#collections and herbarium specimens
x <- c(fo.cleft + cleft.herb, hy.cleft)
n <- c(fo.n + n.herb, hy.n)
z <- prop.test(x,n)
z

##########################################################################
##########################################################################
# RDA of parental populations --------------------------------------------
##########################################################################
##########################################################################

expl <- model.matrix(~d$species)
resp <- d[,c(2,4:10),]
colnames(resp) <- c("species","corolla width","spur length","blade length","blade width","anther exsertion","sepal length","sepal width")

#run RDA only on parental types
RDA <- rda(resp[1:78,-1] ~expl[1:78,],scaling=1)
s <- summary(RDA) #57.21 % of cumulative variance
s
weights <- s$species[,1] #weightings of traits onto first axis
ind <- order(abs(weights),decreasing=TRUE)
weights[ind]
##########################################################################


##########################################################################
##########################################################################
# summary stats ----------------------------------------------------------
##########################################################################
# make trait table -------------------------------------------------------
# include all parental types
# use ind to sort traits in decreasing order of RDA1 loadings
##########################################################################

rnames <- vector()
rnames[seq(1,21,3)] <- names(weights[ind])
rnames[seq(2,20,3)] <- seq(2,20,3) #rownames can't be duplicates
rnames[seq(3,21,3)] <- seq(3,21,3)
v <- data.frame(matrix(NA,nrow=21,ncol=4),row.names=rnames) #dimensions of desired table
colnames(v) <- c("mean","se","95% CI lower","95% CI upper")
sp <- c("A.fo", "A.fla", "hy")
v <- cbind(rep(sp,7),v)
v

#add mean to matrix
for(i in 1:7){
  #for A. formosa
  v[3*i-2,2] <- round(mean(fo.traits[,ind[i]],na.rm=TRUE),2)
  #for A. flavescens
  v[3*i-1,2] <- round(mean(fl.traits[,ind[i]],na.rm=TRUE),2)
  #for hybrids
  v[3*i,2] <- round(mean(hy.traits[,ind[i]],na.rm=TRUE),2)
}

#create length vectors
l.fo <- rep(0,7)
l.fl <- l.fo
l.hy <- l.fo
for(i in 1:7){
  l.fo[i] <- length( fo.traits[,ind[i]])
  l.fl[i] <- length(fl.traits[,ind[i]])
  l.hy[i] <- length(hy.traits[,ind[i]])
}
l.fo #58
l.fl #20
l.hy #32

#add s.e. to matrix
for(i in 1:7){
  #for A. formosa
  v[3*i-2,3] <- round(sd(fo.traits[,ind[i]],na.rm=TRUE)/sqrt(l.fo[i]),2)
  #for A. flavescens
  v[3*i-1,3] <- round(sd(fl.traits[,ind[i]],na.rm=TRUE)/sqrt(l.fl[i]),2)
  #for hybrids
  v[3*i,3] <- round(sd(hy.traits[,ind[i]],na.rm=TRUE)/sqrt(l.hy[i]),2)
}

#create vector of critical values
t.fo <- rep(0,7)
t.fl <- t.fo 
t.hy <- t.fo

for(i in 1:7){
  t.fo[i] <- abs(qt(0.025, df=l.fo[i]-1))
  t.fl[i] <- abs(qt(0.025, df=l.fl[i]-1))
  t.hy[i] <- abs(qt(0.025, df=l.hy[i]-1))
}
t.fo 
t.fl 
t.hy #t vectors already in decreasing order for trait RDA1 loadings
#higher values for A. flavescens = lower sample size

#add confidence intervals
for(i in 1:7){
  #for A. formosa
  v[3*i-2,4] <- round(v[3*i-2,2] - t.fo[i]* v[3*i-2,3],2)
  v[3*i-2,5] <- round(v[3*i-2,2] + t.fo[i]* v[3*i-2,3],2)
  #for A. flavescens
  v[3*i-1,4] <- round(v[3*i-1,2] - t.fl[i]*v[3*i-1,3],2)
  v[3*i-1,5] <- round(v[3*i-1,2] + t.fl[i]*v[3*i-1,3],2)
  #for hybrids
  v[3*i,4] <- round(v[3*i,2] - t.fl[i]*v[3*i,3],2)
  v[3*i,5] <- round(v[3*i,2] + t.fl[i]*v[3*i,3],2)
}
v
##########################################################################
##########################################################################
# multipanel plot --------------------------------------------------------
##########################################################################
##########################################################################

# 7 figures arranged in 4 rows and 2 columns

#pdf("figures/multipanel.traits.pdf")
par(mfrow=c(4,2),bty="l",mar=c(2,5,2,2),oma=c(1,0,0,0))
cols <- c("#F5FF0A","#ff7880","#F50000")
#topleft:anther exertion
set.seed(123)
stripchart(anther.exsertion~species,data=d,method="jitter",vertical=TRUE,
            pch=1, ylab="Anther exsertion (cm)",xaxt="n",las=2, cex.axis = 1.2, cex.lab= 1.2)
set.seed(123)
stripchart(anther.exsertion~species,data=d,method="jitter",vertical=TRUE,
           add=TRUE, pch=16,axes=FALSE,col=cols,las=2)
points(1:3 +.2, tapply(d$anther.exsertion,d$species,mean),pch=20)

#add CI
lines(c(1,1) +0.2, c(v[2,4],v[2,5])) 
lines(c(1.17,1.23) ,y= c(v[2,4],v[2,4]))
lines(c(1.17,1.23) ,y= c(v[2,5],v[2,5]))

lines(c(2,2) +0.2, c(v[3,4],v[3,5])) 
lines(c(2.17,2.23) ,y= c(v[3,4],v[3,4]))
lines(c(2.17,2.23) ,y= c(v[3,5],v[3,5]))

lines(c(3,3)+.2, c(v[1,4],v[1,5]))
lines(c(3.17,3.23) ,y= c(v[1,4],v[1,4]))
lines(c(3.17,3.23) ,y= c(v[1,5],v[1,5]))

#topright: corolla width
set.seed(123)
stripchart(corolla.width~species,data=d,method="jitter",vertical=TRUE,
           pch=1,ylab="Corolla width (cm)",xaxt="n",las=2, cex.axis = 1.2, cex.lab = 1.2)
set.seed(123)
stripchart(corolla.width~species,data=d,method="jitter",vertical=TRUE,
           add=TRUE,pch=16,axes=FALSE,col=cols,las=2)
points(1:3+.2, tapply(d$corolla.width,d$species, mean),pch=20)

#add CI
lines(c(1,1) +0.2, c(v[17,4],v[17,5])) 
lines(c(1.17,1.23) ,y= c(v[17,4],v[17,4]))
lines(c(1.17,1.23) ,y= c(v[17,5],v[17,5]))

lines(c(2,2) +0.2, c(v[18,4],v[18,5])) 
lines(c(2.17,2.23) ,y= c(v[18,4],v[18,4]))
lines(c(2.17,2.23) ,y= c(v[18,5],v[18,5]))

lines(c(3,3)+.2, c(v[16,4],v[16,5]))
lines(c(3.17,3.23) ,y= c(v[16,4],v[16,4]))
lines(c(3.17,3.23) ,y= c(v[16,5],v[16,5]))

#mid1 left: spur length
set.seed(123)
stripchart(spur.length~species,data=d,method="jitter",vertical=TRUE,
           pch=1,ylab="Spur length (cm)",xaxt="n",las=2, cex.axis = 1.2, cex.lab = 1.2)
set.seed(123)
stripchart(spur.length~species,data=d,method="jitter",vertical=TRUE,
           add=TRUE,pch=16,axes=FALSE,col=cols,las=2)
points(1:3+.2, tapply(d$spur.length,d$species, mean),pch=20)

lines(c(1,1) +0.2, c(v[5,4],v[5,5])) 
lines(c(1.17,1.23) ,y= c(v[5,4],v[5,4]))
lines(c(1.17,1.23) ,y= c(v[5,5],v[5,5]))

lines(c(2,2) +0.2, c(v[6,4],v[6,5])) 
lines(c(2.17,2.23) ,y= c(v[6,4],v[6,4]))
lines(c(2.17,2.23) ,y= c(v[6,5],v[6,5]))

lines(c(3,3)+.2, c(v[4,4],v[4,5]))
lines(c(3.17,3.23) ,y= c(v[4,4],v[4,4]))
lines(c(3.17,3.23) ,y= c(v[4,5],v[4,5]))

#mid1 right:blade length
set.seed(123)
stripchart(blade.length~species,data=d,method="jitter",vertical=TRUE,
           pch=1,ylab="Lamina length (cm)",xaxt="n",las=2, cex.axis = 1.2, cex.lab = 1.2)
set.seed(123)
stripchart(blade.length~species,data=d,method="jitter",vertical=TRUE,
           add=TRUE,pch=16,axes=FALSE,col=cols,las=2)
points(1:3+.2, tapply(d$blade.length,d$species, mean),pch=20)

lines(c(1,1) +0.2, c(v[14,4],v[14,5])) 
lines(c(1.17,1.23) ,y= c(v[14,4],v[14,4]))
lines(c(1.17,1.23) ,y= c(v[14,5],v[14,5]))

lines(c(2,2) +0.2, c(v[15,4],v[15,5])) 
lines(c(2.17,2.23) ,y= c(v[15,4],v[15,4]))
lines(c(2.17,2.23) ,y= c(v[15,5],v[15,5]))

lines(c(3,3)+.2, c(v[13,4],v[13,5]))
lines(c(3.17,3.23) ,y= c(v[13,4],v[13,4]))
lines(c(3.17,3.23) ,y= c(v[13,5],v[13,5]))

#mid2 left: blade width
set.seed(123)
stripchart(blade.width~species,data=d,method="jitter",vertical=TRUE,
           pch=1,ylab="Lamina width (cm)",xaxt="n",las=2, cex.axis = 1.2, cex.lab = 1.2)
set.seed(123)
stripchart(blade.width~species,data=d,method="jitter",vertical=TRUE,
           add=TRUE,pch=16,axes=FALSE,col=cols,las=2)
points(1:3+.2, tapply(d$blade.width,d$species, mean),pch=20)

lines(c(1,1)+0.2, c(v[11,4],v[11,5])) 
lines(c(1.17,1.23) ,y= c(v[11,4],v[11,4]))
lines(c(1.17,1.23) ,y= c(v[11,5],v[11,5]))

lines(c(2,2)+0.2, c(v[12,4],v[12,5])) 
lines(c(2.17,2.23) ,y= c(v[12,4],v[12,4]))
lines(c(2.17,2.23) ,y= c(v[12,5],v[12,5]))

lines(c(3,3)+.2, c(v[10,4],v[10,5]))
lines(c(3.17,3.23) ,y= c(v[10,4],v[10,4]))
lines(c(3.17,3.23) ,y= c(v[10,5],v[10,5]))

#mid2 right: sepal length
set.seed(123)
stripchart(sepal.length~species,data=d,method="jitter",vertical=TRUE,
           pch=1,ylab="Sepal length (cm)",xaxt="n",las=2, cex.axis = 1.2, cex.lab = 1.2)
set.seed(123)
stripchart(sepal.length~species,data=d,method="jitter",vertical=TRUE,
           add=TRUE,pch=16,axes=FALSE,col=cols,las=2)
points(1:3+.2, tapply(d$sepal.length,d$species, mean),pch=20)

lines(c(1,1) +0.2, c(v[8,4],v[8,5])) 
lines(c(1.17,1.23) ,y= c(v[8,4],v[8,4]))
lines(c(1.17,1.23) ,y= c(v[8,5],v[8,5]))

lines(c(2,2) +0.2, c(v[9,4],v[9,5])) 
lines(c(2.17,2.23) ,y= c(v[9,4],v[9,4]))
lines(c(2.17,2.23) ,y= c(v[9,5],v[9,5]))

lines(c(3,3) +0.2, c(v[7,4],v[7,5])) 
lines(c(3.17,3.23) ,y= c(v[7,4],v[7,4]))
lines(c(3.17,3.23) ,y= c(v[7,5],v[7,5]))

#bottom left: sepal width
set.seed(123)
stripchart(sepal.width~species,data=d,method="jitter",vertical=TRUE,
           pch=1,ylab="Sepal width (cm)",xaxt="n",las=2, cex.axis = 1.2,cex.lab = 1.2)
set.seed(123)
stripchart(sepal.width~species,data=d,method="jitter",vertical=TRUE,
           add=TRUE,pch=16,axes=FALSE,col=cols,las=2)
points(1:3+.2, tapply(d$sepal.width,d$species, mean),pch=20)

lines(c(1,1) +0.2, c(v[20,4],v[20,5])) 
lines(c(1.17,1.23) ,y= c(v[20,4],v[20,4]))
lines(c(1.17,1.23) ,y= c(v[20,5],v[20,5]))

lines(c(2,2) +0.2, c(v[21,4],v[21,5])) 
lines(c(2.17,2.23) ,y= c(v[21,4],v[21,4]))
lines(c(2.17,2.23) ,y= c(v[21,5],v[21,5]))

lines(c(3,3)+.2, c(v[19,4],v[19,5]))
lines(c(3.17,3.23) ,y= c(v[19,4],v[19,4]))
lines(c(3.17,3.23) ,y= c(v[19,5],v[19,5]))

#bottom right: blank w/ legend
#plot.new()
#savefont <- par(font=3)
#legend("center", inset=c(0,-.5), 
      # legend=c("A. flavescens    ", "A. formosa",""),cex=1.5, 
     #  pt.cex=1.5,xpd=NA, bty="n")
#par(savefont)
#legend("center", inset=c(0,-.5), legend=c("", "","hybrids             "), 
     #  pch=1, cex=1.5,pt.cex=1.5,xpd=NA, bty="n")
#legend("center", inset=c(0,-.5), legend=c("", "","hybrids             "), 
     #  pch=16, cex=1.5,pt.cex=1.5,col=c("gold", "orangered2","hotpink"),
     #  xpd=NA, bty="n")

#dev.off()


###########################################################################
# RDA ---------------------------------------------------------------------
###########################################################################
expl <- model.matrix(~d$species)
resp <- d[,c(2,4:10)]
parents <- droplevels(as.factor(resp[c(1:7,9:78),1]))
colnames(resp) <- c("species","corolla width","spur length","blade length",
                    "blade width","anther exsertion","sepal length","sepal width")

#run for all parental specimens
RDA.all <- rda(resp[c(1:78),-1] ~expl[c(1:78),],scaling=1)
s <- summary(RDA.all) #57.21 % of cumulative variance ***when MK08 is left in
s
weights <- s$species[,1] #weightings of traits onto first axis
ind <- order(abs(weights),decreasing=TRUE)
weights[ind]

s$species[,1] #weightings of traits onto first axis
anova.cca(RDA, permutations = 9999)

#now exclude MK08 (pink sepals)
RDA.train <- rda(resp[c(1:7,9:78),-1] ~expl[c(1:7,9:78),],scaling=1)
s.some <- summary(RDA.train) #56.99% var. explained by RDA1 ***when MK08 is left out
s.some  #22.7% var. explained by PCA1


#predict hybrid scores for each axis
hy.rda <- predict(RDA.train,type="wa",model="CCA",newdata=resp[79:110,-1], scaling=1)
hy.pca <- predict(RDA.train,type="wa",model="CA",newdata=resp[79:110,-1],scaling=1)
mk08.rda <- predict(RDA.train,type="wa",model="CCA",newdata=resp[8,-1],scaling=1)
mk08.pca <- predict(RDA.train,type="wa",model="CA",newdata=resp[8,-1],scaling=1)

#individuals are "sites" traits are "species"
scrs <- scores(RDA.train, display = c("sites", "species"))

# RDA biplot --------------------------------------------------------------
colvec <- c(rep("gold",19), rep("orangered2",58))
pch <- c(rep("24", 19), rep("21",38), rep("23", 20))

#pdf("figures/RDA.collections.pdf")
par(bty="l")
xlim <- with(scrs, range(species[,1], sites[,1])+c(-.4,0))
ylim <- with(scrs, range(species[,2],sites[,2])+c(-.1,.3))
par(bty="l")
plot.new()
plot.window(xlim=xlim, ylim=ylim)
ordiellipse(RDA.train, groups=parents,
            display = "sites", conf = 0.95, label = FALSE,lwd=1.5)
#add MK points
points(scrs$sites[1:19,], pch = 24, col = "black", bg = "#F5FF0A")
#add WG points
points(scrs$sites[20:57,], pch = 23, col = "black", bg = "#F50000")
#add RL points
points(scrs$sites[58:77,], pch = 21, col = "black", bg = "#F50000")

#add MK08 point
#points(y = c(mk08.pca[,1],mk08.pca[,1],mk08.pca[,1]), x=c(mk08.rda, mk08.rda,mk08.rda),
#       pch = c(24,25,11), bg= c("gold1", "gold1", "black"), col = c("gold1","gold1", "black"))
#add PP points
points(hy.pca[,1]~hy.rda, pch = 4, col = "#FF3399", lwd = 2)
#ordipointlabel(RDA,display="species",add=TRUE)
box()
axis(1)
axis(2)
title(xlab="RDA1 (56.99% variance explained",ylab="PCA1 (22.7 % variance explained)" )

#savefont <- par(font=3)
#legend("topright", legend = c("A. flavescens    ","A. formosa",""), bty = "n",col = c("gold","orangered2","hotpink"), pch = 20, cex=.9)
#par(savefont)
#legend("topright",legend =  c("             ","          ","putative hybrids"),bty="n",col=c("gold","orangered2","hotpink"),pch=20,cex=.9)
#dev.off()
system('open "figures/RDA.collections.pdf"')

##################visual verification of model prediction################

expl <- model.matrix(~d$species)
resp <- d[,c(2,4:10)]
parents <- droplevels(as.factor(resp[c(1:7,9:78),1]))
colnames(resp) <- c("species","corolla width","spur length","blade length","blade width",
                    "anther exsertion","sepal length","sepal width")

resp.i <- resp[-8,]
parents.resp <- resp.i[1:77,]
expl.i <- expl[-8,]
parents.expl <- expl.i[1:77,]


######################### PREDICTION ASSESSMENT #############
count <- 0 
count.fla <- 0
count.fo <- 0
################RUN THE NEXT CODE SECTION 1000 times, tally incorrect IDs ################

i <- sample(2:77, 1) #select random individual to remove
parents.i <- droplevels(as.factor(parents.resp[-i,1]))
if(1 <= i & i <=19){col.i = "gold1"} else {col.i = "orangered2"}

RDA.test.i <- rda(parents.resp[-i,-1] ~parents.expl[-i,],scaling=1)

new.x <- predict(RDA.test.i, type = "wa", model = "CCA", 
                   newdata = parents.resp[i,-1], scaling = 1)
new.y <- predict(RDA.test.i, type="wa", model="CA", 
                 newdata=resp[i, -1], scaling=1)

scrs.i <- scores(RDA.test.i, display = c("sites", "species"))

par(bty="l")
xlim <- with(scrs.i, range(species[,1], sites[,1])+c(-.4,0))
ylim <- with(scrs.i, range(species[,2],sites[,2])+c(-.2,.2))
par(bty="l")
plot.new()
plot.window(xlim=xlim, ylim=ylim)
ordiellipse(RDA.test.i, groups=parents.i,
            display = "sites", conf = .99, label = FALSE,lwd=1.5)

points(y = new.y[1], x = new.x, pch = 16, col = col.i)

count <- count + 1
if(1 <= i & i <=19){count.fla <- count.fla + 1 } else {count.fo <- count.fo + 1}

##################################
# totals for points falling outside of .99 confidence ellipse
#ran 1000 times
# fla total = 225.
# fla outside ellipse = 72
# fla fraction = 72/225 = 0.32

# fo total = 775
# fo outside ellipse = 58
# fo fraction = 58/775 = 0.07483871
p.outside <- 72/225 + 58/775
##################################

# calculate exact binomial probability of observing 
# at least 31 individuals outside of parental ellipses

# X~Binom(32, p). Calculate P(X >= 31)

P <- choose(32,30)*(p.outside^30)*(1-p.outside)^2 + choose(32,31)*(p.outside^31)*(1-p.outside) + choose(32,32)*(p.outside^32)
# P << 0.001


######################
#linear discriminant function analysis
library(MASS)

r <- droplevels(resp[c(1:7,9:78),])
z <- lda(species ~ ., data = r)

z$scaling
plot.new()
sc <- predict(z)$x
stripchart(sc[1:19], vertical = FALSE, at = 1, ylim = c(.5,2.5), xlim = range(sc), 
           pch = 24, bg = "gold1", col = "black", method = "jitter")
points(x = sc[20:77], y = jitter(rep(2, 58),amount = .2), 
       pch = c(rep(21, 20), rep(23, 38)), bg = "orangered1", col = "black")

#predict mk08
mk8 <- predict(z, newdata = r[8,2:8])
#add to plot
points(x = rep(mk8$x, 3), y = c(.8,.8,.8), pch = c(24,25,11),
       bg = c("gold1", "gold1", "black"), col = c("gold1","gold1","black"))

#predict hybrids
r[79:110,]
hy.pred <- predict(z, newdata = resp[79:110,2:8])
points(x = hy.pred$x, y = jitter(rep(1.5,32), amount = .2), 
       pch = 4, lwd = 2, col = "hotpink")
axis(2, at = c(.9,2.1), labels = c("A. flavescens", "A. formosa"), lwd.ticks = 0,
     font = 3, padj = 1)
axis(2, at = c(1.5), labels = "hybrid", lwd.ticks = 0, padj = 1)


#add MK points
points(scrs$sites[1:19,], pch = 24, col = "black", bg = "gold1")
#add WG points
points(scrs$sites[20:57,], pch = 23, col = "black", bg = "orangered2")
#add RL points
points(scrs$sites[58:77,], pch = 21, col = "black", bg = "orangered2")

########### V-fold Cross Validation #############
vlda = function(v,formula,data,cl){
  require(MASS)
  grps = cut(1:nrow(data),v,labels=FALSE)[sample(1:nrow(data))]
  pred = lapply(1:v,function(i,formula,data){
    omit = which(grps == i)
    z = lda(formula,data=data[-omit,])
    predict(z,data[omit,])
  },formula,data)
  
  wh = unlist(lapply(pred,function(pp)pp$class))
  table(wh,cl[order(grps)])
}

############################
mis <- vector()
for(i in 1:1000){
xval <- vlda(5,species~.,data=r, cl=r$species)
mis[i] <- sum(xval[row(xval) != col(xval)]) / sum(xval)
}
mean(mis) # 0.001558442


# morphology- altitude cline investigation -----------------------------------------

p18 <- read.csv("/Users/jeff/Projects/Aquilegia/porcupine.valley.phenotypes.csv", stringsAsFactors = FALSE)

pheno <- predict(z.herb, p18[, -c(1:3,12)], prior = c(.5,.5))$x




pch18 <- p18$site
adj <- grep("adj", pch18)
pch18[adj] <- 15
pch18[-adj] <- 16
pch18 <- as.numeric(pch18)
  
plot(pheno ~ p18$elevation, ylim = range(sc), pch = pch18, xlim = c(1400, 2200),
     ylab = "flavescens - formosa discriminant axis", xlab = "elevation (m)")

plot(pheno ~ p18$elevation, ylim = range(-4,4), pch = pch18, xlim = c(1400, 2200),
     ylab = "flavescens - formosa discriminant axis", xlab = "elevation (m)")



#run model for specimens in both valleys
df18 <- cbind.data.frame(elevation = p18$elevation, pheno)
z18 <- lm(pheno ~ elevation, data = df18)

xs <- range(df18$elevation) #endpoints for regression line
ys <- predict(z18, newdata = data.frame(elevation = xs))

newx18 <- seq(min(df18$elevation, na.rm = TRUE), max(df18$elevation, na.rm = TRUE), length.out=100)
py18 <- predict(z18, newdata = data.frame(elevation = newx18),  interval = 'confidence')

#plot regression line and confidence bands
polygon(c(rev(newx18), newx18), c(rev(py18[ ,3]), py18[ ,2]), col = rgb(190/255,190/255,190/255, .6), border = NA)
lines(y = ys, x = xs)


#subset to include only Porcupine valley
pv18 <- subset(p18, site == "porcupine valley transect")
pv.index <- which(p18$site == "porcupine valley transect")
plot(pheno[pv.index] ~ p18$elevation[pv.index], pch = 16, 
     xlim = c(1800, 2150), ylim = c(-3,3), 
     xlab = "Elevation (m)", ylab = "Flavescens - formosa discriminant axis")

#run model for this valley subset
pv18 <- subset(p18, site == "porcupine valley transect")
pvdf18 <- cbind.data.frame(elevation = p18$elevation[pv.index], pheno = pheno[pv.index])
pvz18 <- lm(pheno ~ elevation, data = pvdf18)

#regression line and confidence bands
newxpv <- range(pvdf18$elevation)
newypv <- predict(pvz18, newdata = data.frame(elevation = newxpv))

xconf <-  seq(min(pvdf18$elevation, na.rm = TRUE), max(pvdf18$elevation, na.rm = TRUE), length.out=100)
yconf <- predict(pvz18, newdata = data.frame(elevation = xconf),  interval = 'confidence')
polygon(c(rev(xconf), xconf), c(rev(yconf[ ,3]), yconf[ ,2]), col = rgb(190/255,190/255,190/255, .6), border = NA)
lines(x = newxpv, y = newypv)

anova(pvz18)







# correlation of LDA1 with log.rg -----------------------------------------

#create log(r/g) variable
r <- rowMeans(cbind(d$red.mean1, d$red.mean2), na.rm= TRUE)
r[which(r=="NaN")] <- NA
g <- rowMeans(cbind(d$green.mean1, d$green.mean2), na.rm = TRUE)
g[which(g=="NaN")] <- NA
rg <- r/g
log_rg <- log(r/g)

#visualize with histogram
hist(log_rg[hy])
hist(log_rg[fo])
hist(log_rg[fl], breaks = seq(-0.05, .25, .01))

log_rg
hyb.lda.fl
plot(hyb.lda.fl ~ log_rg[79:110])
cor(y=as.numeric(hyb.lda.fl),x=log_rg[79:110],use="complete.obs")
# -0.2916165

#calculate significance
set.seed(123)
iter <- 9999
x <- log_rg[79:110]
cnt <-0
cor.obs <- cor(x=log_rg[79:110],y=hyb.lda.fl,use="complete.obs")
for(i in 1:iter){
  z <- x
  z <- z[sample(1:length(z))]
  cor.est <- cor(x=z,y=hyb.lda.fl,use="complete.obs")
  
    if(abs(cor.est) >= abs(cor.obs)){cnt <- cnt+1}

}
cnt
p.value <- round(as.numeric((cnt+1)/(iter+1)),digits=4)
p.value #0.1227

############################################################
############################################################
########### Color method validation #########################
#phenotype <- predict(CCA, newdata = resp[which(d$species == "hybrid"),],type = "wa" )
#cor(hybrids, log_rg[hy])

eye <- d[which(!is.na(d[,18])),18] # this variable contains visual classification
img <- log_rg[which(!is.na(d[,18]))] #this variable contains matching log_rg values
plot(img~eye)
cor(eye,img, use = "complete.obs", method = "spearman") #0.87739281

##################################
iter <- 9999
x <- d[hy,4:10]
cnt.1 <-rep(0,7)
cor.obs.1 <- cor(x[,1], x[,])
for(i in 1:iter){
  y <- x
  y[,1] <- y[sample(1:length(y[,1])),1]
  cor.est.1 <- cor(y[,1], y[,])
  for(j in 1:7){
    if(abs(cor.est.1[j]) >= abs(cor.obs.1[j])){cnt.1[j] <- cnt.1[j]+1}
  }
}
p.value.1 <- round(as.numeric((cnt.1)/(iter)), digits=4)
###################################
cnt.2 <- rep(0,7)
cor.obs.2 <- cor(x[,2], x[,])
for(i in 1:iter){
  y <- x
  y[,2] <- y[sample(1:length(y[,2])),2]
  cor.est.2 <- cor(y[,2], y[,])
  for(j in 1:7){
    if(abs(cor.est.2[j])>= abs(cor.obs.2[j])){cnt.2[j] <- (cnt.2[j]+1)}
  }
}
p.value.2 <- round(as.numeric((cnt.2)/(iter)), digits=4)
####################################
cnt.3 <- rep(0,7)
cor.obs.3 <- cor(x[,3], x[,])
for(i in 1:iter){
  y <- x
  y[,3] <- y[sample(1:length(y[,3])),3]
  cor.est.3 <- cor(y[,3], y[,])
  for(j in 1:7){
    if(abs(cor.est.3[j])>= abs(cor.obs.3[j])){cnt.3[j] <- (cnt.3[j]+1)}
  }
}
p.value.3 <- round(as.numeric((cnt.3)/(iter)), digits=4)
#######################################
cnt.4 <- rep(0,7)
cor.obs.4 <- cor(x[,4], x[,])
for(i in 1:iter){
  y <- x
  y[,4] <- y[sample(1:length(y[,4])),4]
  cor.est.4 <- cor(y[,4], y[,])
  for(j in 1:7){
    if(abs(cor.est.4[j])>= abs(cor.obs.4[j])){cnt.4[j] <- (cnt.4[j]+1)}
  }
}
p.value.4 <- round(as.numeric((cnt.4)/(iter)), digits=4)
p.value.4
########################################
cnt.5 <- rep(0,7)
cor.obs.5 <- cor(x[,5], x[,])
for(i in 1:iter){
  y <- x
  y[,5] <- y[sample(1:length(y[,5])),5]
  cor.est.5 <- cor(y[,5], y[,])
  for(j in 1:7){
    if(abs(cor.est.5[j])>= abs(cor.obs.5[j])){cnt.5[j] <- (cnt.5[j]+1)}
  }
}
p.value.5 <- round(as.numeric((cnt.5)/(iter)), digits=4)
p.value.5
########################################
cnt.6 <- rep(0,7)
cor.obs.6 <- cor(x[,6], x[,])
for(i in 1:iter){
  y <- x
  y[,6] <- y[sample(1:length(y[,6])),6]
  cor.est.6 <- cor(y[,6], y[,])
  for(j in 1:7){
    if(abs(cor.est.6[j])>= abs(cor.obs.6[j])){cnt.6[j] <- (cnt.6[j]+1)}
  }
}
p.value.6 <- round(as.numeric((cnt.6)/(iter)), digits=4)
########################################
cnt.7 <- rep(0,7)
cor.obs.7 <- cor(x[,7], x[,])
for(i in 1:iter){
  y <- x
  y[,7] <- y[sample(1:length(y[,7])),7]
  cor.est.7 <- cor(y[,7], y[,])
  for(j in 1:7){
    if(abs(cor.est.7[j])>= abs(cor.obs.7[j])){cnt.7[j] <- (cnt.7[j]+1)}
  }
}
p.value.7 <- round(as.numeric((cnt.7)/(iter), digits=4))
########################################
trait1 <- c(rep("corrola.width",7),rep("spur.length",7),rep("blade.length",7),rep("blade.width",7), rep("anther.exsertion",7),rep("sepal.length",7),rep("sepal.width",7))
trait2 <- rep(c("corrola.width","spur.length", "blade.length", "blade.width", "anther.exsertion","sepal.length", "sepal.width"),7)

correlation <- c(round(cor.obs.1,digits=2), round(cor.obs.2,digits=2), round(cor.obs.3,digits=2), round(cor.obs.4,digits=2),round(cor.obs.5,digits=2), round(cor.obs.6,digits=2),round(cor.obs.7,digits=2))
p.value <- c(p.value.1,p.value.2,p.value.3,p.value.4,p.value.5,p.value.6,p.value.7)
trait.cor <- data.frame(trait1,trait2,correlation,p.value)
head(trait.cor)
trait.cor <- trait.cor[-which(trait.cor[,3] ==1),]
write.csv(trait.cor,"clean_data/trait_correlations.csv",row.names=FALSE)
#########################################################
library(lattice)
histogram(~d$blade.length | d$species,layout=c(1,3), right=FALSE, breaks=seq(0,1,.1))
histogram(~d$spur.length | d$species, layout = c(1,3), right=FALSE)
histogram(~d$corolla.width | d$species, layout = c(1,3), right=FALSE, breaks=seq(1,3,.2))
histogram(~d$anther.exsertion |d$species,layout=c(1,3), right=FALSE, breaks=seq(.4,2,.1))

with(d, hist(spur.length, right=FALSE))
with(d,hist(blade.length, right=FALSE))













#x <- d[hy,4:10]
#iter <- 9999
#for(i in 1:7){
# cnt.i <- rep(0,6)
# cor.obs.i <- cor(x[,i], x[,])
#    for(j in 1:iter){
#      y <- x
#      y[,i] <- x[sample(1:length(x[,i])),i]
#      cor.est.i <- cor(x[,i], x[,])
#      for(k in 1:6){
#        if(abs(cor.est.i[k])>= abs(cor.obs.i[k])){cnt.i[k] <- (cnt.i[k]+1)
#      }
#    }
#  }
#p.value.i <- round(as.numeric((cnt.i+1)/(iter+1)), digits=4)
#p.value.i
#}
















color <- color_data[-which(color_data$standardized =="no"),]
red_mean <- color$red.mean
green_mean <- color$green.mean
rg_ratio <- red_mean / green_mean
log_rg <- log(rg_ratio)
#take average of repeat measurements
names(log_rg) <- color$id
mk <- log_rg[grep("^MK", names(log_rg))]
log_rg <- log_rg[-(grep("^MK", names(log_rg)))]
repeats <- rowMeans(cbind(mk[c(1,3,5,7,9,11,13,15,18)], mk[c(2,4,6,8,10,12,14,16,19)]))
mk <- mk[-c(2,4,6,8,10,12,14,16,19)]
mk[c(1:8,10)] <- repeats
log_rg <- c(log_rg, mk)
clearwater <- log_rg[grep("^WG", names(log_rg))]
kobau <- log_rg[grep("^MK", names(log_rg))]
porcupine <- log_rg[grep("^PP", names(log_rg))]

#correlation between floral phenotype and log_rg
expl <- model.matrix(~traits$species)





#calculate significance of Porcupine Peak floral intermediacy 
library(vegan)
str(traits)
expl <- model.matrix(~traits$species)
resp <- traits[,4:10]
mycca <- cca(resp[1:78,] ~ expl[1:78,], scale= TRUE)
plot(mycca, scaling = 3)


####################









iter <- 20 #permutations
cnt <- 0
for(i in 1:iter){
  x<-traits
  x[,2] <- traits[sample(1:length(traits[,2])),2]
  x <- x[order(x[,2]),]
  fo <- which(x[,2] == "A. formosa")
  fl <- which(x[,2] == "A. flavescens")
  hy <- which(x[,2] == "hybrid")
  e <- model.matrix(~x[,2])
  #fl=00;fo=10,hy=01
  r <- x[4:10]
  rd <- rda(r[-hy,] ~ e[-hy,], scale = TRUE)
  plot(rd)
  scrs <- as.numeric(scores(rd, choices = 1, display = "sites"))
  pred <- predict(rd, newdata = r, type = "wa")
  if( mean(scrs[fl]) < mean(pred) & mean(pred) < mean(scrs[fo]) | mean(scrs[fo]) < mean(pred) & mean(pred) < mean(scrs[fl]) ){
    cnt <- cnt + 1
  }
}
p_value <- round(as.numeric((cnt +1))/(iter=1),digits =4)
p_value
myrda <- rda(resp[traits$species != "hybrid",] ~ expl[traits$species != "hybrid",])
fl <- which(traits$species =="A. flavescens")
fo <- which(traits$species == "A. formosa")
hy <- which(traits$species == "hybrid")

plot(myrda)
scores(myrda)

shit <- rda(resp~ expl)
plot(shit)

