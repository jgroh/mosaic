#miniana

mini <- read.csv("/Users/jeff/Desktop/miniana.csv", stringsAsFactors = FALSE)

mi <- predict(z, newdata= mini[, c(25, 26, 28, 30, 33)], type = "response")
mi2 <- predict(z2, newdata= mini[, c(25, 26, 28, 30, 33)], type = "response")



mini$type <- factor(mini$type, 
       levels = levels(as.factor(mini$type) )[order( tapply( 
         mi, mini$type, mean) ) ] )


plot.new()
par(bty = "l")
plot.window(xlim = c(.1,.85), ylim =c(0,1))
points(y = mi2[which(mini$type == "paratype")], x = jitter(rep(.25, 4), 3), pch = 16)

points(y = mi2[which(mini$type == "type")], x = .5, pch = 16)

points(y = mi2[which(mini$type == "isotype")], x = jitter(rep(.75, 8), 3), pch = c(17,17, rep(16, 6)))

axis(1, labels = c("paratypes", "holotype", "isotypes"), at = c(.25, .5, .75))
axis(2)
box()
mtext(2, text = "Logistic regression hybrid index", line = 3 )
legend('topleft', pch = 17, legend = "annotated by Whittemore")

