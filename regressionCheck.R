# Regression concept check , calc by hand

setwd("/Users/gaeaprimeturman/Desktop")

reg.data <-read.csv("reg.csv")
print(reg.data)

print(var(reg.data$Predictor))

#reg by hand

# plot data
plot(reg.data)
# calculate estimates of slope and intercept using variance covariance matrix
(var.cov.reg.data <- var(reg.data))
(b <- var.cov.reg.data[2,1] / var.cov.reg.data[1,1])
(a <- mean(reg.data$Response) - b * mean(reg.data$Predictor))
# add best fit line to plot
abline(a,b)
# calculate y-hat for each value of the predictor variable
(y.hat <- a + b * reg.data$Predictor)
# calculate the sum of squares regression and residuals
(ss.reg <- sum((y.hat - mean(reg.data$Response))^2))
(ss.res <- sum((reg.data$Response - y.hat)^2))
# you should be able to figure out the remaining code to calculate the
# p-value by hand.
slope <- ss.res/ss.reg
print(paste("Slope is:", slope))


#var and cov both give the same answer
print(var(reg.data))

print(cov(reg.data))


#Regression all @ once in one line
fit<-lm(Response~Predictor, data=reg.data)
print(summary(fit))

anova(fit)
names(fit) 
fit$residuals #Gives the residuals
fit$coefficients #Gives the residuals

#Both same
cor.test(reg.data$Predictor, reg.data$Response)
cor.test(reg.data$Response, reg.data$Predictor)
#Same

cor.test(reg.data$Predictor, reg.data$Response, alt = "less")
cor.test(reg.data$Predictor, reg.data$Response, method = "spearman")

#Create some more fake data
reg.data2 <- data.frame(reg.data, Predictor2 = rpois(12, 3), Response2 = sort(rnorm(12, 23)), CatPredictor = gl(3, 4, labels = c("Low", "Med", "High")))
#Plot the all the fake data
plot(reg.data2)

plot(Response ~ Predictor, data = reg.data, pch = 16)
abline(fit, col = "grey")

i <- round(fit$coef[[1]], 2) #i for the intercept
s <- round(fit$coef[[2]], 2) #s for the slope
text(2, 13, paste("y=", i, "+", s, "*x"), pos = 4)

plot(Response ~ Predictor, data = reg.data, pch = 16)

corInfo <- cor.test(reg.data$Predictor, reg.data$Response)
text(2, 13, paste("r = ", signif(corInfo[["estimate"]], 3), "\nP = ", signif(corInfo[["p.value"]], 3), sep = ""), pos = 4) 
#\n returns to next line, #duh R, duh.

#########################>>>FRUITY DATA####################################################>>>

cat("\n")

print("fruity data")
print("predicting diameter from weight")
fruity.data<- read.csv("fruit.csv")
print(fruity.data)

fruit_fit<-lm(Response~Predictor,data=fruity.data)
summary(fruit_fit)
print(fruit_fit)

print(anova(fruit_fit))

print(cor.test(fruity.data$Response, fruity.data$Predictor))

plot(Response ~ Predictor, data = fruity.data, pch = 16)

corInfo <- cor.test(fruity.data$Predictor, fruity.data$Response)
text(2, 13, paste("r = ", signif(corInfo[["estimate"]], 3), "\nP = ", signif(corInfo[["p.value"]], 3), sep = ""), pos = 4) #\n returns to next line

