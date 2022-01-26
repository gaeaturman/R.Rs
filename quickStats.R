#### stats R examp.
setwd("/Users/gaeaprimeturman/Desktop")
se<- function(x) {sd(x)/sqrt(length(x))}


# Calculating p-val for an F-ratio of 4.1 w/ 3 and 12 df's between and within treatments..  
pvL = 1-pf(4.1, 3, 12)


flower.data <- read.csv("snowgrowth.csv")


(mean.flower <- tapply(flower.data$y, flower.data$x, mean))

# check means
print(mean.flower)

low <- c(72.2, 104.5, 65, 75.8, 78.9, 77.5, 82.1, 58, 107.6,
         119, 55.4, 109.5, 94.4, 70.4, 60.5, 67.6, 62,
         106.4, 56, 79.1, 42.9, 124, 79.3, 101.5, 114.9,
         27.4, 73.7, 76.5, 50.8, 91.5, 84.4, 58.1, 105.9,
         87.9, 58.1, 56.2, 103.7, 88.5, 83, 106.2, 86, 63.1,
         48.1, 65.3, 85.3, 73.3, 62.4, 87.5, 57.6, 82.5)
med <- c(46.4, 92.4, 68.7, 58.5, 82.5, 63.4, 60.9, 57.4,
         60.8, 81.4, 96.5, 59.3, 82.6, 67.5, 70.6, 47.5,
         63.5, 54.9, 89.3, 75.7, 85.9, 81.8, 74.8, 55.9, 93,
         116.3, 57.5, 69.6, 103.3, 55.5, 86.4, 53.7, 75.3,
         61.4, 73, 90.1, 63.4, 58.3, 28.9, 53.2, 36.3, 57,
         91.4, 84.1, 77.6, (63.5 + value15), (64.5 + value15),
         (65.5 + value15), (66.5 + value15), (67.5 + value15))

high <- c(61, 39.5, 75.4, 49.1, 90.2, 52.2, 30.7, 116.8,
          62.2, 79.7, 47.8, 76, 77.8, 79, 82.6, 70.5, 112.7,
          43.8, 62, 83.4, 39, 64.5, 47, 81.5, 49.7, 92.9,
          80.1, 60.7, 72.7, 43.1, 77.1, 93.7, 55.9, 80.8,
          106.2, 79.3, 83.9, 66.4, 66.7, 41.5, 67.5, 57.8,
          74, 79.7, 89.8, (58.5 + value25), (59.5 + value25),
          (60.5 + value25), (61.5 + value25), (62.5 + value25))

low.var <- var(low)
med.var <- var(med)
high.var <- var(high)
print(paste("variances for High, Med, and Low treatment respectively are: ", low.var,", ", med.var,", ", high.var,"."))
  
low.sd <- sd(low)
med.sd <- sd(med)
high.sd <- sd(high)
print(paste("standard deviation for High, Med, and Low treatment respectively are: ", low.sd,", ", med.sd,", ", high.sd,"."))  
  
low.se <- se(low)
med.se <- se(med)
high.se <- se(high)
print(paste("standard error for High, Med, and Low treatment respectively are: ", low.se,", ", med.se,", ", high.se,"."))


plot(y ~ x, flower.data, xlab = "Snow Level", ylab = "Plant height (mm)")

#Plot 95% confidence intervals
flower_means <- with(flower.data, tapply(flower.data$y, flower.data$x, mean))
flower_se <- with(flower.data, tapply(flower.data$y, flower.data$x, function(x) sd(x)/sqrt(length(x))))


x_vals <- barplot(flower_means,
                  ylim = c(0, 50),
                  ylab = "Plant Height (mm)"
)
arrows(x_vals, flower_means + flower_se, x_vals, flower_means - flower_se, angle = 90, code = 3)
box()


print(shapiro.test(flower.data$y))



(n.flower <- tapply(flower.data$y, flower.data$x, length))
(N.flower <- length(flower.data$y))
(k.flower <- length(unique(flower.data$x)))
(df.groups.flower <- k.flower - 1)
(df.error.flower <- N.flower - k.flower)
sum_squares <- function(x) {
  return(sum((x - mean(x))^2))
}
(grand.mean.flower <- mean(flower.data$y))
(ss.groups.flower <- sum(n.flower*(mean.flower - grand.mean.flower)^2))
(ss.error.flower <- sum(tapply(flower.data$y, flower.data$x, sum_squares)))
(ms.groups.flower <- ss.groups.flower / df.groups.flower)
(ms.error.flower <- ss.error.flower / df.error.flower)
(f.flower <- ms.groups.flower / ms.error.flower)
print(f.flower)
print(1 - pf(f.flower, 2, 9))


fit_lm <- lm(y ~ x, flower.data)
print(anova(fit_lm))

fit_aov <- aov(y ~ x, flower.data)
summary(fit_aov)



flw.lm.data <- subset(flower.data, x == "Low" | x == "Med", select = c(x,y))
flw.mh.data <- subset(flower.data, x == "Med" | x == "High", select = c(x,y))
flw.lh.data <- subset(flower.data, x == "Low" | x == "High", select = c(x,y))

print("Analysis of treatment combinations...")
cat("\n")
fit_lm <- lm(y ~ x, flw.lm.data)
print("For Low and Med")
print(anova(fit_lm))
cat("\n")

fit_lm <- lm(y ~ x, flw.mh.data)
print("For Med and High")
print(anova(fit_lm))
cat("\n")

fit_lm <- lm(y ~ x, flw.lh.data)
print("For Low and High")
print(anova(fit_lm))
