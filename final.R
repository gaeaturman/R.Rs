g.test <- function(x, y = NULL, correct="williams",
                   p = rep(1/length(x), length(x)), simulate.p.value = FALSE, B = 2000)
  #can also use correct="none" or correct="yates"
{
  DNAME <- deparse(substitute(x))
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x)) {
    if (min(dim(x)) == 1) 
      x <- as.vector(x)
  }
  if (!is.matrix(x) && !is.null(y)) {
    if (length(x) != length(y)) 
      stop("x and y must have the same length")
    DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    OK <- complete.cases(x, y)
    x <- as.factor(x[OK])
    y <- as.factor(y[OK])
    if ((nlevels(x) < 2) || (nlevels(y) < 2)) 
      stop("x and y must have at least 2 levels")
    x <- table(x, y)
  }
  if (any(x < 0) || any(is.na(x))) 
    stop("all entries of x must be nonnegative and finite")
  if ((n <- sum(x)) == 0) 
    stop("at least one entry of x must be positive")
  #If x is matrix, do test of independence
  if (is.matrix(x)) {
    #Test of Independence
    nrows<-nrow(x)
    ncols<-ncol(x)
    if (correct=="yates"){ # Do Yates' correction?
      if(dim(x)[1]!=2 || dim(x)[2]!=2) # check for 2x2 matrix
        stop("Yates' correction requires a 2 x 2 matrix")
      if((x[1,1]*x[2,2])-(x[1,2]*x[2,1]) > 0)
      {
        x[1,1] <- x[1,1] - 0.5
        x[2,2] <- x[2,2] - 0.5
        x[1,2] <- x[1,2] + 0.5
        x[2,1] <- x[2,1] + 0.5
      }
      else
      {
        x[1,1] <- x[1,1] + 0.5
        x[2,2] <- x[2,2] + 0.5
        x[1,2] <- x[1,2] - 0.5
        x[2,1] <- x[2,1] - 0.5
      }
    }
    
    sr <- apply(x,1,sum)
    sc <- apply(x,2,sum)
    E <- outer(sr,sc, "*")/n
    # are we doing a monte-carlo?
    # no monte carlo GOF?
    if (simulate.p.value){
      METHOD <- paste("Log likelihood ratio (G-test) test of independence\n\t with simulated p-value based on", B, "replicates")
      tmp <- .C("gtestsim", as.integer(nrows), as.integer(ncols),
                as.integer(sr), as.integer(sc), as.integer(n), as.integer(B),
                as.double(E), integer(nrows * ncols), double(n+1),
                integer(ncols), results=double(B), PACKAGE= "ctest")
      g <- 0
      for (i in 1:nrows){
        for (j in 1:ncols){
          if (x[i,j] != 0) g <- g + x[i,j] * log(x[i,j]/E[i,j])
        }
      }
      STATISTIC <- G <- 2 * g
      PARAMETER <- NA
      PVAL <- sum(tmp$results >= STATISTIC)/B
    }
    else {
      # no monte-carlo
      # calculate G
      g <- 0
      for (i in 1:nrows){
        for (j in 1:ncols){
          if (x[i,j] != 0) g <- g + x[i,j] * log(x[i,j]/E[i,j])
        }
      }
      q <- 1
      if (correct=="williams"){ # Do Williams' correction
        row.tot <- col.tot <- 0    
        for (i in 1:nrows){ row.tot <- row.tot + 1/(sum(x[i,])) }
        for (j in 1:ncols){ col.tot <- col.tot + 1/(sum(x[,j])) }
        q <- 1+ ((n*row.tot-1)*(n*col.tot-1))/(6*n*(ncols-1)*(nrows-1))
      }
      STATISTIC <- G <- 2 * g / q
      PARAMETER <- (nrow(x)-1)*(ncol(x)-1)
      PVAL <- 1-pchisq(STATISTIC,df=PARAMETER)
      if(correct=="none")
        METHOD <- "Log likelihood ratio (G-test) test of independence without correction"
      if(correct=="williams")
        METHOD <- "Log likelihood ratio (G-test) test of independence with Williams' correction"
      if(correct=="yates")
        METHOD <- "Log likelihood ratio (G-test) test of independence with Yates' correction"
    }
  }
  else {
    # x is not a matrix, so we do Goodness of Fit
    METHOD <- "Log likelihood ratio (G-test) goodness of fit test"
    if (length(x) == 1) 
      stop("x must at least have 2 elements")
    if (length(x) != length(p)) 
      stop("x and p must have the same number of elements")
    E <- n * p
    
    if (correct=="yates"){ # Do Yates' correction
      if(length(x)!=2)
        stop("Yates' correction requires 2 data values")
      if ( (x[1]-E[1]) > 0.25) {
        x[1] <- x[1]-0.5
        x[2] <- x[2]+0.5
      }
      else if ( (E[1]-x[1]) > 0.25){
        x[1] <- x[1]+0.5
        x[2] <- x[2]-0.5
      }
    }
    names(E) <- names(x)
    g <- 0
    for (i in 1:length(x)){
      if (x[i] != 0) g <- g + x[i] * log(x[i]/E[i])
    }
    q <- 1
    if (correct=="williams"){ # Do Williams' correction
      q <- 1+(length(x)+1)/(6*n)
    }
    STATISTIC <- G <- 2*g/q
    PARAMETER <- length(x) - 1
    PVAL <- pchisq(STATISTIC, PARAMETER, lower = FALSE)
  }
  names(STATISTIC) <- "Log likelihood ratio statistic (G)"
  names(PARAMETER) <- "X-squared df"
  names(PVAL) <- "p.value"
  structure(list(statistic=STATISTIC,parameter=PARAMETER,p.value=PVAL,
                 method=METHOD,data.name=DNAME, observed=x, expected=E),
            class="htest")
}







setwd("/Users/gaeaprimeturman/Desktop")

wID<- 01269599

v1<- 100-95
v2<- 100-99

##Q1
# ANOVA or t-test, for both 

cell.data <- read.csv("cellsStuff.csv")

fit_lm <- lm(y ~ x, cell.data)
F<- (anova(fit_lm)[["F value"]])
print(paste("1a. test statistic value is:", F)) ## 1(a)

print(paste("1b. p-value is:", (1 - pf(F,2,147)))) ## 1(b)

x<- c(0.197, 0.363, 0.545, 0.121, 0.182, 0.08, 0.47, 0.131, 0.094, (v1/100))
y<- c(142.9, 13.4, 54.3, 118.9, 50.4, 109, 7, 91.3, 122.7, v2)

cor.data <- data.frame(x,y)
x.m<- mean(x)
y.m<- mean(y)

sp.xy<- sum((x - x.m) * (y - y.m))
s.xy<- sp.xy / (length(cor.data$x) - 1)

(r<- s.xy / (sd(x) * sd(y)))
print(paste("2a. correlation coefficient is:", r)) ## 2(a)

lm.cor.data <- lm(y ~ x, data = cor.data)
print("2b. get the p-val below:")
print(summary(lm.cor.data)) ## 2(b) {get p-val; ~ 0.1937}

(antibody.mean <- mean(cor.data$x))
(fluor.mean <- mean(cor.data$y))

(var.cov.cor.data <- cov(cor.data))

print("The slope of the linear regression line is:")
print(b.cor.data <- var.cov.cor.data[1,2] / var.cov.cor.data[1,1]) ##slope
a.cor.data <- fluor.mean - b.cor.data*antibody.mean

cat("\n")
print("The standard error of the slope is:")
print("49.9762") ## see residual standard error from summary

print("t value is: -1.419")
cat<- sqrt((1-(r^2))/(10-2))

print("t val is:")
print(r/cat)

value1<- v1
value2<- v2

##### Q3
# t test 2 sample

pis <- c(0.107, 0.11, 0.173, 0.198, 0.11, 0.138, 0.23,
         0.088, 0.158, 0.087, 0.167, 0.162, 0.22, 0.129,
         (0.13 + value1/2000)) 
der <- c(0.102, 0.113, 0.174, 0.127, 0.192, 0.083, 0.176,
         0.16, 0.109, 0.116, 0.111, 0.089, 0.073, 0.132,
         (0.13 - value2/2000))
print("Q3 value of t and p correct, pis der")
print(t.test(pis,der,var.equal = T, alternative = "greater"))

#### Q4
print("Prob of 119 successes")
print(dbinom(119, size=120, prob=0.25))

val44<- ceiling((value1/5))
print("Prob of at least one success")

printVal<- (dbinom(0, size=25, prob=0.25) + dbinom(1, size=25, prob=0.25))
print(printVal)


#### Q5
mean5<- (value2/3)
variance<- (value2/3)

exactly<- (ceiling(value1/2))
print("Prob of exactly 3 successes:")
print(dpois(3, lambda=mean5))

print("prob of 1 or less")
print(ppois((ceiling(value2/4)), lambda=mean5))


#### Q6
Ltum<- c(419, 401, 417, 414)
Stum<- c(90, (50+value1), (50+value2), (value1+value2))
GrandList<- Ltum+Stum

ok<- Reduce("+", GrandList) 
GrandSum<- (ok)
print(GrandSum)

#exp probs
p.larger<- (Reduce("+", Ltum))/GrandSum
p.smaller<- (Reduce("+", Stum))/GrandSum

p.no<- (419+90)/GrandSum
p.doc<- (401+(50+value1))/GrandSum
p.dex<- (417+(50+value2))/GrandSum
p.cyc<- (414+(value1+value2))/GrandSum

#exp combo freqs
f.Lno<- (p.larger*p.no*GrandSum)
f.Ldoc<- (p.larger*p.doc*GrandSum)
f.Ldex<- (p.larger*p.dex*GrandSum)
f.Lcyc<- (p.larger*p.cyc*GrandSum)

f.Sno<- (p.smaller*p.no*GrandSum)
f.Sdoc<- (p.smaller*p.doc*GrandSum)
f.Sdex<- (p.smaller*p.dex*GrandSum)

print("6a exp for smaller tumors cyc treatment")
print(f.Scyc<- (p.smaller*p.cyc*GrandSum))

no <- c(smaller = 90, larger = 419)
doc <- c(smaller = 55, larger = 401)
dox <- c(samller = 51, larger = 417)
cyc <- c(smaller = 6, larger = 414)
  
(obs.freq <- matrix(c(no,doc, dox, cyc),2,
              dimnames=list(activity=c("smaller","larger"),
              treatment=c("no treatment","doc treatment", "dox treatment", "cyc treatment"))))

print("P val for G test")
print(g.test(obs.freq,correct="none"))

### Q7

#(a)
groups<- (ceiling(1 + (value1/8)))
chi<- (2*sqrt(value2))
print("Pval is (a)")
print(pChi <- (1 - pchisq(chi, 1)))

chi2<-(value1/8)
print("Pval is (b)")
print(pChi <- (1 - pchisq(chi2, 2)))

tVal<- ((-value2)/30)
n3<- (value1+10)
p3<- (1-pt(tVal,n3-1))
print("Pval is (c)")
print(p3)


### Q8 paired t-test

add <- c(7.4, 7.1, 17.3, 21.7, 14.1, 18.7, 10, 7.1, 14.1,
         20.2, 19.2, 22.6, 10.1, 17.4, 9.8, 10.1, 19.4, 17.4,
         12.6, 15.4, 7.1, 20.1, 10.9, 16.5, 14.9, 18.3, 8.8,
         16.8, 16.9, 15.1, 18.2, 18.9, 11.3, (5 + value1/10),(5 + value2/10))
amb <- c(10.1, 5.4, 14.4, 18, 15.9, 17.9, 12.8, 6.4, 16.1,
         22, 21.8, 24, 10.9, 19.1, 8.8, 12.2, 22.7, 18.7,
         13.8, 19, 7.7, 18.4, 11.2, 14.2, 15.6, 20.4, 14.1,
         17.4, 16.2, 19.9, 15.4, 23.3, 11.6, (5 + value2/10),(5 + value1/10))
diff <- add-amb

print("8a, df in table,  = 34")
print(t.test(diff,mu=0,var.equal = T))
print("8b p < 0.05, 1/significant")
