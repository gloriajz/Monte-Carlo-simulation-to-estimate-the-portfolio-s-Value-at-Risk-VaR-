# Implementation of necessary libraries
install.packages("VineCopula") 
library(VineCopula)
install.packages("fGarch")
library(fGarch)
install.packages("KScorrect")
library(KScorrect)
library(stats)
install.packages("ADGofTest")
library(ADGofTest)
install.packages("goftest") 
library(goftest)
library(fBasics)

# Read datasets
CAC_40 <- read.csv("CAC_40.csv")
DAX_30 <- read.csv("DAX.csv")

# Create dataframes
CAC_missing <- data.frame(CAC_40)
DAX_missing <- data.frame(DAX_30)

# Deal with "null" rows
CAC_missing[CAC_missing == "null"] = NA
CAC <- CAC_missing[complete.cases(CAC_missing), ]
DAX_missing[DAX_missing == "null"] = NA
DAX <- DAX_missing[complete.cases(DAX_missing), ]

# Extract the adjusted closing price column
price_CAC <- as.numeric(CAC$Adj.Close)
price_DAX <- as.numeric(DAX$Adj.Close)

# log-returns of the adjusted closing prices
log_returns_CAC <- diff(log(price_CAC), lag=1)
log_returns_DAX <- diff(log(price_DAX), lag=1)

# Test for normality
jarqueberaTest(log_returns_CAC)
jarqueberaTest(log_returns_DAX)

# Building AR models: the Box - Jenkins approach
# Step 1: Identification

# returns cac
par(mfrow=c(2,2))
acf(log_returns_CAC,lag.max="10",col="green", lwd=2)
pacf(log_returns_CAC,lag.max="10",col="green", lwd=2)
acf(log_returns_CAC ^2, col="red", lwd=2)
par(mfrow=c(1,1))

# returns dax
par(mfrow=c(2,2))
acf(log_returns_DAX,lag.max="10",col="green", lwd=2)
pacf(log_returns_DAX,lag.max="10",col="green", lwd=2)
acf(log_returns_DAX^2, col="red", lwd=2)
par(mfrow=c(1,1))

# Step 2: Estimation : conditional distribution for our model

# First, we estimate for cac40
cac1.1=garchFit(formula=~arma(1,0)+garch(1,1),data=log_returns_CAC,trace=
                        F,cond.dist="norm")
cac1.2=garchFit(formula=~arma(1,0)+garch(1,1),data=log_returns_CAC,trace=
                        F,cond.dist="snorm")
cac1.3=garchFit(formula=~arma(1,0)+garch(1,1),data=log_returns_CAC,trace=
                        F,cond.dist="ged")
cac1.4=garchFit(formula=~arma(1,0)+garch(1,1),data=log_returns_CAC,trace=
                        F,cond.dist="sged")
cac1.5=garchFit(formula=~arma(1,0)+garch(1,1),data=log_returns_CAC,trace=
                        F,cond.dist="std")
cac1.6=garchFit(formula=~arma(1,0)+garch(1,1),data=log_returns_CAC,trace=
                        F,cond.dist="sstd")
ic=rbind(cac1.1@fit$ics,cac1.2@fit$ics,cac1.3@fit$ics,cac1.4@fit$ics,cac1.5
         @fit$ics,cac1.6@fit$ics)
rownames(ic)<-c("norm", "snorm","ged", "sged", "std", "sstd" )
print(ic)

# skewed t distribution is the best, as the BIC is the lowest


# Now, we estimate for dax30
dax1.1 = garchFit(formula=~arma(3,0)+garch(1,1), data=log_returns_DAX, trace=
                          F, cond.dist="norm")
dax1.2 = garchFit(formula=~arma(3,0)+garch(1,1), data=log_returns_DAX, trace=
                          F, cond.dist="snorm")
dax1.3 = garchFit(formula=~arma(3,0)+garch(1,1), data=log_returns_DAX, trace=
                          F, cond.dist="ged")
dax1.4 = garchFit(formula=~arma(3,0)+garch(1,1), data=log_returns_DAX, trace=
                          F, cond.dist="sged")
dax1.5 = garchFit(formula=~arma(3,0)+garch(1,1), data=log_returns_DAX, trace=
                          F, cond.dist="std")
dax1.6 = garchFit(formula=~arma(3,0)+garch(1,1), data=log_returns_DAX, trace=
                          F, cond.dist="sstd")

ic = rbind(dax1.1@fit$ics, dax1.2@fit$ics, dax1.3@fit$ics, dax1.4@fit$ics, 
           dax1.5@fit$ics, dax1.6@fit$ics)
rownames(ic) <- c("norm", "snorm", "ged", "sged", "std", "sstd")
print(ic)

# again, skewed t distribution is the best, as the BIC is the lowest

# We now try to find the best AR + GARCH model combination
# we will consider the significant lags found in step 1 and up to GARCH(3,3)

# first for cac40
cac2.1=garchFit(formula=~arma(1,0)+garch(1,1),data=log_returns_CAC,trace=
                        F,cond.dist="sstd")
cac2.2=garchFit(formula=~arma(1,0)+garch(1,2),data=log_returns_CAC,trace=
                        F,cond.dist="sstd")
cac2.3=garchFit(formula=~arma(1,0)+garch(1,3),data=log_returns_CAC,trace=
                        F,cond.dist="sstd")
cac2.4=garchFit(formula=~arma(1,0)+garch(2,1),data=log_returns_CAC,trace=
                        F,cond.dist="sstd")
cac2.5=garchFit(formula=~arma(1,0)+garch(2,2),data=log_returns_CAC,trace=
                        F,cond.dist="sstd")
cac2.6=garchFit(formula=~arma(1,0)+garch(2,3),data=log_returns_CAC,trace=
                        F,cond.dist="sstd")
cac2.7=garchFit(formula=~arma(1,0)+garch(3,1),data=log_returns_CAC,trace=
                        F,cond.dist="sstd")

ic=rbind(cac2.1@fit$ics,cac2.2@fit$ics,cac2.3@fit$ics,cac2.4@fit$ics,
         cac2.5@fit$ics,cac2.6@fit$ics,cac2.7@fit$ics)
rownames(ic)<-c("AR(1)garch(1,1)","AR(1)garch(1,2)","AR(1)garch(1,3)",
                "AR(1)garch(2,1)","AR(1)garch(2,2)","AR(1)garch(2,3)",
                "AR(1)garch(3,1)")

print(ic)

# here, we can clearly see that the optimal combination is AR(1)+GARCH(1,1) --> lowest BIC

# now for dax30
dax2.1 = garchFit(formula=~arma(3,0)+garch(1,1), data=log_returns_DAX, trace=
                          F, cond.dist="sstd")
dax2.2 = garchFit(formula=~arma(3,0)+garch(1,2), data=log_returns_DAX, trace=
                          F, cond.dist="sstd")
dax2.3 = garchFit(formula=~arma(3,0)+garch(1,3), data=log_returns_DAX, trace=
                          F, cond.dist="sstd")
dax2.4 = garchFit(formula=~arma(3,0)+garch(2,1), data=log_returns_DAX, trace=
                          F, cond.dist="sstd")
dax2.5 = garchFit(formula=~arma(3,0)+garch(3,1), data=log_returns_DAX, trace=
                          F, cond.dist="sstd")

ic = rbind(dax2.1@fit$ics, dax2.2@fit$ics, dax2.3@fit$ics,
             dax2.4@fit$ics, dax2.5@fit$ics)

rownames(ic)<-c("AR(3)garch(1,1)","AR(3)garch(1,2)","AR(3)garch(1,3)",
                  "AR(3)garch(2,1)","AR(3)garch(3,1)")

print(ic)

# here again, we can clearly see that the optimal combination is AR(3)+GARCH(1,1) --> lowest BIC

# Step 3: Model checking
# We will use the residual diagnostics to check our model
# We compute the residuals and consider their correlation and autocorrelation functions

# residuals for cac40
res_cac <- residuals(cac2.1, standardize=TRUE)
par(mfrow=c(2,1))
pacf(res_cac, col="green", lag.max="10", lwd=2)
pacf(res_cac^2, col="red",lag.max="10", lwd=2)
par(mfrow=c(1,1))

#residuals for dax30
res_dax <- residuals(dax2.1, standardize=TRUE)
par(mfrow=c(2,1))
pacf(res_dax, col="green", lag.max="10", lwd=2)
pacf(res_dax^2, col="red",lag.max="10", lwd=2)
par(mfrow=c(1,1))

# in both cases, it seems we have managed to explain most of the correlation with our model
# let's also use the Ljung-Box test to test for autocorrelation between residuals

Box.test(res_cac, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res_dax, lag = 10, type = c("Ljung-Box"), fitdf = 3)

# we want to check if the skewed Student-t is a suitable conditional distribution
# we perform a probability integral transform on the residuals and testing for uniformity

M1<-psstd(res_cac, 0,1, nu=coef(cac2.1)
          ["shape"],xi=coef(cac2.1)["skew"])
hist(M1)
M2<-psstd(res_dax, 0,1, nu=coef(dax2.1)
          ["shape"],xi=coef(dax2.1)["skew"])
hist(M2)

# Kolmogorov-Smirnov and Anderson-Darling tests
KStest1<-LcKS(u1, cdf = "punif")
KStest1$p.value

ADtest1<-ad.test(u1, null="punif")
ADtest1$p.value

# We pass the test for uniformity, so we can proceed to copula modelling
# Using BiCopSelect function fit various copulas to the dataset and select
# the copula that provides the best fit based on the AIC criterion

model=BiCopSelect(M1, M2, familyset=NA, selectioncrit="AIC", indeptest=TRUE, level=0.05,se = TRUE)
model

# The model of best fit is BB1 copula with the estimated parameter theta=0.25.

# Estimate the Value-at-Risk using the Monte Carlo simulation 
# approach based on copula theory
# Value-at-Risk using Monte Carlo simulation
N=10000
# Using function BiCopSim simulate from a BB1 copula with estimated 
# parameter theta=0.25.
u_sim=BiCopSim(N, family=model$family, model$par,  model$par2)
# Assuming marginal models are N(0,1)
# Apply the component-wise Inverse Probability Integral Transform (IPIT).
res1_sim=qnorm(u_sim[,1], mean = 0, sd = 1) 
res2_sim=qnorm(u_sim[,2], mean = 0, sd = 1) 

# "res1_sim" and "res2_sim" are i.i.d.
# Re-introduce autocorrelation and GARCH effects observed in data
alpha1 <- coef(model1)["alpha1"]
beta1 <- coef(model1)["beta1"]
omega1 <- coef(model1)["omega"]
epsi1 <- rnorm(10000, 0, 1)
c1 <- coef(model1)["mu"]
sigma1 <-model1@sigma.t
sigma1_1252 <- omega1 + alpha1*(res1[1251])^2+beta1*(sigma1[1251])^2
y1simulated <- c1 + alpha1*res1[1251]+sigma1_1252*res1_sim

alpha2 <- coef(model2)["alpha1"]
beta2 <- coef(model2)["beta1"]
omega2 <- coef(model2)["omega"]
epsi2 <- rnorm(10000, 0, 1)
c2 <- coef(model2)["mu"]
sigma2 <-model2@sigma.t
sigma2_1252 <- omega2 + alpha2*(res2[1251])^2+beta2*(sigma2[1251])^2
y2simulated <- c2 + alpha1*res2[1251]+sigma2_1252*res2_sim

# Compute portfolio log-returns
portsim <- matrix(0, nrow = N, ncol = 1)
varsim <- matrix(0, nrow = 1, ncol = 2)

portsim=log(1+((exp(y1simulated)-1)+(exp(y2simulated)-1))*(1/2))

# The estimated 99% and 95% Value-at-Risk estimates are:
varsim=quantile(portsim,c(0.01,0.05))
print(varsim)

