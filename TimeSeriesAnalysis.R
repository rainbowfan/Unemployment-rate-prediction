##Time Series Analysis
##Unemployment Data in Yolo County California

UR = read.table('CAYOLO3URN.txt', header = TRUE)
#The data has no zero value or missing value.

####initial data transformaiton####
#use the data from 1990 to 2013 to estimate
UR.data = UR[1:288,]
UR.predict = UR[289:300,]

date = UR.data[,1]
rate = UR.data[,2]

#Plot the data 
plot(rate, type = 'l', xaxt = 'n', main = 'Time Series on "Unemployment rate"', xlab = 'Dates')
axis(1, at = seq(1, length(rate), by = 24), labels = date[seq(1, length(rate), by = 24)])

#data transformation#
#Use Boxcox
library(FitAR)
BoxCox.ts(rate)
#lambda = 0.076

#transform
lambda = .076
rate = lambda^(-1)*(rate^lambda - 1)

#Compare by plots
par(mfrow = c(1,2))
plot(UR.data[,2], type = 'l', xaxt = 'n', main = 'Data before transformation', xlab = 'Dates')
axis(1, at = seq(1, length(rate), by = 24), labels = date[seq(1, length(rate), by = 24)])

plot(rate, type = 'l', xaxt = 'n', main = 'Data after transformation', xlab = 'Dates')
axis(1, at = seq(1, length(rate), by = 24), labels = date[seq(1, length(rate), by = 24)])
par(mfrow = c(1,1))

####Analyze the "smooth"components####
#First, deseasonlize the data
#Use Moving average estimation
#d = 12, we have q = 6, N = 25

#Step 1: Filtering
data.matrix = matrix(rate, ncol = 12, byrow = T)
two.sided.filter = filter(rate, sides = 2, c(0.5, rep(1, 11), 0.5)/12)
filter.matrix = matrix(two.sided.filter, ncol = 12, byrow = T) 
mu.matrix = data.matrix - filter.matrix  #This is the filterded matrix

#Step 2: Seasonal estimation
mu.k = colMeans(mu.matrix, na.rm = T)
sk = mu.k - mean(mu.k)
sk.matrix2 = matrix(rep(sk, 24), ncol = 12, byrow = T) #seasonal componenets

deseasonalized2 = as.numeric(t(data.matrix - sk.matrix2))

par(mfrow = c(1,2))
plot(as.numeric(t(sk.matrix2)), type = 'l', main = 'Seasonal components', ylab = '', 
     xlab = 'Months')
plot(as.numeric(t(deseasonalized2)), type = 'l', main = 'Deseasonalized data', ylab = '',
     xlab = 'Months')
par(mfrow = c(1,1))

#Then, detrend the deseasonalized data
#Method 1: polynomial 
n = length(rate)
design.matrix2 = data.frame(y = deseasonalized2, Intercept = rep(1, n)) 
var.name2 = names(design.matrix2)
fit2 = vector('list', 20)
Rsquared2 = numeric(20)
Adj.Rsquared2 = numeric(20)
for(p in 1:20){
  design.matrix2 = cbind(design.matrix2, (1:n)^p)   
  var.name2 = c(var.name2, paste('t', p, sep = ''))  
  colnames(design.matrix2) = var.name2
  fit2[[p]] = lm(y ~ . - 1, data = design.matrix2)  
  Rsquared2[p] = summary(fit2[[p]])$r.squared
  Adj.Rsquared2[p] = summary(fit2[[p]])$adj.r.squared
}

print(Adj.Rsquared2)
plot(Adj.Rsquared2, main = 'Adj.Rsquared2 vs. model size', xlab = 'p')
#4, 6, 9, 11th

#Can also use AIC/BIC
model.selection2 = sapply(fit2, function(i){
  return(c(AIC(i), BIC(i)))
})
rownames(model.selection2) = c('AIC', 'BIC')
print(model.selection2)

#plot out the AIC/BIC vs different model sizes and see how the AIC/BIC behaves. 
par(mfrow = c(1,2))
plot(model.selection2[1,], main = 'AIC vs. model size', xlab = 'p')
plot(model.selection2[2,], main = 'BIC vs. model size', xlab = 'p')
par(mfrow = c(1,1))
#2, 4, 6, 9, 11th

#Look at the 11th model
summary(fit2[[2]])
summary(fit2[[4]])
summary(fit2[[6]])
summary(fit2[[9]])
summary(fit2[[11]])

#Compare by plots to see which one looks reasonable towards the end of the observation period
plot(rate, type = 'l', xaxt = 'n', main = 'Trend estimates for different orders', xlab = 'Dates')
axis(1, at = seq(1, length(rate), by = 24), labels = date[seq(1, length(rate), by = 24)])
lines(fit2[[4]]$fitted, col = 2)
lines(fit2[[6]]$fitted, col = 3)
lines(fit2[[9]]$fitted, col = 4)
lines(fit2[[11]]$fitted, col = 6)
legends = paste('p = ', c(4,6,9,11), sep = '')
legend('bottomright', legend = legends, col = c(2,3,4,6), lty = 1, cex = 0.9)
#4th best

plot(deseasonalized2 - fit2[[4]]$fitted, type = 'l', 
     main = 'Residuals after 4th polynomial detrending', ylab = '', xlab = 'Months')

#Seems like there is still some trend left
#need to detrend again
#Use moving average
detrended = deseasonalized2 - fit2[[4]]$fitted

one.sided.filters = vector('list', 9)
a = seq(0.1, 0.9, by = 0.1)
for(q in 1:9){
  one.sided.filters[[q]][1] = detrended[1]
  for(j in 2:length(rate)){
    one.sided.filters[[q]][j] = one.sided.filters[[q]][j-1]*(1-a[q]) + detrended[j]*a[q]
  }
}

#Plot residuals
par(mfrow = c(3,3))
plot(detrended - one.sided.filters[[1]], main = 'One sided filter: a = 0.1', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(detrended - one.sided.filters[[2]], main = 'One sided filter: a = 0.2', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(detrended - one.sided.filters[[3]], main = 'One sided filter: a = 0.3', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(detrended - one.sided.filters[[4]], main = 'One sided filter: a = 0.4', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(detrended - one.sided.filters[[5]], main = 'one sided filter: a = 0.5', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(detrended - one.sided.filters[[6]], main = 'One sided filter: a = 0.6', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(detrended - one.sided.filters[[7]], main = 'One sided filter: a = 0.7', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(detrended - one.sided.filters[[8]], main = 'one sided filter: q = 0.8', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(detrended - one.sided.filters[[9]], main = 'One sided filter: a = 0.9', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])
par(mfrow = c(1,1))
##Looks like the residual plot does not have an obvious trend from a = 0.5

#Plot the filters
par(mfrow = c(1,2))
plot(detrended, type = 'l', xaxt = 'n', main = 'One sided moving average filter', xlab = 'Dates')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])
lines(one.sided.filters[[5]], col = 2)
plot(detrended - one.sided.filters[[5]], main = 'Residuals after redetrending', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])
par(mfrow = c(1,1))

trend =  one.sided.filters[[5]] +  fit2[[4]]$fitted

res1 = detrended - one.sided.filters[[5]]

##method2: moving average
one.sided.filters2 = vector('list', 9)
a = seq(0.1, 0.9, by = 0.1)
for(q in 1:9){
  one.sided.filters2[[q]][1] = deseasonalized2[1]
  for(j in 2:length(rate)){
    one.sided.filters2[[q]][j] = one.sided.filters2[[q]][j-1]*(1-a[q]) + deseasonalized2[j]*a[q]
  }
}

#plot residulas
par(mfrow = c(3,3))
plot(deseasonalized2 - one.sided.filters2[[1]], main = 'One sided filter: a = 0.1', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(deseasonalized2 - one.sided.filters2[[2]], main = 'One sided filter: a = 0.2', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(deseasonalized2 - one.sided.filters2[[3]], main = 'One sided filter: a = 0.3', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(deseasonalized2 - one.sided.filters2[[4]], main = 'One sided filter: a = 0.4', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(deseasonalized2 - one.sided.filters2[[5]], main = 'one sided filter: a = 0.5', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(deseasonalized2 - one.sided.filters2[[6]], main = 'One sided filter: a = 0.6', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(deseasonalized2 - one.sided.filters2[[7]], main = 'One sided filter: a = 0.7', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(deseasonalized2 - one.sided.filters2[[8]], main = 'one sided filter: q = 0.8', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(deseasonalized2 - one.sided.filters2[[9]], main = 'One sided filter: a = 0.9', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])
par(mfrow = c(1,1))

par(mfrow = c(1,2))
plot(deseasonalized2, type = 'l', xaxt = 'n', main = 'One sided moving average filter', xlab = 'Dates')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])
lines(one.sided.filters2[[5]], col = 2)
plot(deseasonalized2 - one.sided.filters2[[5]], main = 'Residuals after redetrending', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])
par(mfrow = c(1,1))

res2 = deseasonalized2 - one.sided.filters2[[5]]

###Analyzing the "rough" component
#compare res1 vs. res2
par(mfrow = c(1,2))
plot(res1, main = 'Residuals for method 1', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(res2, main = 'Residuals for method 2', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])
par(mfrow = c(1,1))

####analyze the residuals#####
###1. The sample ACF
acf(res1)
acf(res2)

res1.acf = acf(res1)$acf[2:20]
res2.acf = acf(res2)$acf[2:20]
bounds = c(-1.96/sqrt(n), 1.96/sqrt(n))
sum(res1.acf < bounds[2] & res1.acf > bounds[1]) #9/19 are within the bounds
sum(res2.acf < bounds[2] & res2.acf > bounds[1]) #10/19 are within the bounds
#The white noise assumption might be rejected


###2: The Portmanteau test
teststat.1 = numeric(20)
pvals.1 = numeric(20)
teststat.2 = numeric(20)
pvals.2 = numeric(20)
for(i in 1:20){
  test1 = Box.test(res1, lag = i, type = 'Ljung')
  teststat.1[i] = test1$statistic  
  pvals.1[i] = test1$p.value       
  
  test2 = Box.test(res2, lag = i, type = 'Ljung')
  teststat.2[i] = test2$statistic
  pvals.2[i] = test2$p.value
}
#Comparing p-values
pvals.1 < 0.05
pvals.2 < 0.05
#The hypothesis of i.i.d. residuals is rejected at level 0.05. 


###3: rank test
mu.pi = 1/4*n*(n-1)
sigma.pi = sqrt(1/72*n*(n-1)*(2*n+5))
alpha = 0.05
z = qnorm(1 - alpha/2)
#Find Pi for residual 1
Pi.1 = 0
for(j in 1:(n-1)){
  for(i in (j+1):n){
    if(res1[i] > res1[j]) Pi.1 = Pi.1 + 1
  }
}
P1 = abs(Pi.1 - mu.pi)/sigma.pi
P1 > z  #Comparing test statistics with critical value
#Find Pi for residual 2
Pi.2 = 0
for(j in 1:(n-1)){
  for(i in (j+1):n){
    if(res2[i] > res2[j]) Pi.2 = Pi.2 + 1
  }
}
P2 = abs(Pi.2 - mu.pi)/sigma.pi
P2 > z 

#Both results in a false at alpha = 0.05, 
#Thus the assumption that the residuals are i.i.d. cannot be rejected. 
#This says there is no more trend in the data

###### Analyze the "rough" component ######
###res1
lag.plot(res1, lags = 12, layout = c(3,4))

ar1 = vector('list', 20)
aic1 = numeric(20)
for (j in 1:20){
  ar1[[j]] = arima(res1, order = c(j,0,0), method = 'ML')
  aic1[j] = ar1[[j]]$aic
}
aic1

ma1 = vector('list', 20)
aic2 = numeric(20)
for (j in 1:20){
  ma1[[j]] = arima(res1, order = c(0,0,j), method = 'ML')
  aic2[j] = ma1[[j]]$aic
}
aic2

arma1 = matrix(NA, ncol = 10, nrow = 10)

for(i in 1:10){
  for(j in 1:10){
    arma1[i,j] = arima(res1, order = c(i,0,j), method = "ML")$aic
  }
}
arma1

aic3 = numeric(100)
pq = numeric(100)
for(i in 1:10){
  for (j in 1:10){
    aic3[(i-1)*10 + j] = arma1[i,j]
    pq[(i-1)*10 + j] = i + j
  }
}

limit = c(min(aic1, aic2, aic3), max(aic1, aic2, aic3))
par(mfrow = c(1,3))
plot(aic1, ylim = limit, xlab = 'p', main = 'AIC of the fitted AR(p) mode')
plot(aic2, ylim = limit, xlab = 'q', main = 'AIC of the fitted MA(q) model')
plot(aic3 ~ pq, xlab = 'p+q', ylim = limit, main = 'AIC of the fitted ARMA(p,q) model')
points(x = 9, y = -1238.515, col = 'red', pch = 19)
par(mfrow = c(1,1))
# the red point represent arma(5,4) <- best!

###res2
lag.plot(res2, lags = 12, layout = c(3,4))

ar2 = vector('list', 20)
aic1.2 = numeric(20)
for (j in 1:20){
  ar2[[j]] = arima(res2, order = c(j,0,0), method = 'ML')
  aic1.2[j] = ar2[[j]]$aic
}
aic1.2

ma2 = vector('list', 20)
aic2.2 = numeric(20)
for (j in 1:20){
  ma2[[j]] = arima(res2, order = c(0,0,j), method = 'ML')
  aic2.2[j] = ma2[[j]]$aic
}
aic2.2

arma2 = matrix(NA, ncol = 10, nrow = 10)

for(i in 1:10){
  for(j in 1:10){
    arma2[i,j] = arima(res2, order = c(i,0,j), method = "ML")$aic
  }
}
arma2

aic3.2 = numeric(100)
pq.2 = numeric(100)
for(i in 1:10){
  for (j in 1:10){
    aic3.2[(i-1)*10 + j] = arma2[i,j]
    pq.2[(i-1)*10 + j] = i + j
  }
}

limit = c(min(aic1.2, aic2.2, aic3.2), max(aic1.2, aic2.2, aic3.2))
par(mfrow = c(1,3))
plot(aic1.2, ylim = limit, xlab = 'p', main = 'AIC of the fitted AR(p) mode')
plot(aic2.2, ylim = limit, xlab = 'q', main = 'AIC of the fitted MA(q) mode')
plot(aic3.2 ~ pq, xlab = 'p+q', ylim = limit, main = 'AIC of the fitted ARMA(p,q) mode')
points(x = 9, y = -1227.302, col = 'red', pch = 19)
par(mfrow = c(1,1))
# the red point represent arma(5,4) <- best!

#### Check if residuals conform to white noise
#Method 1: The sample ACF
resid1 = res1 - fitted(arima(res1, order = c(5, 0, 4), method = "ML"))
resid2 = res2 - fitted(arima(res2, order = c(5, 0, 4), method = 'ML'))
resid1.acf = acf(resid1)
resid2.acf = acf(resid2)

#Check how many within bound
resid1.acf = resid1.acf$acf[2:20]  #The first acf value is of lag 0, which we ignore here
resid2.acf = resid2.acf$acf[2:20]  #The first acf value is of lag 0, which we ignore here
bounds = c(-1.96/sqrt(n), 1.96/sqrt(n))
sum(resid1.acf < bounds[2] & resid1.acf > bounds[1]) #18/19 are within the bounds
sum(resid2.acf < bounds[2] & resid2.acf > bounds[1]) #18/19 are within the bounds

#Method 2: The Portmanteau test
Q.1 = cumsum(resid1.acf^2) * n
Q.2 = cumsum(resid2.acf^2) * n

Q.1 > qchisq(0.95, df = 1:19)
Q.2 > qchisq(0.95, df = 1:19)

#The hypothesis of i.i.d. residuals is accepted at level 0.05. 

#Method 3: Rank test
mu.pi = 1/4*n*(n-1)
sigma.pi = sqrt(1/72*n*(n-1)*(2*n+5))
alpha = 0.05
z = qnorm(1 - alpha/2)
#Find Pi for residual 1
Pi.1 = 0
for(j in 1:(n-1)){
  for(i in (j+1):n){
    if(resid1[i] > resid1[j]) Pi.1 = Pi.1 + 1
  }
}
P1 = abs(Pi.1 - mu.pi)/sigma.pi
P1 > z  #Comparing test statistics with critical value
#Find Pi for residual 2
Pi.2 = 0
for(j in 1:(n-1)){
  for(i in (j+1):n){
    if(resid2[i] > resid2[j]) Pi.2 = Pi.2 + 1
  }
}
P2 = abs(Pi.2 - mu.pi)/sigma.pi
P2 > z  #Comparing test statistics with critical value

#Both results in a false at alpha = 0.05,
#thus the assumption that the residuals are i.i.d. cannot be rejected. 
#This says there is no more trend in the data

par(mfrow = c(1,2))
plot(resid1, main = 'Residuals for method 1', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])

plot(resid2, main = 'Residuals for method 2', ylab = 'Residuals', xlab = 'Date', xaxt = 'n', type = 'l')
axis(1, at = seq(1, length(rate), by = 48), labels = date[seq(1, length(rate), by = 48)])
par(mfrow = c(1,1))


###### Predict future values #####
#deseasonlize: moving average
#detrend: one-sided moving average
#residuals: arma(5,4) 

seasonality = sk.matrix2[1,]

plot(one.sided.filters2[[5]][265:288], xaxt = 'n', 
     main = 'One sided moving average filter (last 24 points)', xlab = 'Dates')
axis(1, at = seq(1, 24, by = 12), labels = date[seq(265, 288, by = 12)])

#We decide to use the last 12 point to predict m289
#m13
Y = one.sided.filters2[[5]][277:288]
matrix = data.frame(y = Y, Intercept = rep(1, 12)) 
name = names(matrix)
model = vector('list', 5)
R = numeric(5)
AR = numeric(5)
for(p in 1:5){
  matrix = cbind(matrix, (1:12)^p)   #Each time attach a new column
  name = c(name, paste('t', p, sep = ''))  #Append a new column name
  colnames(matrix) = name
  model[[p]] = lm(y ~ . - 1, data = matrix)  
  R[p] = summary(model[[p]])$r.squared
  AR[p] = summary(model[[p]])$adj.r.squared
}
AR
plot(AR)
model13 = model[[3]]
m13 = predict(model13, newdata = data.frame(Intercept = 1, t1 = 13, t2 = 13^2, t3 = 13^3))
m13

#m14
Y = c(one.sided.filters2[[5]][278:288], m13)
matrix = data.frame(y = Y, Intercept = rep(1, 12)) 
name = names(matrix)
model = vector('list', 5)
R = numeric(5)
AR = numeric(5)
for(p in 1:5){
  matrix = cbind(matrix, (1:12)^p)   #Each time attach a new column
  name = c(name, paste('t', p, sep = ''))  #Append a new column name
  colnames(matrix) = name
  model[[p]] = lm(y ~ . - 1, data = matrix)  
  R[p] = summary(model[[p]])$r.squared
  AR[p] = summary(model[[p]])$adj.r.squared
}
AR
plot(AR)

model14 = model[[3]]
m14 = predict(model14, newdata = data.frame(Intercept = 1, t1 = 13, t2 = 13^2, t3 = 13^3))

#m15
Y = c(one.sided.filters2[[5]][279:288], m13, m14)
matrix = data.frame(y = Y, Intercept = rep(1, 12)) 
name = names(matrix)
model = vector('list', 5)
R = numeric(5)
AR = numeric(5)
for(p in 1:5){
  matrix = cbind(matrix, (1:12)^p)   #Each time attach a new column
  name = c(name, paste('t', p, sep = ''))  #Append a new column name
  colnames(matrix) = name
  model[[p]] = lm(y ~ . - 1, data = matrix)  
  R[p] = summary(model[[p]])$r.squared
  AR[p] = summary(model[[p]])$adj.r.squared
}
AR
plot(AR)

model15.2 = model[[2]]
model15.3 = model[[3]]
m15.2 = predict(model15.2, newdata = data.frame(Intercept = 1, t1 = 13, t2 = 13^2))
m15.2
m15.3 = predict(model15.3, newdata = data.frame(Intercept = 1, t1 = 13, t2 = 13^2, t3 = 13^3))
m15.3
m15 = m15.2
#m15.2

#m16
Y = c(one.sided.filters2[[5]][280:288], m13, m14, m15)

matrix = data.frame(y = Y, Intercept = rep(1, 12)) 
name = names(matrix)
model = vector('list', 5)
R = numeric(5)
AR = numeric(5)
for(p in 1:5){
  matrix = cbind(matrix, (1:12)^p)   #Each time attach a new column
  name = c(name, paste('t', p, sep = ''))  #Append a new column name
  colnames(matrix) = name
  model[[p]] = lm(y ~ . - 1, data = matrix)  
  R[p] = summary(model[[p]])$r.squared
  AR[p] = summary(model[[p]])$adj.r.squared
}
AR
plot(AR)

model16 = model[[2]]
m16 = predict(model16, newdata = data.frame(Intercept = 1, t1 = 13, t2 = 13^2))
m16

##m17
Y = c(one.sided.filters2[[5]][281:288], m13, m14, m15, m16)
matrix = data.frame(y = Y, Intercept = rep(1, 12)) 
name = names(matrix)
model = vector('list', 5)
R = numeric(5)
AR = numeric(5)
for(p in 1:5){
  matrix = cbind(matrix, (1:12)^p)   #Each time attach a new column
  name = c(name, paste('t', p, sep = ''))  #Append a new column name
  colnames(matrix) = name
  model[[p]] = lm(y ~ . - 1, data = matrix)  
  R[p] = summary(model[[p]])$r.squared
  AR[p] = summary(model[[p]])$adj.r.squared
}
AR
plot(AR)

model17 = model[[2]]
m17 = predict(model17, newdata = data.frame(Intercept = 1, t1 = 13, t2 = 13^2))
m17

#m18
Y = c(one.sided.filters2[[5]][282:288], m13, m14, m15, m16, m17)

matrix = data.frame(y = Y, Intercept = rep(1, 12)) 
name = names(matrix)
model = vector('list', 5)
R = numeric(5)
AR = numeric(5)
for(p in 1:5){
  matrix = cbind(matrix, (1:12)^p)   #Each time attach a new column
  name = c(name, paste('t', p, sep = ''))  #Append a new column name
  colnames(matrix) = name
  model[[p]] = lm(y ~ . - 1, data = matrix)  
  R[p] = summary(model[[p]])$r.squared
  AR[p] = summary(model[[p]])$adj.r.squared
}
AR
plot(AR)

model18 = model[[2]]
m18 = predict(model18, newdata = data.frame(Intercept = 1, t1 = 13, t2 = 13^2, t3 = 13^3))
m18

#m19
Y = c(one.sided.filters2[[5]][283:288], m13, m14, m15, m16, m17, m18)

matrix = data.frame(y = Y, Intercept = rep(1, 12)) 
name = names(matrix)
model = vector('list', 5)
R = numeric(5)
AR = numeric(5)
for(p in 1:5){
  matrix = cbind(matrix, (1:12)^p)   #Each time attach a new column
  name = c(name, paste('t', p, sep = ''))  #Append a new column name
  colnames(matrix) = name
  model[[p]] = lm(y ~ . - 1, data = matrix)  
  R[p] = summary(model[[p]])$r.squared
  AR[p] = summary(model[[p]])$adj.r.squared
}
AR
plot(AR)

model19 = model[[2]]
m19 = predict(model19, newdata = data.frame(Intercept = 1, t1 = 13, t2 = 13^2, t3 = 13^3))
m19

#m20
Y = c(one.sided.filters2[[5]][284:288], m13, m14, m15, m16, m17, m18, m19)
matrix = data.frame(y = Y, Intercept = rep(1, 12)) 
name = names(matrix)
model = vector('list', 5)
R = numeric(5)
AR = numeric(5)
for(p in 1:5){
  matrix = cbind(matrix, (1:12)^p)   #Each time attach a new column
  name = c(name, paste('t', p, sep = ''))  #Append a new column name
  colnames(matrix) = name
  model[[p]] = lm(y ~ . - 1, data = matrix)  
  R[p] = summary(model[[p]])$r.squared
  AR[p] = summary(model[[p]])$adj.r.squared
}
AR
plot(AR)

model20 = model[[2]]
m20 = predict(model20, newdata = data.frame(Intercept = 1, t1 = 13, t2 = 13^2, t3 = 13^3))
m20

#m21
Y = c(one.sided.filters2[[5]][285:288], m13, m14, m15, m16, m17, m18, m19, m20)
matrix = data.frame(y = Y, Intercept = rep(1, 12)) 
name = names(matrix)
model = vector('list', 5)
R = numeric(5)
AR = numeric(5)
for(p in 1:5){
  matrix = cbind(matrix, (1:12)^p)   #Each time attach a new column
  name = c(name, paste('t', p, sep = ''))  #Append a new column name
  colnames(matrix) = name
  model[[p]] = lm(y ~ . - 1, data = matrix)  
  R[p] = summary(model[[p]])$r.squared
  AR[p] = summary(model[[p]])$adj.r.squared
}
AR
plot(AR)

model21 = model[[2]]
m21 = predict(model21, newdata = data.frame(Intercept = 1, t1 = 13, t2 = 13^2, t3 = 13^3))
m21

#m22
Y = c(one.sided.filters2[[5]][286:288], m13, m14, m15, m16, m17, m18, m19, m20, m21)
matrix = data.frame(y = Y, Intercept = rep(1, 12)) 
name = names(matrix)
model = vector('list', 5)
R = numeric(5)
AR = numeric(5)
for(p in 1:5){
  matrix = cbind(matrix, (1:12)^p)   #Each time attach a new column
  name = c(name, paste('t', p, sep = ''))  #Append a new column name
  colnames(matrix) = name
  model[[p]] = lm(y ~ . - 1, data = matrix)  
  R[p] = summary(model[[p]])$r.squared
  AR[p] = summary(model[[p]])$adj.r.squared
}
AR
plot(AR)

model22 = model[[2]]
m22 = predict(model22, newdata = data.frame(Intercept = 1, t1 = 13, t2 = 13^2, t3 = 13^3))
m22

#m23
Y = c(one.sided.filters2[[5]][287:288], m13, m14, m15, m16, m17, m18, m19, m20, m21, m22)
matrix = data.frame(y = Y, Intercept = rep(1, 12)) 
name = names(matrix)
model = vector('list', 5)
R = numeric(5)
AR = numeric(5)
for(p in 1:5){
  matrix = cbind(matrix, (1:12)^p)   #Each time attach a new column
  name = c(name, paste('t', p, sep = ''))  #Append a new column name
  colnames(matrix) = name
  model[[p]] = lm(y ~ . - 1, data = matrix)  
  R[p] = summary(model[[p]])$r.squared
  AR[p] = summary(model[[p]])$adj.r.squared
}
AR
plot(AR)

model23 = model[[2]]
m23 = predict(model23, newdata = data.frame(Intercept = 1, t1 = 13, t2 = 13^2, t3 = 13^3))
m23

#m24
Y = c(one.sided.filters2[[5]][288], m13, m14, m15, m16, m17, m18, m19, m20, m21, m22, m23)
matrix = data.frame(y = Y, Intercept = rep(1, 12)) 
name = names(matrix)
model = vector('list', 5)
R = numeric(5)
AR = numeric(5)
for(p in 1:5){
  matrix = cbind(matrix, (1:12)^p)   #Each time attach a new column
  name = c(name, paste('t', p, sep = ''))  #Append a new column name
  colnames(matrix) = name
  model[[p]] = lm(y ~ . - 1, data = matrix)  
  R[p] = summary(model[[p]])$r.squared
  AR[p] = summary(model[[p]])$adj.r.squared
}
AR
plot(AR)

model24 = model[[2]]
m24 = predict(model24, newdata = data.frame(Intercept = 1, t1 = 13, t2 = 13^2, t3 = 13^3))
m24

trend.predict = c(m13, m14, m15, m16, m17, m18, m19, m20, m21, m22, m23, m24)
plot(trend.predict)

fit.res = arima (res2, c(5, 0, 4)) 
predict.res = predict(fit.res, n.ahead = 12)        
plot(res2, type = "l", main = 'Residual prediction', ylab = 'Residual')                              
lines(predict.res$pred, col ='red') 

predict = trend.predict + seasonality + predict.res$pred 
predict.t = (1 + lambda*predict)^(1 / lambda)

#True value VS predictted values
plot(UR.predict$VALUE, xaxt = 'n', ylim = c(0,12))
axis(1, at = seq(1, 12, by = 1), labels = UR[,1][seq(289, 300, by = 1)])

points(as.numeric(predict.t), col = 2, ylim = c(0 ,12), pch=19)
