rm(list = ls())

require(tseries)
require(forecast)
require(dplyr)
demo.data = read.csv("C:/Users/Alvin/Desktop/FORECASTS/Biomath Brown Bag/demo_data.csv")

train.length = 148
forecast.length = 15

train.data = demo.data[1: train.length, ]
test.data = demo.data[(train.length+1) : (train.length + forecast.length), ]

a = ccf(train.data$Cases, train.data$X1)
X1.lag = max(forecast.length, a$lag[which(a$acf == max(a$acf))])

a = ccf(train.data$Cases, train.data$X2)
X2.lag = max(forecast.length, a$lag[which(a$acf == max(a$acf))])

a = ccf(train.data$Cases, train.data$X3)
X3.lag = max(forecast.length, a$lag[which(a$acf == max(a$acf))])

a = ccf(train.data$Cases, train.data$X4)
X4.lag = max(forecast.length, a$lag[which(a$acf == max(a$acf))])

#lagged.data = demo.data
#lagged.data$X1.lag = lag(lagged.data$X1, X1.lag)
#lagged.data$X2.lag = lag(lagged.data$X2, X2.lag)
#lagged.data$X3.lag = lag(lagged.data$X3, X3.lag)
#lagged.data$X4.lag = lag(lagged.data$X4, X4.lag)

lagged.data$lag.cases = lag(lagged.data$Cases, 1)

train.lag.data.ml = lagged.data[1:train.length, ]
test.lag.data.ml = lagged.data[(train.length + 1):(train.length + forecast.length), ]

RF.model = randomForest(Cases ~ X1.lag + lag.cases+
                          X2.lag  + X3.lag + X4.lag, 
                        data = na.omit(train.lag.data.ml))


RF.pred = data.frame(matrix(NA, forecast.length, 4))
colnames(RF.pred) = c("point", "lower", "upper", "sd")

for(i in 1 : forecast.length){
  ##' here we alter the lag cases field to be the prediction from the previous week
  if(i > 1){
    test.lag.data.ml$lag.cases[i] = RF.pred[i-1, 1]
  }
  p = predict(RF.model, test.lag.data.ml[i,], predict.all = T)
  RF.pred[i,1] = as.numeric(p$aggregate)
  RF.pred[i,2] = quantile(p$individual, 0.025)
  RF.pred[i,3] = quantile(p$individual, 0.975)
  RF.pred[i,4] = sd(p$individual)
}

plot(c(tail(train.lag.data$Cases,20), test.lag.data$Cases), type = 'o', ylim = c(0, max(upper.arima, upper.lm.lag)), pch = 16, 
     ylab = 'cases', bty = 'n')

points(21:(20+forecast.length), RF.pred$point, col = 'purple', pch = 16, type = 'o')
polygon(x = c(21:(20+forecast.length),rev(21:(20+forecast.length))),y = c(RF.pred$lower, rev(RF.pred$upper)),
        col=adjustcolor('purple', 0.25), border=F)
