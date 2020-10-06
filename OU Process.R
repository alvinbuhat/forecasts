setwd("~/Downloads")

library(lubridate)

case <- read.csv("Daily and Cumulative Confirmed Cases from DOH COVID Data Drop_ 20200731 - 04 Case Information.csv")
date=case$Date[1:148]
Xd=case$Daily[1:148]
Xc=case$Cumulative[1:148]

plot(Xc, type="l")
plot(Xd, type="l")


# Parameter estimation and bias correction for diffusion processes by Tang and
n=length(Xd)-1
summ1=0
for (i in 2:(n+1)){
  summ1 = summ1+ Xd[i]*Xd[i-1]
}
gamma1 = ((1/n)*summ1-(1/n^2)*sum(Xd[2:(n+1)])*sum(Xd[1:n]))/((1/n)*sum(Xd[1:n]^2)-(1/n^2)*(sum(Xd[1:n]))^2)

summ2=0
for (i in 2:(n+1)){
  summ2 = summ2 + Xd[i]-gamma1*Xd[i-1]
}
gamma2 = n^(-1)*summ2/(1-gamma1)

summ3 = 0
for (i in 2:(n+1)){
  summ3 = summ3 + (Xd[i]-gamma1*Xd[i-1]-gamma2*(1-gamma1))^2
}
gamma3 = n^(-1)*summ3

alpha_est = -log(gamma1)
print(alpha_est)
beta_est = gamma2
print(beta_est)
sigma_est = sqrt(2*alpha_est*gamma3/(1-gamma1^2))
print(sigma_est)

# Y are the forecast values
Yd=rep(0, 15)
Yd[1]=Xd[148]*exp(-alpha_est) + beta_est*(1-exp(-alpha_est))
for (i in 2:15){
  Yd[i] = Yd[i-1]*exp(-alpha_est) + beta_est*(1-exp(-alpha_est))
}

Yd[1]=sum(Xd)+Yd[1]
for (i in 2:15){
  Yd[i] = Yd[i-1] + Yd[i]
}

# Parameter estimation and bias correction for diffusion processes by Tang and
n=length(Xc)-1
summ1=0
for (i in 2:(n+1)){
  summ1 = summ1+ Xc[i]*Xc[i-1]
}
gamma1 = ((1/n)*summ1-(1/n^2)*sum(Xc[2:(n+1)])*sum(Xc[1:n]))/((1/n)*sum(Xc[1:n]^2)-(1/n^2)*(sum(Xc[1:n]))^2)

summ2=0
for (i in 2:(n+1)){
  summ2 = summ2 + Xc[i]-gamma1*Xc[i-1]
}
gamma2 = n^(-1)*summ2/(1-gamma1)

summ3 = 0
for (i in 2:(n+1)){
  summ3 = summ3 + (Xc[i]-gamma1*Xc[i-1]-gamma2*(1-gamma1))^2
}
gamma3 = n^(-1)*summ3

alpha_est = -log(gamma1)
print(alpha_est)
beta_est = gamma2
print(beta_est)
sigma_est = sqrt(2*alpha_est*gamma3/(1-gamma1^2))
print(sigma_est)

# Y are the forecast values
Yc=rep(0, 15)
Yc[1]=Xc[148]*exp(-alpha_est) + beta_est*(1-exp(-alpha_est))
for (i in 2:15){
  Yc[i] = Yc[i-1]*exp(-alpha_est) + beta_est*(1-exp(-alpha_est))
}

dates=rep(0,15)
dates[1] = as.Date("2020-08-01")
for (i in 2:15){
  dates[i] = dates[i-1]+days(1)
}
Actual = rep(0, 15)
  Actual[1] = 97421
  Actual[2] = 102348
  Actual[3] = 105503
  Actual[4] = 111490
  Actual[5] = 114913
  Actual[6] = 118409
  Actual[7] = 121733
  Actual[8] = 125907
  Actual[9] = 128948
  Actual[10] = 135841
  Actual[11] = 138761
  Actual[12] = 143124
  Actual[13] = 147063
  Actual[14] = 153219
  Actual[15] = 157526

plot(Actual, type="l", ylab="cumulative cases")
lines(Yd, type="l", col="red")
lines(Yc, type="l", col="blue")

M=length(Xc)+15
X=rep(0,M)
X[1:148]=Xc
X[149:163]=Actual
XX=rep(0,M)
XX[1:148]=rep(NA,148)
XX[149:163]=Yd
dates = seq(as.Date("2020/3/6"), as.Date("2020/8/15"), "days")
dat = data.frame(dates,X,XX)
plot(dates, X, type="l", xlab="date", ylab="cumulative cases", lwd=2)
lines(dates, XX, type="l", col="red", lwd=2, lty=2)
legend("topleft", legend=c("actual values", "forecast values"), col=c("black", "red"), lty=1:2, cex=0.8, bg="lightblue")
