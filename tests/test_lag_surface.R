library("TCNtools")
library("pracma")


data(Lambda) # we load a vector containing the attenuation length into the environment
rho = 2.75 # we also define the density (g/cm3)
data(prm) # we load a matrix containing the production/decay parameters into the environment

altitude = 2252 # elevation in m
latitude = 45.97 # latitude in degrees
longitude = 6.96 # longitude
A = 0.939*exp(-2.6*rho/160)
#P = atm_pressure(alt=altitude,model="stone2000") # compute atmospheric pressure at site
P = atm_pressure(alt=altitude,lat=latitude,lon=longitude,model="era40") # compute atmospheric pressure at site
S = scaling_st(P,latitude) # compute the scaling parameters according to Stone (2000)

data = data.frame(t1=seq(0,50e3,length.out=2000)) # 
data$t2 = max(data$t1) - data$t1 # time BP
data$vdm = get_vdm(data$t2,model="lsd")
data$rc = vdm2rc(data$vdm,latitude,model="elsasser54")
data$lm = scaling_lm(P,data$rc)

plot(data$t2,data$lm,type="l")
abline(h=S$Nneutrons,col="red")
abline(h=mean(data$lm),lty=2)


prm2 = prm[,"Be10"]
prm2[1]=4.11
data$C1 = solv_conc_eul(0,0,data$t1,0,prm2,S*A,Lambda) # compute concentration eulerian
#data$Cobs1 = max(data$C1) - data$C1
data$C2 = cumtrapz(data$t1,(prm2[1]*data$lm*A+(prm[2,"Be10"] + prm[3,"Be10"])*A*S$Nmuons)*(exp(-prm["lambda",'Be10']*data$t1)))[,1] 
#       + (prm[2,"Be10"] + prm[3,"Be10"])*A*S$Nmuons*(exp(-prm["lambda",'Be10']*data$t1))*data$t1 

nr=nrow(data)
data$C2obs =  NA
data$C1obs = NA
for (i in 1:(nr-1)){
  data$C1obs[i] = solv_conc_eul(0,0,data$t1[1:(nr-i)],0,prm2,S*A,Lambda)[nr-i] # compute concentration eulerian
  data$C2obs[i] = cumtrapz(data$t1[i:nr],(prm2[1]*data$lm[i:nr]*A+(prm[2,"Be10"] + prm[3,"Be10"])*A*S$Nmuons)*(exp(-prm["lambda",'Be10']*data$t1[i:nr])))[nr-i,1] 
  
}

# lines(data$t1,data$C2,col="red")

# 
data$Smu  = rep(as.numeric(S[2]),nrow(data)) # scaling muons
data$z = 0 # we stay at the surface
data$Psp10 = prm2[1]*exp(-1*data$z/Lambda["Lspal"]) # spallation
data$Pmu10 = prm["Pstop",'Be10']*exp(-1*data$z/Lambda["Lstop"]) + prm["Pfast",'Be10']*exp(-1*data$z/Lambda["Lfast"]) # muons
data$C3 = solv_conc_lag(data$t1,data$z,0,data$Psp10,data$Pmu10,prm["lambda",'Be10'],cbind(data$lm,data$Smu)*A,final=FALSE)
#data$Cobs2 = max(data$C2) - data$C2
test = solv_conc_lag(data$t1,data$z,0,data$Psp10,data$Pmu10,prm["lambda",'Be10'],cbind(data$lm,data$Smu)*A,final=TRUE)

plot(data$t1,data$C1,type="l")
lines(data$t1,data$C2,col="green")
lines(data$t1,data$C3,col="blue")


plot(data$t2,data$C1obs,type="l")
lines(data$t2,data$C2obs,col="green")

# 
# compute age
res = data.frame(Cobs=rnorm(n=200e3,26.8e4,1.3e4))
res$A1 = approx(data$C1obs,data$t2,res$Cobs)$y
plot(density(res$A1))
res$A2 = approx(data$C2obs,data$t2,res$Cobs)$y
lines(density(res$A2),col="green")
















library(pracma) # useful library containing the cumtrapz function for trapezoidal integration
Pspal = prm[1,nuc]*S$Nneutrons*Ss # scaled spallation production rate in at/g/y (st scaling)
Pstop = prm[2,nuc]*S$Nmuons*Ss # scaled stopped muons  production rate in at/g/y
Pfast = prm[3,nuc]*S$Nmuons*Ss # scaled fast muons  production rate in at/g/y
lambda = prm[4,nuc] # radioactive decay (1/y)
data$Prod_st = Pspal + Pstop + Pfast

data$C1 = solv_conc_eul(0,0,data$t1,0,prm,S*Ss,Lambda) 
data$C1b = cumtrapz(data$t1,data$Prod_st*exp(-lambda*data$t1))
data$C1b = exp(-lambda*data$t1)*cumtrapz(data$t1,data$Prod_st*exp(lambda*data$t1))

data$Prod_lm = prm[1,"Be10"]*data$lm*Ss+(prm[2,"Be10"] + prm[3,"Be10"])*Ss*S$Nmuons
data$C2 = exp(-lambda*data$t1)*cumtrapz(data$t1,data$Prod_lm*exp(lambda*data$t1))
data$C2b = cumtrapz(data$t1,data$Prod_lm*exp(-lambda*data$t1))
data$C2c = NA
for (i in 2:nrow(data)){
  data$C2c[i] = trapz(data$t1[1:i],data$Prod_lm[1:i]*exp(-lambda*(data$t1[i]-data$t1[1:i])))
}

#data$C2 = cumtrapz(data$t1,data$Prod_lm*(exp(-prm["lambda",'Be10']*(data$t1))))[,1] 
plot(data$t1,data$C1,type="l",xlab="Time since start of exposure (a)",ylab="Concentration (at/g)",col="cyan4")
lines(data$t1,data$C2,col="darkgoldenrod2")
lines(data$t1,data$C1b,col="red",lty=2)
lines(data$t1,data$C2b,col="pink",lty=2)
lines(data$t1,data$C2c,col="green",lty=2)

legend("topleft",c("st","lm"),lty=1,col=c("cyan4", "darkgoldenrod2"))