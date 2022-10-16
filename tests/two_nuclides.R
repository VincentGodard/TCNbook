library("TCNtools")


data(Lambda) # we load a vector containing the attenuation length into the environment
data(prm) # we load a matrix containing the production/decay parameters into the environment
rho = 2.7 # we also define the density (g/cm3)
altitude = 2000 # elevation in m
latitude = 43 # latitude in degrees
P = atm_pressure(alt=altitude,model="stone2000") # compute atmospheric pressure at site
S = scaling_st(P,latitude) # compute the scaling parameters according to Stone (2000)


N1 = "Be10" # longer half-life
N2 = "Al26" # shorter half-life

data = data.frame(t = seq(0,7e6,by=2000))
ero = 50 * rho * 100 / 1e6 # in g/cm2/a

C1_0 = solv_conc_eul(0,ero,Inf,0,prm[,N1],S,Lambda)
C2_0 = solv_conc_eul(0,ero,Inf,0,prm[,N2],S,Lambda)

data$C1 = solv_conc_eul(Inf,0,data$t,C1_0,prm[,N1],S,Lambda)
data$C2 = solv_conc_eul(Inf,0,data$t,C2_0,prm[,N2],S,Lambda)
plot(data$t,data$C2/data$C1,type="l",lwd=2,col="darkorange",xlab="Time (a)",ylab=paste(N2,"/",N1))

tmp = tnp_curves(prm[,N1],prm[,N2],Lambda,S,rho,dlim = c(0.01, 1000),alim = c(1000, 1e+09))
ss_ero = tmp[[1]]
cst_exp = tmp[[2]]

plot(NA,xlim=range(cst_exp$C1,ss_ero$C1),ylim=range(1,cst_exp$C2/cst_exp$C1,ss_ero$C2/ss_ero$C1),log="x",
     xlab=paste(N1,"(at/g)"),ylab=paste(N2,"/",N1))
lines(cst_exp$C1,cst_exp$C2/cst_exp$C1,lty=2,col="khaki4") # constant exposure, dashed line
lines(ss_ero$C1,ss_ero$C2/ss_ero$C1,col="khaki4") # steady-state erosion, solid line
lines(data$C1,data$C2/data$C1,lwd=2,col="darkorange")

A = (1:10)*1e6 # increments in age (a)
E = c(1,2,5,10,20,50,100,200,500,1000,2000) # increments in denudation rate (m/Ma)
res = tnp_burial(A,E,prm[,N1],prm[,N2],Lambda,S,rho,n=100) # compute array

plot(NA,xlim=range(1e4,2e5),ylim=range(0.01,cst_exp$C2/cst_exp$C1,ss_ero$C2/ss_ero$C1),log="x",
     xlab=paste(N1,"(at/g)"),ylab=paste(N2,"/",N1))
lines(res[[1]]$C1,res[[1]]$C2/res[[1]]$C1,col="grey") # constant age
lines(res[[2]]$C1,res[[2]]$C2/res[[2]]$C1,col="black") # constant denudation
lines(cst_exp$C1,cst_exp$C2/cst_exp$C1,lty=2,col="khaki4") # constant exposure, dashed line
lines(ss_ero$C1,ss_ero$C2/ss_ero$C1,col="khaki4") # steady-state erosion, solid line

lines(data$C1,data$C2/data$C1,lwd=2,col="darkorange")

data = read.table("data/sartegou2020late.dat",header=T)

el = tnp_ellipse(data$C10*1e3, data$C10_e*1e3, data$C26*1e3, data$C26_e*1e3,confidence=0.68)
for (i in 1:length(el)) {polygon(el[[i]],col="pink",border="grey")}
text(data$C10*1e3,data$C26/data$C10,data$num,cex=0.5)

# C10 = 13
# C26 = 38
# points(C10*1000,C26/C10)

