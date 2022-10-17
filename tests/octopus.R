library("terra")
library("TCNtools")

data(Lambda)
data(prm)
rho = 2.65

octopus = vect("/home/vincent/data/chantiers/monde/vecteurs/octopus_v2/octopus_v2.gpkg")
octopus = as.data.frame(octopus)
octopus = octopus[octopus$ELEV_AVE>0 & 
                  octopus$AREA<100 & octopus$AREA>20 &
                  octopus$EBE_MMKYR<2000 & octopus$EBE_MMKYR>0  ,]
write.table(octopus,"data/gis/octopus.dat",col.names = T,row.names = F)



octopus$P = atm_pressure(alt=octopus$ELEV_AVE,model="stone2000")
S = scaling_st(octopus$P,octopus$Y_WGS84)
octopus = cbind(octopus,S)

Pspal = prm["Pspal","Be10"]*octopus$Nneutrons
Pstop = prm["Pstop","Be10"]*octopus$Nmuons
Pfast = prm["Pfast","Be10"]*octopus$Nmuons
lambda = prm["lambda","Be10"]
C = octopus$BE10NC

octopus$E1 = (Pspal*Lambda[1])/(rho*C)*10*1000

plot(octopus$EBE_MMKYR,octopus$E1,pch=21,bg="lightblue",cex=0.5,
     xlab="Octopus denudation rate (mm/ka)",ylab="Recalculated denudation rate (mm/ka)")
abline(0,1,col="red",lty=2)
ER1 = mean(((octopus$EBE_MMKYR-octopus$E1)/octopus$EBE_MMKYR)^2)

octopus$E2 = (Pspal*Lambda[1] + Pstop*Lambda[2] + Pfast*Lambda[3])/(rho*C)*10*1000
plot(octopus$EBE_MMKYR,octopus$E2,pch=21,bg="lightblue",cex=0.5)
abline(0,1,col="red",lty=2)
ER2 = mean(((octopus$EBE_MMKYR-octopus$E2)/octopus$EBE_MMKYR)^2)


plot(octopus$EBE_MMKYR,octopus$E2,log="xy",pch=21,bg="lightblue",cex=0.5)
abline(0,1,col="red",lty=2)


compute_C <- function(ero,prm,S,Lambda){
  C = solv_conc_eul(0,ero,Inf,0,prm,S,Lambda) # compute concentration
  return(C)
}
fun_opt <-function(ero,prm,S,Lambda,Cmes){
  Cmod =  compute_C(ero,prm,S,Lambda)
  return(abs(Cmod-Cmes))
}

octopus$E3 = NA
for (i in 1:nrow(octopus)){
  print(i)
  res = optimize(fun_opt,c(0,3000)/10*rho,prm[,"Be10"],c(octopus$Nneutrons[i],octopus$Nmuons[i]),
                 Lambda,octopus$BE10NC[i])
  octopus$E3[i] = res$minimum*10/rho*1000
}

plot(octopus$EBE_MMKYR,octopus$E3,log="xy",pch=21,bg="lightblue",cex=0.5)
abline(0,1,col="red",lty=2)
ER3 = mean(((octopus$EBE_MMKYR-octopus$E3)/octopus$EBE_MMKYR)^2)


# library(sf)
# library(leaflet)
# 
# m <- leaflet() %>%
#   addTiles() %>%
#   addProviderTiles("OpenStreetMap",group = "OpenStreetMap") %>%
#   addProviderTiles("Stamen.Terrain",group = "Stamen.Terrain") %>%
#   addProviderTiles("Esri.WorldImagery",group = "Esri.WorldImagery") %>%
#   addLayersControl(baseGroups = c("OpenStreetMap","Stamen.Terrain", "Esri.WorldImagery"),position = "topleft") %>%
#   addMarkers(lng = octopus$X_WGS84, lat = octopus$Y_WGS84)
# 


