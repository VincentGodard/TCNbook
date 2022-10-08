library("terra")
library("TCNtools")

bas = vect("data/gis/marsyandi.gpkg")
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
# plot ----
for (name in c("Khudi","Chudi","Marsyandi")){
dem = rast(paste0("data/gis/",name,".tif"))
slope <- terra::terrain(dem, "slope", unit="radians")
aspect <- terra::terrain(dem, "aspect", unit="radians")
hill <- terra::shade(slope, aspect, 40, 270)
plot(hill, col=grey(0:100/100), legend=FALSE,axes=F,buffer=T,main=name)
plot(dem, col=rainbow(25, alpha=0.35), add=TRUE,axes=F,plg=list(cex=0.5),legend="bottom")
sbar(label="km")
lines(bas)
#
dem = project(dem,"epsg:4326")
data =  as.data.frame(dem,xy=T)
colnames(data)[3]<-"z"
write.table(data,paste0("data/gis/",name,".dat"),col.names = T,row.names = F)
}

# importing the data
#path_to_data = "data/gis/"
path_to_data ="https://raw.githubusercontent.com/VincentGodard/TCNbook/main/data/gis/"
names_basins = c("Marsyandi","Chudi","Khudi")
cols = c("deepskyblue", "khaki3", "darkolivegreen3")
basins = list()
for (i in 1:length(names_basins)){
  basins[[i]] = read.table(paste0(path_to_data,names_basins[i],".dat"),header=T)
}
  
# plotting
for (i in 1:length(names_basins)){
  if (i == 1){ # we initiate a plot at the first iteration
    plot(density(basins[[i]]$z),ylim=c(0,0.0022),xlab="Elevation (m)",col=cols[i],main="")
  }else{
    lines(density(basins[[i]]$z),col=cols[i])
  }
}
legend("topright",names_basins,cex=0.5,lty=1,col=cols)

# computing scaling pixel by pixel
for (i in 1:length(names_basins)){
basins[[i]]$P = atm_pressure(alt=basins[[i]]$z,model="stone2000") 
st = scaling_st(basins[[i]]$P,basins[[i]]$y)
basins[[i]] = cbind(basins[[i]],st)
}

# plotting result
for (i in 1:length(names_basins)){
  if (i == 1){ # we initiate a plot at the first iteration
    plot(density(basins[[i]]$Nneutrons),ylim=c(1e-4,2.5),xlab="Spallation scaling factor",col=cols[i],main="",log="y")
  }else{
    lines(density(basins[[i]]$Nneutrons),col=cols[i])
  }
}
legend("topright",names_basins,cex=0.5,lty=1,col=cols)

# computing averages values
avg_st = data.frame(names=names_basins)
avg_st$S1 = NA
avg_st$S2 = NA
for (i in 1:nrow(avg_st)){
  avg_st$S1[i] = mean(basins[[i]]$Nneutrons)
  P2 = atm_pressure(alt=mean(basins[[i]]$z),model="stone2000") 
  st = scaling_st(P2,mean(basins[[i]]$y))
  avg_st$S2[i] = st$Nneutrons
}
print(avg_st)

plot(avg_st$S2,avg_st$S1,pch=21,bg=cols,
     xlab="S2 : Scaling from average elevation",ylab="S1 : Scaling pixel by pixel",
     xlim=range(0,avg_st$S2,avg_st$S1),ylim=range(0,avg_st$S2,avg_st$S1))
grid()
text(avg_st$S2,avg_st$S1,avg_st$names,cex=0.5,pos=1)
abline(0,1,lty=2)


test = download.file("https://github.com/VincentGodard/TCNbook/blob/main/data/gis/Chudi.dat",destfile="/tmp/test.dat")
test = read.table("https://raw.githubusercontent.com/VincentGodard/TCNbook/main/data/gis/Chudi.dat",header=T)
