library(sf)
library(leaflet)
df = data.frame(X =  6.96065,Y = 45.96670,label="ARG-16-9")
  
m <- leaflet(df) %>%
  addTiles() %>%
  addProviderTiles("OpenStreetMap",group = "OpenStreetMap") %>%
  addProviderTiles("Stamen.Terrain",group = "Stamen.Terrain") %>%
  addProviderTiles("Esri.WorldImagery",group = "Esri.WorldImagery") %>%
  addLayersControl(baseGroups = c("OpenStreetMap","Stamen.Terrain", "Esri.WorldImagery"),position = "topleft") %>%
  addLabelOnlyMarkers(data = df,
                      lng = ~X, lat = ~Y, label = ~label, 
                      labelOptions = labelOptions(noHide = TRUE, direction = 'auto', textOnly = FALSE)) 
m