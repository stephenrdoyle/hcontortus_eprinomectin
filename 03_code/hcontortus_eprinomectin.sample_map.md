



https://www.data.gouv.fr/fr/datasets/contours-des-departements-francais-issus-d-openstreetmap/#/community-resources
"Map of French departments (simplified Rhône dept.)"

This data comes from crowdsourcing carried out by contributors to the OpenStreetMap project and is under ODbL license which requires identical sharing and the obligatory attribution statement must be " © OpenStreetMap contributors under ODbL license " in accordance with http: //osm.org/copyrigh


library(ggplot2)
library(dplyr)
library(rgdal)

france.shp <- readOGR(dsn = "map_fr_dept_remaked.shx.shp", stringsAsFactors = F)

highlight.map <- france.shp[france.shp$nom == "Pyrénées-Atlantiques",]

sample_data <- read.table("map_data.txt", header=T, sep="\t")

ggplot() + 
    geom_polygon(data = france.shp, aes(x = long, y = lat, group = group), colour = "white", fill = "grey") + 
    coord_map(xlim=c(-5.2,8),  ylim=c(42,51)) + 
    theme_void() + 
    geom_polygon(data = highlight.map, aes(x = long, y = lat, group = group), colour = "white", fill = "cornflowerblue") +
    geom_point(data=sample_data, aes(longitude, latitude), size=1)

ggsave("france_sample_map.pdf")

wget https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/raster/HYP_HR_SR.zip
unzip HYP_HR_SR.zip


https://eriqande.github.io/rep-res-eeb-2017/plotting-spatial-data-with-ggplot.html

library(ggplot2)
library(dplyr)
library("raster")
library(sp)
library(ggspatial)
library(rgdal)
library(ggrepel)






library(mapview)
library(sf)

sample_data <- read.table("map_data.txt", header=T, sep="\t")


#creat shapefiles
data_shapefile <- st_as_sf(sample_data,
    coords = c("longitude", "latitude"),
    crs= "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# Render map
map <- mapview(data_shapefile,
                            map.types = "OpenStreetMap",
                            zcol = "eprinomectin_response",
                            layer.name = "eprinomectin_response", # 'layer.name' is legend title
                            use.layer.names = "eprinomectin_response",
                            legend = TRUE,
                            col.regions=list("red", "blue"))
#view map 
mapshot(map, file="france_sample_map_zoom.pdf", remove_controls = c("zoomControl", "layersControl", "homeButton","drawToolbar", "easyButton"))