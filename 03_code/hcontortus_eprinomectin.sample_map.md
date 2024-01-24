# hcontortus eprinomectin analysis: maps of sampling sites

### author: Stephen Doyle

- want to generate some maps to show the sampling sites, and importantly get a sense of the relative distance between sites
- I think this would be useful for understanding shared variation between sampling sites and whether resistance is evolving independently or spreading




## Map of France
- found some shape files that might be useful of regions of France 
    - https://www.data.gouv.fr/fr/datasets/contours-des-departements-francais-issus-d-openstreetmap/#/community-resources   
    - "Map of French departments (simplified Rhône dept.)"
    - has some coditions of use, which might be important for publishing
        - This data comes from crowdsourcing carried out by contributors to the OpenStreetMap project and is under ODbL license which requires identical sharing and the obligatory attribution statement must be " © OpenStreetMap contributors under ODbL license " in accordance with http: //osm.org/copyright

```R
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

```


## Zoomed in view of the sampling sites
- want a zoomed in view, which will better show distance between sites, and perhaps some landscape information
- tried a few different things, before getting some help from Shannon S using OpenStreetMap and "mapview" package

```R
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

```