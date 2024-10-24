## Fig 1. Locations of collection sites

library(ggmap)

## Fig 1A. Map of study jurisdiction

sites <- read.table("catch-basin-microbes/input-files/metadata.txt", sep = "\t", header = TRUE)
sites <- subset(sites,sample_control=="sample")
ll_means <- sapply(sites[91:92], mean, na.rm=TRUE) #point to lat and lon columns
map_region <- get_map(location = ll_means,  maptype = "terrain", source = "google", zoom = 6) 
region_point <- ggmap(map_region) + 
  geom_point(data = subset(sites, UCB=="41"),
             color = "red", size = 6, shape=8, alpha=0.7)+
  labs(x="Longitude", y="Latitude")
region_point

## Fig 1B. Distribution of basins by basin type and flow group

map_og <- get_map(location = ll_means,  maptype = "satellite", source = "google", zoom = 15) 
cook <- read.table("input-files/WQ_Dips_2021.csv",sep=",",header=TRUE)
basins_flowgroup <- ggmap(map_og) + 
  geom_jitter(data = subset(cook, Date=="2021-06-11"), #Date=="6/11/21"
             mapping = aes(x=lon, y=lat, color= Flowgroup_coarse, shape=combined_separate),
             size=2, position = position_jitter(width = 0.0005, height = 0.0005)) + 
  scale_color_discrete()+ #name = "Flow group" ) + #labels = c("Combined", "Separated")
  labs(x="Longitude", y="Latitude", shape="Basin type",color="Flow group")
basins_flowgroup
