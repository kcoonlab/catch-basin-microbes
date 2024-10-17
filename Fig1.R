library("ggmap")

sites <- read.table("metadata.txt", sep = "\t", header = TRUE)
sites <- subset(sites,sample_control=="sample")
ll_means <- sapply(sites[91:92], mean, na.rm=TRUE) #point to lat and lon columns. may be different col number with different iterations of metadata file
map_og <- get_map(location = ll_means,  maptype = "satellite", source = "google", zoom = 15) 

map_region <- get_map(location = ll_means,  maptype = "terrain", source = "google", zoom = 6) 
region_point <- ggmap(map_region) + 
  geom_point(data = subset(sites,UCB=="41"),
             color = "red", size = 6, shape=8, alpha=0.7)+
  labs(x="Longitude", y="Latitude")
region_point

cook <- read.table("WQ_Dips_2021.csv",sep=",",header=TRUE)

basins_flowgroup <- ggmap(map_og) + 
  geom_jitter(data = subset(cook, Date=="2021-06-11"), #Date=="6/11/21"
             mapping = aes(x=lon, y=lat, color= Flowgroup_coarse, shape=combined_separate),
             size=2, position = position_jitter(width = 0.0005, height = 0.0005)) + 
  scale_color_discrete()+ #name = "Flow group" ) + #labels = c("Combined", "Separated")
  labs(x="Longitude", y="Latitude", shape="Basin type",color="Flow group")
basins_flowgroup
