# getting lat long for basins
# https://community.intuiface.com/t/export-google-earth-pro-longitude-latitude-into-excel-asset/702
grep -E 'name|coordinates' CBs_sampled_2021.kml > coordinates.txt

library("ggmap")
sites <- read.table("metadata.txt", sep = "\t", header = TRUE)
sites <- subset(sites,sample_control=="sample")
ll_means <- sapply(sites[91:92], mean, na.rm=TRUE) #point to lat and lon columns. may be different col number with different iterations of metadata file
ll_means
array(ll_means)
map_og <- get_map(location = ll_means,  maptype = "satellite", source = "google", zoom = 15) 

novariables <- ggmap(map_og) + 
  geom_point(data = sites, color = "red", size = 1)
novariables

cook <- read.table("WQ_Dips_2020_2021.csv",sep=",",header=TRUE)
cook$Date <- as.Date(cook$Date)
cook$DO_orig <- as.numeric(as.character(cook$DO_orig))
cook$DO <- as.numeric(as.character(cook$DO))
cook <- subset(cook, Year=="2021")

cook[,31:32]
ll_means_weekly <- sapply(cook[31:32], mean, na.rm=TRUE)
ll_means_weekly
array(ll_means_weekly)
map_og_weekly <- get_map(location = ll_means_weekly,  maptype = "satellite", source = "google", zoom = 15) 
ggmap(map_og_weekly)

basins_flowgroup <- ggmap(map_og) + 
  geom_jitter(data = subset(cook, Date=="2021-06-11"), #Date=="6/11/21"
             mapping = aes(x=lon, y=lat, color= Flowgroup_coarse, shape=combined_separate),
             size=2, position = position_jitter(width = 0.0005, height = 0.0005)) + 
  scale_color_discrete()+ #name = "Flow group" ) + #labels = c("Combined", "Separated")
  labs(x="Longitude", y="Latitude", shape="Basin type",color="Flow group")
basins_flowgroup

region_point <- ggmap(map_region) + 
  geom_point(data = subset(sites,UCB=="41"),
             color = "red", size = 6, shape=8, alpha=0.7)+
  labs(x="Longitude", y="Latitude")
region_point