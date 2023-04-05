---
title: "Poster Envr Data"
author: "Madeline Eppley"
date: "4/5/2023"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Set working directory
```{r}
setwd("~/Desktop/RA_Lotterhos/presentations_posters/Data for Poster 2023")
```

Install packages
```{r}
install.packages("sdmpredictors")
install.packages("leaflet")
install.packages("raster")

install.packages("rgeos", repos="http://R-Forge.R-project.org", type="source")
install.packages("rgdal", repos="http://R-Forge.R-project.org", type="source")
library(devtools)
install_github("r-spatial/sf", configure.args = "--with-proj-lib=/usr/local/lib/")

install.packages("sf")
install.packages("maps")
install.packages("rasterVis")
install.packages("ggrepel")
install.packages("colourvalues")
install.packages("ggsci")
```

```{r}
library(data.table)
library(sdmpredictors)
library(leaflet)
library(ggplot2)
library(rgdal)
library(raster)
library(LEA)
library(sf)
library(maps)
library(sp)
library(rasterVis)
library(scales)
library(ggrepel)
library(colourvalues)
library(stats)
```

Load the Bio-ORACLE Marine Data. All layers can be seen with list_layers() function. Store selected layers. 
```{r}
list_layers()

oyster.envr <- load_layers(layercodes = c("BO_salinity", "BO_sstmean", "BO_sstmax"))

salinity <- load_layers("BO_salinity")
sstmean <- load_layers("BO_sstmean")
sstmax <- load_layers("BO_sstmax")

layer_stats("BO_sstmax")
```

Read in Seascape Samples (Oyster Phenotype) Data from CSV file. 
```{r}
library(readr)
metadata<- read_csv("SeascapeSamples - Phenotype-Lat Version 3.20.23.csv")
View(metadata)
```

View list of unique sites and bind to the environmental data
```{r}
sites <- cbind("lat" = metadata$lat, "lon" = metadata$lon)
head(sites)

envr_sites <- data.frame(Name = sites, extract(oyster.envr, sites))
head(envr_sites)
```

Create environment site matrix data
```{r}
envr_sites_matrix<-as.matrix(envr_sites[,-c(1,2)])

envr_sites_matrix_nas<-na.omit(envr_sites_matrix)

write.env(envr_sites_matrix_nas, "sites_environ_matrix.env")
```

Find the latitude and longitude of the Seascape Samples range. Add some buffer room so the furthest sites crop inside of the map. 
```{r}
lon.range <- range(envr_sites$Name.lon)
lat.range <- range(envr_sites$Name.lat)

print(lon.range)
print(lat.range)

range.extent <- extent(-100, -60, 24, 49)
```

Crop the environmental data and plot it with an color scheme that corresponds to ocean surface temperature. This makes a nice map of our environmental data!
```{r}
sst_max.crop<-crop(sstmax, range.extent)

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

plot(sst_max.crop,col=my.colors(1000), axes=FALSE, box=FALSE)
title(cex.sub = 1.25, sub ="Maximum temperature at the sea surface (ºC)")
```

Let's plot the site lat/long points on our environmental data map. To do this, we can use two methods. The first is to plot with GEO files - Raster & Raster Brick. We need to make a data frame with just our lat/long and site codes data so we can assign a CRS to match the CRS of our cropped environmental data. This is so our points on the map will only show up once (not by x N of individuals) and have the same spatial scale. Export the data frame with the same Coordinate Reference System (CRS) as the Raster map file. Plot the coordinates on the map to make sure they align!
```{r}
#Check the CRS of the enviornmental data
crs(sst_max.crop)

df <- data.frame(longitude = c(-89.49761, -79.061842, -94.90361111, -97.174, -80.92203, -81.02173, -90.01914, -96.2233, -75.88222, -97.359167, -81.2957724, -81.363169, -91.70716, -90.911362, -79.8976693, -93.3179, -70.911544, -71.453, -71.437883, -88.31833, -81.43638, -70.83657, -74.913017, -75.516944, -65.910547, -65.9007577, -93.91806, -76.562783, -76.63554, -75.54978, -62.8874544, -76.97669, -69.53203214, -69.91531434), latitude = c(29.64379, 33.5235297, 29.5475, 28.096, 32.01721, 31.98995, 29.42362, 28.6936, 37.87694, 26.559722, 31.4177723, 31.453512, 29.4196, 29.239925, 32.7523094, 29.84499, 43.053746, 41.505, 41.55897, 30.37278, 30.44003, 42.751733, 39.073592, 39.434722, 43.812486, 43.82746, 29.785, 37.063817, 37.14888, 35.79389, 46.1459253, 38.33696, 44.02452236, 43.83136167), site = c("LKF", "OYL", "GLB", "CPB", "FPL", "SKR", "CUB", "TPL", "CUS", "LGM", "SAP", "MER", "NRF", "SLK", "GCV", "CAL", "GBN", "NWR", "DKV", "CPR", "KPF", "SWI", "CSH", "HPC", "SLB", "EEL", "SLT", "WRS", "DWS", "NPS", "PEI", "LCP", "DYC", "ORC"))

#Match our data frame to the same CRS by creating a Spatial Points Data Frame from the coordinates
my_coords <- df

# create a SpatialPointsDataFrame from the coordinates
my_coords_sp <- SpatialPointsDataFrame(coords=my_coords[,c("longitude","latitude")],
                                        data=my_coords,proj4string=CRS("+proj=longlat +datum=WGS84"))

# plot the raster brick with the latitude and longitude points
plot(sst_max.crop,col=my.colors(1000), axes=FALSE, box=FALSE)
title(cex.sub = 1.25, sub ="Maximum temperature at the sea surface (ºC)")
plot(my_coords_sp, add=TRUE, pch=20, col="navy")
text(my_coords_sp$longitude, my_coords_sp$latitude, labels = my_coords_sp$site, pos = 3, cex = 0.75)

```


We can convert data the other way as well. Before, we turned the coordinates into a Raster format. This time, we can turn our raster data into a data frame! Start by converting the Raster Brick into a data frame, then plot it along with the site coordinates in lat/long format in our 'my_coords' data frame. While we do this, let's also assign population colors so that we can see a nice gradient of latitude based on color in our map. 
```{r}
# Convert the RasterBrick object to a data frame
r_df <- as.data.frame(sst_max.crop, xy=TRUE)
coord_df <- as.data.frame(my_coords_sp)

r_df.nas <- na.omit(r_df)

#Assign color based on population
print(sites.unique <- unique(coord_df$site))
print(n.sites <- length(sites.unique))

coord_df <- coord_df[order(coord_df$latitude),]
coord_df$color <- factor(coord_df$site, levels = coord_df$site, labels = heat.colors(34))


# Add the my_coords_sp data to the plot, use the df with NA's removed. Use population colors for the map. 
ggplot() +
  geom_raster(data = r_df.nas, aes(x = x, y = y, fill = BO_sstmax)) +
  scale_fill_gradientn(colours = my.colors(1000), name = "Sea Surface Temperature (ºC)") +
  geom_point(data = coord_df, aes(x = longitude, y = latitude), color = coord_df$color, size = 4) +
  geom_label(data = coord_df, aes(x = longitude, y = latitude), label = coord_df$site, size = 4, hjust = 0, vjust = 0, fill = coord_df$color) +
  coord_equal() + 
  theme_classic() +
  theme(axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "black")) + 
  labs(x = NULL, y = NULL, 
    title = "Seascape Genomics Sampling Range",
    subtitle = "Eastern Oyster (C. virginica)", 
    caption = "Marine Environmental Data: Bio-ORACLE, 2014") + 
  theme(plot.title = element_text(hjust = 0.9), plot.subtitle = element_text(hjust = 0.9), plot.caption = element_text(hjust = 1)) + 
  theme(legend.position = "bottom", legend.background = element_rect(fill = "black")) + 
  theme(axis.text = element_text(color = "white"),
        axis.title = element_text(color = "white"),
        plot.caption = element_text(color = "white"),
        plot.subtitle = element_text(color = "white"),
        plot.title = element_text(color = "white"))

```

Analyze the phenotypic metadata. Start with creating a data matrix of only the numeric values. This will later include disease risk scores and month of collection. 
```{r}
# Create a matrix of the numeric variables
df_matrix <- data.matrix(metadata[, c("lat", "lon", "Shell_length_cm", "Shell_width_cm")])

#normalize
data.norm <- scale(df_matrix)
head(data.norm)

# Run the PCA
pca_results <- princomp(data.norm)


# Print the results
print(pca_results)
summary(pca_results)
plot(pca_results)

#Geographic position explains the most amount of variance, then shell length and width
pca_results$loadings[, 1:2]

#Save the scores and check the names of the PCA loadings to label the graph
pcaplot <- as.data.frame(pca_results$scores)
names(pcaplot)


# Plot the PCA results
ggplot(data = pcaplot, aes(x = Comp.1, y = Comp.3)) +
  geom_point(color = metadata$Color) +
  theme(axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "black")) +
  theme(axis.text = element_text(color = "white"),
        axis.title = element_text(color = "white"),
        plot.title = element_text(color = "white")) +
  labs(x = "PC1 (50.9%)", y = "PC3 (9.4%)", 
    title = "PCA of Phenotype Data")
  

```

Plot shell length by latitude relationship
```{r}
#Shell length and shell width by latitude
plot(metadata$Shell_width_cm ~ metadata$lat)
plot(metadata$Shell_length_cm ~ metadata$lat)

ggplot(data = metadata, aes(x = Shell_length_cm, y = lat)) +
  geom_point(color = metadata$Color) +
  xlab("Shell Length (cm)") + 
  ylab("Latitude") + 
  ggtitle("Shell Length by Latitude") +
  theme(panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "black"), panel.grid = element_blank(), 
        axis.text = element_text(color = "white"),
        axis.title = element_text(color = "white"),
        plot.title = element_text(color = "white")) 

```

Plot various other relationships to interpret metadata
```{r}
plot(metadata$Pea_crab ~ metadata$lat)
```

```{r}
plot(metadata$Ripe_Gonad ~ metadata$lat)
```

```{r}
plot(metadata$Boring_sponge ~ metadata$lat)
```
```{r}
plot(metadata$Ripe_Gonad ~ metadata$Shell_length_cm)
```

```{r}
#Pea crab and shell length relationship
plot(metadata$Pea_crab ~ metadata$Shell_length_cm)
```


```{r}
#y-axis male = 1, female = 0
metadata.nas<- metadata[!is.na(metadata$Sex),]
plot(metadata.nas$Sex ~ metadata.nas$lat)
```


```{r}
plot(metadata$Shell_length_cm ~ metadata$Month)

```

```{r}
plot(metadata$Sex ~ metadata$Shell_length_cm)
```

```{r}
boxplot(Shell_length_cm~Sex,data=metadata, main="Sex by Length",
   xlab="Sex, M=1 F=0", ylab="Length in cm")

```


```{r}
boxplot(Shell_length_cm ~ Pea_crab*Boring_sponge, data=metadata, notch=TRUE,
  col=(c("gold","darkgreen")),
  main="Length of Oyster as an indicator of health with parasites", xlab="Pea Crab and Boring Sponge Prevalence")


```

```{r}
#Summary of all metadata
summary(metadata)

#Check for normality

#Normal distribution
hist(metadata$Shell_length_cm)
shapiro.test(metadata$Shell_length_cm)

log.shell <- log10(metadata$Shell_length_cm)
hist(log.shell)
shapiro.test(log.shell)

#Bimodal
hist(metadata$lat)
shapiro.test(metadata$lat)

log.lat <- log10(metadata$lat)
hist(log.lat)
shapiro.test(log.lat)

```

```{r}
metadata.m <- as.matrix(metadata, rownames.force = NA)

metadata.df <- data.frame(metadata)

lm.lengthbyheight = lm(metadata.df$Shell_length_cm ~ metadata.df$lat)

plot(Shell_length_cm ~ lat, data=metadata)                     
summary(lm.lengthbyheight)

par(mfrow=c(2,2))
plot(lm.lengthbyheight)


# With log transformed data
lm.log.lengthbyheight = lm(log.shell ~ log.lat)
summary(lm.log.lengthbyheight)
par(mfrow=c(2,2))


plot(lm.log.lengthbyheight)


```

Create a Disease Risk Index by writing a for loop to add points to our data frame based on the amount of disease present in the system. The risk score values are assigned based on severity of disease, resultant mortality, and co-infection. The boring sponge is the most severe disease that can cause the highest mortality of the macrofauna, so we will assign it a value of 1.5. The pea crab will be assigned a value of 1, resulting in less mortality than the boring sponge. Hinge worms will be assigned a value of 0.5. If any individual is infected with more than one disease, it will receive a co-infection score of 0.5 added to the overall score.  
```{r}
#Let's start with a test data frame to show the functionality of the for loop. 
metadata.test <- data.frame(
  Label_ind = c("Sample1", "Sample2", "Sample3"),
  Pea_crab = c(TRUE, FALSE, TRUE),
  Boring_sponge = c(FALSE, TRUE, TRUE),
  riskscore = 0
)

# Loop through each row in the dataframe and add to the riskscore column based on the conditions that are present.
for (i in 1:nrow(metadata.test)) {
  if (metadata.test$Pea_crab[i]) {
    metadata.test$riskscore[i] <- metadata.test$riskscore[i] + 1
  }
  if (metadata.test$Boring_sponge[i]) {
    metadata.test$riskscore[i] <- metadata.test$riskscore[i] + 1.5
  }
}

# View the updated dataframe
metadata.test
```

Now for the real data frame:
```{r}
metadata$riskscore <- c(0)

# Loop through each row in the dataframe and add to the riskscore column based on the conditions that are present.
for (i in 1:nrow(metadata)) {
  if (metadata$Pea_crab[i]) {
    metadata$riskscore[i] <- metadata$riskscore[i] + 1
  }
  if (metadata$Boring_sponge[i]) {
    metadata$riskscore[i] <- metadata$riskscore[i] + 1.5
  }
}

# View the updated dataframe
metadata

#Plot the results
plot(metadata$lat, metadata$riskscore)

# Plot the results with the color scheme
ggplot(data = metadata, aes(x = lat, y = riskscore)) +
  geom_point(color = metadata$Color) +
 #geom_count()
  xlab("Latitude") + 
  ylab("Disease Risk") + 
  ggtitle("Disease Risk Score by Latitude") +
  theme(panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "black"), panel.grid = element_blank(), 
        axis.text = element_text(color = "white"),
        axis.title = element_text(color = "white"),
        plot.title = element_text(color = "white")) 
```

