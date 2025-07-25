

### projections and distribution maps for abiotic variables (from ESMs)
### extract, bias correct and plot ESM data
### need to bring in all .nc files
### produces values, dfs, etc used in Plankton SDMs file


library(SDMTools)
library(ncdf4)
library(ggOceanMaps)
library(raster)
library(cowplot)
library(sf)
library(ggplot2)
library(tmap)
library(tmaptools)
library(marmap)
library(rnaturalearth)
library(maps)
library(writexl)
library(rnaturalearthdata)
library(scales)
library(purrr)
library(dplyr)
library(RColorBrewer)
library(readxl)
library(gridExtra)
library(grid)

#get country shape file ready 
country <- ne_countries(scale = "medium", returnclass = "sf")

#bring in data for bias correction
historicalclimatedata <- read_excel("historicalclimatedata.xlsx")

#### SALINITY 

#sst historical sal.. set wd to desktop to get data in r 

#extract and get means for each ESM for each value type

#anomalies
surfsal_anomaly1 = raster("Netcdf Files/salESM.nc",
                         level = 260, 
                         varname = "anomaly1")
surfsal_anomaly2 = raster("Netcdf Files/salESM.nc",
                         level = 260, 
                         varname = "anomaly2")
surfsal_anomaly3 = raster("Netcdf Files/salESM.nc",
                         level = 260, 
                         varname = "anomaly3")

surfsal_anomaly.df1 = raster::as.data.frame(surfsal_anomaly1, xy = TRUE)
surfsal_anomaly.df2 = raster::as.data.frame(surfsal_anomaly2, xy = TRUE)
surfsal_anomaly.df3 = raster::as.data.frame(surfsal_anomaly3, xy = TRUE)

surfsal_anomaly.df1 <-surfsal_anomaly.df1 %>%
  rename(Sea.Surface.Salinity1 = Sea.Surface.Salinity)
surfsal_anomaly.df2 <-surfsal_anomaly.df2 %>%
  rename(Sea.Surface.Salinity2 = Sea.Surface.Salinity)
surfsal_anomaly.df3 <-surfsal_anomaly.df3 %>%
  rename(Sea.Surface.Salinity3 = Sea.Surface.Salinity)

combined_df <- cbind(
  surfsal_anomaly.df1,
  surfsal_anomaly.df2["Sea.Surface.Salinity2"],
  surfsal_anomaly.df3["Sea.Surface.Salinity3"]
)

combined_df$mean.Sea.Surface.Salinity <- rowMeans(combined_df[, c("Sea.Surface.Salinity1", "Sea.Surface.Salinity2", "Sea.Surface.Salinity3")], na.rm = TRUE)

surfsal_anomaly.df <- combined_df %>%
  group_by(x, y) %>%
  summarize(Sea_Surface_Salinity = mean(mean.Sea.Surface.Salinity, na.rm = TRUE))

#historical values

surfsal_histclim1 = raster("Netcdf Files/salESM.nc",
                          level = 260, 
                          varname = "histclim1")
surfsal_histclim2 = raster("Netcdf Files/salESM.nc",
                           level = 260, 
                           varname = "histclim2")
surfsal_histclim3 = raster("Netcdf Files/salESM.nc",
                           level = 260, 
                           varname = "histclim3")

surfsal_histclim.df1 = raster::as.data.frame(surfsal_histclim1, xy = TRUE)
surfsal_histclim.df2 = raster::as.data.frame(surfsal_histclim2, xy = TRUE)
surfsal_histclim.df3 = raster::as.data.frame(surfsal_histclim3, xy = TRUE)

surfsal_histclim.df1 <-surfsal_histclim.df1 %>%
  rename(Sea.Surface.Salinity1 = Sea.Surface.Salinity)
surfsal_histclim.df2 <-surfsal_histclim.df2 %>%
  rename(Sea.Surface.Salinity2 = Sea.Surface.Salinity)
surfsal_histclim.df3 <-surfsal_histclim.df3 %>%
  rename(Sea.Surface.Salinity3 = Sea.Surface.Salinity)

combined_df <- cbind(
  surfsal_histclim.df1,
  surfsal_histclim.df2["Sea.Surface.Salinity2"],
  surfsal_histclim.df3["Sea.Surface.Salinity3"]
)

combined_df$mean.Sea.Surface.Salinity <- rowMeans(combined_df[, c("Sea.Surface.Salinity1", "Sea.Surface.Salinity2", "Sea.Surface.Salinity3")], na.rm = TRUE)

surfsal_histclim.df <- combined_df %>%
  group_by(x, y) %>%
  summarize(Sea_Surface_Salinity = mean(mean.Sea.Surface.Salinity, na.rm = TRUE))


#merge anomalies with historical values
surfsal_proj.df <- full_join(surfsal_anomaly.df, surfsal_histclim.df, by = c("x", "y"))
surfsal_proj.df<-na.omit(surfsal_proj.df)

# Merge the data frames by x and y, and then sum the anomalies to the historical values
surfsal_proj.df <- surfsal_proj.df %>%
  na.omit() %>%
  rename(Sea.Surface.Sal.Anom = Sea_Surface_Salinity.x,
         Sea.Surface.Sal.Hist = Sea_Surface_Salinity.y) %>%
  mutate(Sea.Surface.Sal.Proj = Sea.Surface.Sal.Anom + Sea.Surface.Sal.Hist)


surfsal_proj.df <- surfsal_proj.df %>%
  mutate(x = ifelse(x > 180, x - 360, x))

#calculcate bias.. which is 1.11... apply after calculcating procted salinities
surfsal_bias<-mean(surfsal_proj.df$Sea.Surface.Sal.Hist)-mean(historicalclimatedata$SALSURF, na.rm = TRUE)
print(surfsal_bias)


# Bias correct the columns in surfsal_proj.df
surfsal_proj.df <- surfsal_proj.df %>%
  mutate(Sea.Surface.Sal.Proj = Sea.Surface.Sal.Proj - surfsal_bias,
         Sea.Surface.Sal.Hist = Sea.Surface.Sal.Hist - surfsal_bias)


# Calculate distances to land for each observation in the salinity data
surfsal_proj.df<-ggOceanMaps::dist2land(surfsal_proj.df)


#set lat and lon bounds according to the observational data
lat_bounds <- range(springplankton4$lat)
lon_bounds <- range((springplankton4$lon))

# Filter the salinity data to include only observations within the bounds of zoopdfspring and values outside 100 km from shore
surfsal_proj.df_filtered <- surfsal_proj.df %>%
  filter(ldist > 100, 
         x >= min(lon_bounds) & x <= max(lon_bounds) & 
           y >= min(lat_bounds) & y <= max(lat_bounds)) 

springplankton4$SALSURF

# filter and remove values outside of gps bounds of interest 
lon_min <- -95
lon_max <- -90
lat_min <- 28
lat_max <- 29

surfsal_proj.df_filtered2 <- surfsal_proj.df_filtered %>%
  filter(!(x >= lon_min & x <= lon_max & y >= lat_min & y <= lat_max))



#### CHLOROPHYLL
 
#anomalies
chlor_anomaly1 = raster("Netcdf Files/chlorESM2.nc",
                          level = 260, 
                          varname = "anomaly1")
chlor_anomaly2 = raster("Netcdf Files/chlorESM2.nc",
                          level = 260, 
                          varname = "anomaly2")
chlor_anomaly3 = raster("Netcdf Files/chlorESM2.nc",
                          level = 260, 
                          varname = "anomaly3")

chlor_anomaly.df1 = raster::as.data.frame(chlor_anomaly1, xy = TRUE)
chlor_anomaly.df2 = raster::as.data.frame(chlor_anomaly2, xy = TRUE)
chlor_anomaly.df3 = raster::as.data.frame(chlor_anomaly3, xy = TRUE)

chlor_anomaly.df1 <-chlor_anomaly.df1 %>%
  rename(Chlorophyll1 =Mass.Concentration.of.Total.Phytoplankton.Expressed.as.Chlorophyll.in.Sea.Water)
chlor_anomaly.df2 <-chlor_anomaly.df2 %>%
  rename(Chlorophyll2 =Mass.Concentration.of.Total.Phytoplankton.Expressed.as.Chlorophyll.in.Sea.Water)
chlor_anomaly.df3 <-chlor_anomaly.df3 %>%
  rename(Chlorophyll3 =Mass.Concentration.of.Total.Phytoplankton.Expressed.as.Chlorophyll.in.Sea.Water)

combined_df <- cbind(
  chlor_anomaly.df1,
  chlor_anomaly.df2["Chlorophyll2"],
  chlor_anomaly.df3["Chlorophyll3"])

combined_df$mean.Chlorophyll <- rowMeans(combined_df[, c("Chlorophyll1", "Chlorophyll2", "Chlorophyll3")], na.rm = TRUE)

chlor_anomaly.df <- combined_df %>%
  group_by(x, y) %>%
  summarize(Chlorophyll = mean(mean.Chlorophyll, na.rm = TRUE))

chlor_anomaly.df<- chlor_anomaly.df %>%
  mutate(PP.Anom = Chlorophyll / 1000)

# projections
chlor_proj1 = raster("Netcdf Files/chlorESM2.nc",
                        level = 260, 
                        varname = "histclim1")
chlor_proj2 = raster("Netcdf Files/chlorESM2.nc",
                        level = 260, 
                        varname = "histclim2")
chlor_proj3 = raster("Netcdf Files/chlorESM2.nc",
                        level = 260, 
                        varname = "histclim3")

chlor_proj.df1 = raster::as.data.frame(chlor_proj1, xy = TRUE)
chlor_proj.df2 = raster::as.data.frame(chlor_proj2, xy = TRUE)
chlor_proj.df3 = raster::as.data.frame(chlor_proj3, xy = TRUE)

chlor_proj.df1 <-chlor_proj.df1 %>%
  rename(Chlorophyll1 = Mass.Concentration.of.Total.Phytoplankton.Expressed.as.Chlorophyll.in.Sea.Water)
chlor_proj.df2 <-chlor_proj.df2 %>%
  rename(Chlorophyll2 = Mass.Concentration.of.Total.Phytoplankton.Expressed.as.Chlorophyll.in.Sea.Water)
chlor_proj.df3 <-chlor_proj.df3 %>%
  rename(Chlorophyll3 = Mass.Concentration.of.Total.Phytoplankton.Expressed.as.Chlorophyll.in.Sea.Water)

combined_df <- cbind(
  chlor_proj.df1,
  chlor_proj.df2["Chlorophyll2"],
  chlor_proj.df3["Chlorophyll3"])

combined_df$mean.Chlorophyll <- rowMeans(combined_df[, c("Chlorophyll1", "Chlorophyll2", "Chlorophyll3")], na.rm = TRUE)

chlor_hist.df <- combined_df %>%
  group_by(x, y) %>%
  summarize(Chlorophyll = mean(mean.Chlorophyll, na.rm = TRUE))

chlor_hist.df<- chlor_hist.df %>%
  mutate(PP.Hist = Chlorophyll / 1000)

#merge anomalies with historical values
chlor_proj.df <- full_join(chlor_anomaly.df, chlor_hist.df, by = c("x", "y"))
chlor_proj.df<-na.omit(chlor_proj.df)

# Merge the data frames by x and y, and then sum the anomalies to the historical values
chlor_proj.df <- chlor_proj.df %>%
  na.omit() %>%
  mutate(PP.Proj =PP.Anom + PP.Hist)


chlor_proj.df <- chlor_proj.df %>%
  mutate(x = ifelse(x > 180, x - 360, x))

# calculcate bias.. 
chlor_bias<-mean(chlor_proj.df$PP.Hist)-mean(historicalclimatedata$CHLORSURF, na.rm = TRUE)
print(chlor_bias)


# Bias correct the columns in surfsal_proj.df
chlor_proj.df <- chlor_proj.df %>%
  mutate(PP.Proj = PP.Proj - chlor_bias,
         PP.Hist = PP.Hist - chlor_bias)


# Calculate distances to land for each observation in the salinity data
chlor_proj.df<-ggOceanMaps::dist2land(chlor_proj.df)

# set lat and lon bounds according to the observational data
lat_bounds <- range(springplankton4$lat)
lon_bounds <- range((springplankton4$lon))

# Filter the salinity data to include only observations within the bounds of zoopdfspring and values outside 100 km from shore
chlor_proj.df_filtered <- chlor_proj.df %>%
  filter(ldist > 100, 
         x >= min(lon_bounds) & x <= max(lon_bounds) & 
           y >= min(lat_bounds) & y <= max(lat_bounds)) 

# filter and remove values outside of gps bounds of interest 
lon_min <- -95
lon_max <- -90
lat_min <- 28
lat_max <- 29

chlor_proj.df_filtered2 <- chlor_proj.df_filtered %>%
  filter(!(x >= lon_min & x <= lon_max & y >= lat_min & y <= lat_max))


 #make custom theme for maps
theme_ray_map <- function() {
  theme(
    legend.position = c(.45, -.48), 
    legend.direction = "horizontal", 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    legend.background = element_rect(fill = "transparent"), 
    text = element_text(size = 7),
    legend.text = element_text(color = "black", size = 7),
    plot.margin=unit(c(-.5,.5,-.1,.1),"cm"),
    panel.background = element_rect(fill = "grey90",color = "black"),
    panel.grid = element_line(color = "grey90", size = 0.8),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5, linetype = "solid")
  )}

PPprojplot<-ggplot() + 
  geom_raster(data = chlor_proj.df_filtered2, aes(x = x, y = y, fill = PP.Proj), interpolate = FALSE,na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 24)) +
  theme_ray_map() + 
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    breaks = pretty_breaks(n = 5),  
    guide = guide_colorbar(
      title = "Chlorophyll (μg/L) Future Period (2070-2099)",
      title.position = "top",
      barwidth = 8,      
      barheight = .2           
    )
  ) 
PPprojplot

PPhistplot<-ggplot() + 
  geom_raster(data = chlor_proj.df_filtered2, aes(x = x, y = y, fill = PP.Hist), interpolate = FALSE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 24)) +
  theme_ray_map() + 
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    breaks = pretty_breaks(n = 5),  
    guide = guide_colorbar(
      title = "Chlorophyll (μg/L) Historical Period (1985-2014)",
      title.position = "top",
      barwidth = 8,      
      barheight = .2           
    )
  ) 
PPhistplot

PPanomplot<-ggplot() + 
  geom_raster(data = chlor_proj.df_filtered2, aes(x = x, y = y, fill = PP.Anom), interpolate = FALSE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 24)) +
  theme_ray_map() + 
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    breaks = pretty_breaks(n = 5),  
    guide = guide_colorbar(
      title = "Chlorophyll Anomaly (μg/L)",
      title.position = "top",
      barwidth = 8,      
      barheight = .2           
    )
  ) 
PPanomplot


# Plot the bias corrected, filtered projected data
salprojplot<-ggplot() + 
  geom_raster(data = surfsal_proj.df_filtered2, aes(x = x, y = y, fill = Sea.Surface.Sal.Proj), interpolate = FALSE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 24)) +
  theme_ray_map() + 
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    breaks = pretty_breaks(n = 5),  
    guide = guide_colorbar(
      title = "Salinity Future Period (2070-2099)",
      title.position = "top",
      barwidth = 8,      
      barheight = .2           
    ))
salprojplot

# Plot the bias corrected, filtered projected data
salhistplot<-ggplot() + 
  geom_raster(data = surfsal_proj.df_filtered2, aes(x = x, y = y, fill = Sea.Surface.Sal.Hist), interpolate = FALSE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 24)) +
  theme_ray_map() + 
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    breaks = pretty_breaks(n = 5),  
    guide = guide_colorbar(
      title = "Salinity Historical Period (1985-2014)",
      title.position = "top",
      barwidth = 8,      
      barheight = .2           
    )
  ) 
salhistplot

salanomplot<-ggplot() + 
  geom_raster(data = surfsal_proj.df_filtered2, aes(x = x, y = y, fill = Sea.Surface.Sal.Anom), interpolate = FALSE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 24)) +
  theme_ray_map() + 
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    breaks = pretty_breaks(n = 5),  
    guide = guide_colorbar(
      title = "Salinity Anomaly",
      title.position = "top",
      barwidth = 8,      
      barheight = .2           
    )
  ) 
salanomplot

salandchlor<-grid.arrange(arrangeGrob(salhistplot, salprojplot, salanomplot, 
                                    PPhistplot, PPprojplot, PPanomplot,
                                    ncol=3, nrow=2),heights=c(.7,.5))


ggsave("salandchlor.png", salandchlor, dpi = 250, bg = "white",
       width = 2000,
       height = 1400,
       units = "px") 

#projections
sst_proj1 = raster("Netcdf Files/sstESM.nc",
                     level = 260, 
                     varname = "histclim1")
sst_proj2 = raster("Netcdf Files/sstESM.nc",
                     level = 260, 
                     varname = "histclim2")
sst_proj3 = raster("Netcdf Files/sstESM.nc",
                     level = 260, 
                     varname = "histclim3")
sst_proj4 = raster("Netcdf Files/sstESM.nc",
                   level = 260, 
                   varname = "histclim4")


sst_proj.df1 = raster::as.data.frame(sst_proj1, xy = TRUE)
sst_proj.df2 = raster::as.data.frame(sst_proj2, xy = TRUE)
sst_proj.df3 = raster::as.data.frame(sst_proj3, xy = TRUE)
sst_proj.df4 = raster::as.data.frame(sst_proj4, xy = TRUE)

sst_proj.df1 <-sst_proj.df1 %>%
  rename(sst1 =Sea.Surface.Temperature)
sst_proj.df2 <-sst_proj.df2 %>%
  rename(sst2= Sea.Surface.Temperature)
sst_proj.df3 <-sst_proj.df3 %>%
  rename(sst3 = Sea.Surface.Temperature)
sst_proj.df4 <-sst_proj.df4 %>%
  rename(sst4 = Sea.Surface.Temperature)

combined_df <- cbind(
  sst_proj.df1,
  sst_proj.df2["sst2"],
  sst_proj.df3["sst3"],
  sst_proj.df4["sst4"])


combined_df$mean.sst <- rowMeans(combined_df[, c("sst1", "sst2", "sst3")], na.rm = TRUE)

sst_hist.df <- combined_df %>%
  group_by(x, y) %>%
  summarize(sst = mean(mean.sst, na.rm = TRUE))


#anomalies 
sst_anom1 = raster("Netcdf Files/sstESM.nc",
                   level = 260, 
                   varname = "anomaly1")
sst_anom2 = raster("Netcdf Files/sstESM.nc",
                   level = 260, 
                   varname = "anomaly2")
sst_anom3 = raster("Netcdf Files/sstESM.nc",
                   level = 260, 
                   varname = "anomaly3")
sst_anom4 = raster("Netcdf Files/sstESM.nc",
                   level = 260, 
                   varname = "anomaly4")


sst_anom.df1 = raster::as.data.frame(sst_anom1, xy = TRUE)
sst_anom.df2 = raster::as.data.frame(sst_anom2, xy = TRUE)
sst_anom.df3 = raster::as.data.frame(sst_anom3, xy = TRUE)
sst_anom.df4 = raster::as.data.frame(sst_anom4, xy = TRUE)

sst_anom.df1 <-sst_anom.df1 %>%
  rename(sst1 =Sea.Surface.Temperature)
sst_anom.df2 <-sst_anom.df2 %>%
  rename(sst2= Sea.Surface.Temperature)
sst_anom.df3 <-sst_anom.df3 %>%
  rename(sst3 = Sea.Surface.Temperature)
sst_anom.df4 <-sst_anom.df4 %>%
  rename(sst4 = Sea.Surface.Temperature)

combined_df <- cbind(
  sst_anom.df1,
  sst_anom.df2["sst2"],
  sst_anom.df3["sst3"],
  sst_anom.df4["sst4"])

combined_df$mean.sst <- rowMeans(combined_df[, c("sst1", "sst2", "sst3")], na.rm = TRUE)

sst_anom.df <- combined_df %>%
  group_by(x, y) %>%
  summarize(sst = mean(mean.sst, na.rm = TRUE))



# merge anomalies with historical values
sst_proj.df <- full_join(sst_anom.df, sst_hist.df, by = c("x", "y"))
sst_proj.df<-na.omit(sst_proj.df)

sst_proj.df <-sst_proj.df %>%
  rename(sst.Anom = sst.x)  %>%
  rename(sst.Hist = sst.y)


# Merge the data frames by x and y, and then sum the anomalies to the historical values
sst_proj.df <- sst_proj.df %>%
  na.omit() %>%
  mutate(sst.Proj =sst.Anom + sst.Hist)

sst_proj.df <- sst_proj.df %>%
  mutate(x = ifelse(x > 180, x - 360, x))


# Calculate distances to land for each observation in the salinity data
sst_proj.df<-ggOceanMaps::dist2land(sst_proj.df)

# Filter the salinity data to include only observations within the bounds of zoopdfspring and values outside 100 km from shore
sst_proj.df_filtered <- sst_proj.df %>%
  filter(ldist > 100, 
         x >= min(lon_bounds) & x <= max(lon_bounds) & 
           y >= min(lat_bounds) & y <= max(lat_bounds)) 

# filter and remove values outside of gps bounds of interest 
sst_proj.df_filtered2 <- sst_proj.df_filtered %>%
  filter(!(x >= lon_min & x <= lon_max & y >= lat_min & y <= lat_max))

# calculcate bias... apply after calculcating procted salinities
sst_bias<-mean(sst_proj.df_filtered2$sst.Hist)-mean(historicalclimatedata$TEMPSURF,na.rm = TRUE)
print(sst_bias)

# Bias correct the columns in surfsal_proj.df
sst_proj.df_filtered2 <- sst_proj.df_filtered2%>%
  mutate(sst.Proj = sst.Proj - sst_bias,
         sst.Hist = sst.Hist - sst_bias)


#botton temp

#anomalies 
bt_anom1 = raster("Netcdf Files/btESM.nc",
                   level = 260, 
                   varname = "anomaly1")
bt_anom2 = raster("Netcdf Files/btESM.nc",
                   level = 260, 
                   varname = "anomaly2")
bt_anom3 = raster("Netcdf Files/btESM.nc",
                   level = 260, 
                   varname = "anomaly3")

bt_anom.df1 = raster::as.data.frame(bt_anom1, xy = TRUE)
bt_anom.df2 = raster::as.data.frame(bt_anom2, xy = TRUE)
bt_anom.df3 = raster::as.data.frame(bt_anom3, xy = TRUE)

bt_anom.df1 <-bt_anom.df1 %>%
  rename(bt1 =Sea.Water.Temperature.at.200m)
bt_anom.df2 <-bt_anom.df2 %>%
  rename(bt2= Sea.Water.Temperature.at.200m)
bt_anom.df3 <-bt_anom.df3 %>%
  rename(bt3 = Sea.Water.Temperature.at.200m)

combined_df <- cbind(
  bt_anom.df1,
  bt_anom.df2["bt2"],
  bt_anom.df3["bt3"])

combined_df$mean.bt <- rowMeans(combined_df[, c("bt1", "bt2", "bt3")], na.rm = TRUE)

bt_anom.df <- combined_df %>%
  group_by(x, y) %>%
  summarize(bt.Anom = mean(mean.bt, na.rm = TRUE))

#historical and projections
bt_proj1 = raster("Netcdf Files/btESM.nc",
                  level = 260, 
                  varname = "histclim1")
bt_proj2 = raster("Netcdf Files/btESM.nc",
                  level = 260, 
                  varname = "histclim2")
bt_proj3 = raster("Netcdf Files/btESM.nc",
                  level = 260, 
                  varname = "histclim3")

bt_proj.df1 = raster::as.data.frame(bt_proj1, xy = TRUE)
bt_proj.df2 = raster::as.data.frame(bt_proj2, xy = TRUE)
bt_proj.df3 = raster::as.data.frame(bt_proj3, xy = TRUE)

bt_proj.df1 <-bt_proj.df1 %>%
  rename(bt1 =Sea.Water.Temperature.at.200m)
bt_proj.df2 <-bt_proj.df2 %>%
  rename(bt2= Sea.Water.Temperature.at.200m)
bt_proj.df3 <-bt_proj.df3 %>%
  rename(bt3 = Sea.Water.Temperature.at.200m)

combined_df <- cbind(
  bt_proj.df1,
  bt_proj.df2["bt2"],
  bt_proj.df3["bt3"])

combined_df$mean.bt <- rowMeans(combined_df[, c("bt1", "bt2", "bt3")], na.rm = TRUE)

bt_proj.df <- combined_df %>%
  group_by(x, y) %>%
  summarize(bt.Hist = mean(mean.bt, na.rm = TRUE))

# merge anomalies with historical values
bt_proj.df <- full_join(bt_anom.df, bt_proj.df, by = c("x", "y"))
bt_proj.df<-na.omit(bt_proj.df)

# Merge the data frames by x and y, and then sum the anomalies to the historical values
bt_proj.df <- bt_proj.df %>%
  na.omit() %>%
  mutate(bt.Proj =bt.Anom + bt.Hist)

bt_proj.df <- bt_proj.df %>%
  mutate(x = ifelse(x > 180, x - 360, x))

# Calculate distances to land for each observation in the salinity data
bt_proj.df<-ggOceanMaps::dist2land(bt_proj.df)

bt_proj.df_filtered <- bt_proj.df %>%
  filter(ldist > 100, 
         x >= min(lon_bounds) & x <= max(lon_bounds) & 
           y >= min(lat_bounds) & y <= max(lat_bounds)) 

# filter and remove values outside of gps bounds of interest 
bt_proj.df_filtered2 <- bt_proj.df_filtered %>%
  filter(!(x >= lon_min & x <= lon_max & y >= lat_min & y <= lat_max))

# calculcate bias.. .. apply after calculcating procted salinities
bt_bias<-mean(bt_proj.df_filtered2$bt.Hist)-mean(historicalclimatedata$TEMPMAX, na.rm = TRUE)
print(bt_bias)

# Bias correct the columns in surfsal_proj.df
bt_proj.df_filtered2<- bt_proj.df_filtered2 %>%
  mutate(bt.Proj = bt.Proj - bt_bias,
         bt.Hist = bt.Hist - bt_bias)


#bottom temp plota
btprojplot<-ggplot() + 
  geom_raster(data = bt_proj.df_filtered2, aes(x = x, y = y, fill = bt.Proj), interpolate = FALSE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 24)) +
  theme_ray_map() + 
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    breaks = pretty_breaks(n = 5),  
    guide = guide_colorbar(
      title = "Temp at Depth (°C) Future Period (2070-2099)",
      title.position = "top",
      barwidth = 8,      
      barheight = .2           
    ))
btprojplot

bthistplot<-ggplot() + 
  geom_raster(data = bt_proj.df_filtered2, aes(x = x, y = y, fill = bt.Hist), interpolate = FALSE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 24)) +
  theme_ray_map() + 
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    breaks = pretty_breaks(n = 5),  
    guide = guide_colorbar(
      title = "Temp at Depth (°C) Historical Period (1985-2014)",
      title.position = "top",
      barwidth = 8,      
      barheight = .2           
    ))
bthistplot

btanomplot<-ggplot() + 
  geom_raster(data = bt_proj.df_filtered2, aes(x = x, y = y, fill = bt.Anom), interpolate = FALSE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 24)) +
  theme_ray_map() + 
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    breaks = pretty_breaks(n = 5),  
    guide = guide_colorbar(
      title = "Temp at Depth (°C) Anomaly",
      title.position = "top",
      barwidth = 8,      
      barheight = .2           
    ))
btanomplot


#sst plots
sstprojplot<-ggplot() + 
  geom_raster(data = sst_proj.df_filtered2, aes(x = x, y = y, fill = sst.Proj), interpolate = FALSE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 24)) +
  theme_ray_map() + 
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    breaks = pretty_breaks(n = 5),  
    guide = guide_colorbar(
      title = "SST (°C) Future Period (2070-2099)",
      title.position = "top",
      barwidth = 8,      
      barheight = .2           
    ))
sstprojplot

ssthistplot<-ggplot() + 
  geom_raster(data = sst_proj.df_filtered2, aes(x = x, y = y, fill = sst.Hist), interpolate = FALSE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 24)) +
  theme_ray_map() + 
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    breaks = pretty_breaks(n = 5),  
    guide = guide_colorbar(
      title = "SST (°C) Historical Period (1985-2014)",
      title.position = "top",
      barwidth = 8,      
      barheight = .2           
    ))
ssthistplot

sstanomplot<-ggplot() + 
  geom_raster(data = sst_proj.df_filtered2, aes(x = x, y = y, fill = sst.Anom), interpolate = FALSE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 24)) +
  theme_ray_map() + 
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    breaks = pretty_breaks(n = 5),  
    guide = guide_colorbar(
      title = "SST (°C) Anomaly",
      title.position = "top",
      barwidth = 8,      
      barheight = .2           
    ))
sstanomplot

tempplots<-grid.arrange(arrangeGrob(ssthistplot, sstprojplot, sstanomplot, 
                                      bthistplot, btprojplot, btanomplot,
                                      ncol=3, nrow=2),heights=c(.7,.5))


ggsave("tempplots.png", tempplots, dpi = 250, bg = "white",
       width = 2000,
       height = 1400,
       units = "px") 
