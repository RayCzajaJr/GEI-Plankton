

### make distribution maps for each taxon, for each model type, for each time period
### uses output/values from ESM Projections and MAPs file


library(geosphere)
library(mgcv)
library(gridExtra)
library(grid)
library(dplyr)
library(ggplot2)

#get projections into one df 
dfs <- list(sst_proj.df_filtered2,  bt_proj.df_filtered2, surfsal_proj.df_filtered2, chlor_proj.df_filtered2)

#merge DataFrames for model based predictions
projections.df <- reduce(dfs, full_join, by = c("x", "y"))

#create a new df for forecast predictions
newdataproj_full <- data.frame(
  TEMPSURF = projections.df$sst.Proj,
  TEMPMAX = projections.df$bt.Proj,
  CHLORSURF = projections.df$PP.Proj,
  SALSURF = projections.df$Sea.Surface.Sal.Proj,
  lon = projections.df$x,
  lat = projections.df$y
  )

newdataproj_full <- newdataproj_full %>%
  mutate(TIME_EMIL = 142.9091)

#now historical
newdatahist_full <- data.frame(
  TEMPSURF = projections.df$sst.Hist,
  TEMPMAX = projections.df$bt.Hist,
  CHLORSURF = projections.df$PP.Hist,
  SALSURF = projections.df$Sea.Surface.Sal.Hist,
  lon = projections.df$x,
  lat = projections.df$y)

newdatahist_full <- newdatahist_full %>%
  mutate(TIME_EMIL = 142.9091)


##### ZOOPLANKTON 

#make projections for calanoids
CalGam <- gam(Calanoids ~ s(TEMPSURF, k=3, m=1) + s(SALSURF, k=3) + s(TEMPMAX, k=3, m=1) +s(CHLORSURF, k=3, m=1)+ s(lon, lat) +s(TIME_EMIL, bs="cs"),
               data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")

datahist_cala <- newdatahist_full[, !colnames(newdatahist_full) %in% "distance"]
predicted_abundances_hist_cala <- predict(CalGam, newdata = datahist_cala, type = "response")

projections.df$calanoid_projections_hist <- predicted_abundances_hist_cala
mean(projections.df$calanoid_projections_hist)

dataproj_cala <- newdataproj_full[, !colnames(newdataproj_full) %in% "distance"]
predicted_abundances_proj_cala <- predict(CalGam, newdata = dataproj_cala , type = "response")

projections.df$calanoid_projections <- predicted_abundances_proj_cala
mean(projections.df$calanoid_projections)

projections.df<- projections.df %>%
  mutate(calanoid_anom = calanoid_projections -calanoid_projections_hist)

#make theme for sdms 
theme_ray_SDM_zoop <- theme(
  legend.position = "right",
  legend.direction = "vertical",
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.background = element_rect(fill = "transparent"),
  legend.text = element_text(size = 5),
  legend.title = element_text(size = 5),
  text = element_text(size = 5),
  legend.key.height = unit(1, "lines"),
  panel.background = element_rect(fill = "grey90",color = "black"),
  panel.grid = element_line(color = "grey90", size = 0.8),
  panel.border = element_rect(color = "black", fill = NA, size = 0.5, linetype = "solid"))

common_mincala <- min(c(projections.df$calanoid_projections, projections.df$calanoid_projections_hist), na.rm = TRUE)
common_maxcala <- max(c(projections.df$calanoid_projections, projections.df$calanoid_projections_hist), na.rm = TRUE)
common_limitscala <- c(common_mincala, common_maxcala)

projcalaSDM<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = calanoid_projections), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitscala, oob = scales::squish) +
  theme_ray_SDM_zoop+
  guides(fill = guide_colorbar(barheight = unit(1.5, "cm"), barwidth = unit(0.4, "cm")))+
  labs(fill = "Calanoid\nAbundance")
projcalaSDM

projcalaSDM <- ggdraw(projcalaSDM) +
  draw_label("Future Period (2070-2099)", x = 0.25, y = 0.86, size = 4.5, color = "aliceblue")

histcalaSDM <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = calanoid_projections_hist), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitscala, oob = scales::squish) +
  theme_ray_SDM_zoop+
  guides(fill = guide_colorbar(barheight = unit(1.5, "cm"), barwidth = unit(0.4, "cm")))+
  labs(fill = "Calanoid\nAbundance") 
histcalaSDM

histcalaSDM <- ggdraw(histcalaSDM) +
  draw_label("Historical Period (1985-2014)", x = 0.25, y = 0.86, size = 4.5, color = "aliceblue")

calaSDMs<-grid.arrange(histcalaSDM, projcalaSDM,nrow=1)

center_hist_cala <- c(
  weighted.mean(projections.df$x, projections.df$calanoid_projections_hist, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$calanoid_projections_hist, na.rm = TRUE)
)
center_hist_cala

center_proj_cala <- c(
  weighted.mean(projections.df$x, projections.df$calanoid_projections, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$calanoid_projections, na.rm = TRUE)
)
center_proj_cala 

distance_meters_cala <- distHaversine(center_hist_cala, center_proj_cala )
print(distance_km_cala <- distance_meters_cala / 1000)

print(initial_bearing_cala <- bearing(center_hist_cala, center_proj_cala))


projections.df$calanoid_percentchange<-
  ((projections.df$calanoid_projections-projections.df$calanoid_projections_hist)/projections.df$calanoid_projections_hist)*100
projections.df$calanoid_percentchange

##### Caldocerans
CladGam <- gam(Cladocerans~ s(TEMPSURF, k=3) + s(SALSURF, k=3)+ s(lon, lat),
               data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")

datahist_clad <- newdatahist_full[, !colnames(newdatahist_full) %in% c("CHLORSURF","TEMPMAX", "TIME_EMIL")]
predicted_abundances_hist_clad <- predict(CladGam, newdata = datahist_clad, type = "response")

projections.df$cladoceran_projections_hist <- predicted_abundances_hist_clad
mean(projections.df$cladoceran_projections_hist)

dataproj_clad <- newdataproj_full[, !colnames(newdataproj_full) %in% c("CHLORSURF","TEMPMAX", "TIME_EMIL")]
predicted_abundances_proj_clad <- predict(CladGam, newdata = dataproj_clad , type = "response")

projections.df$cladoceran_projections <- predicted_abundances_proj_clad
mean(projections.df$cladoceran_projections)

projections.df<- projections.df %>%
  mutate(cladoceran_anom = cladoceran_projections -cladoceran_projections_hist)

common_minclad <- min(c(projections.df$cladoceran_projections, projections.df$cladoceran_projections_hist), na.rm = TRUE)
common_maxclad <- max(c(projections.df$cladoceran_projections, projections.df$cladoceran_projections_hist), na.rm = TRUE)
common_limitsclad <- c(common_minclad, common_maxclad)

projcladSDM<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = cladoceran_projections), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"),limits = c(0, 200),oob = scales::squish) +  
  theme_ray_SDM_zoop+
  guides(fill = guide_colorbar(barheight = unit(1.5, "cm"), barwidth = unit(0.4, "cm")))+
  labs(fill = "Cladoceran\nAbundance*")
projcladSDM

projcladSDM <- ggdraw(projcladSDM) +
  draw_label("Future Period (2070-2099)", x = 0.25, y = 0.86, size = 4.5, color = "aliceblue")

histcladSDM <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = cladoceran_projections_hist), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu")) +
  guides(fill = guide_colorbar(barheight = unit(1.5, "cm"), barwidth = unit(0.4, "cm")))+
  theme_ray_SDM_zoop+
  labs(fill = "Cladoceran\nAbundance*")
histcladSDM

histcladSDM <- ggdraw(histcladSDM) +
  draw_label("Historical Period (1985-2014)", x = 0.25, y = 0.86, size = 4.5, color = "aliceblue")

cladSDMs<-grid.arrange(histcladSDM, projcladSDM,nrow=1)

center_hist_clad <- c(
  weighted.mean(projections.df$x, projections.df$cladoceran_projections_hist, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$cladoceran_projections_hist, na.rm = TRUE)
)
center_hist_clad

center_proj_clad <- c(
  weighted.mean(projections.df$x, projections.df$cladoceran_projections, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$cladoceran_projections, na.rm = TRUE)
)
center_proj_clad

distance_meters_clad <- distHaversine(center_hist_clad, center_proj_clad)
print(distance_km_clad <- distance_meters_clad / 1000)
print(initial_bearing_clad <- bearing(center_hist_clad, center_proj_clad))

projections.df$cladoceran_percentchange<-
  +   ((projections.df$cladoceran_projections-projections.df$cladoceran_projections_hist)/projections.df$cladoceran_projections_hist)*100

projections.df$cladoceran_percentchange

#larvaceans
LarvGam<- gam(Larvaceans~ s(TEMPSURF, k=3, m=1) + s(TEMPMAX, k=3, m=1) + s(lon, lat),
               data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")

#make projections for larvaceans
datahist_larv <- newdatahist_full[, !colnames(newdatahist_full) %in% c("CHLORSURF", "SALSURF", "TIME_EMIL")]
predicted_abundances_hist_larv <- predict(LarvGam, newdata = datahist_larv, type = "response")

projections.df$larvacean_projections_hist <- predicted_abundances_hist_larv
mean(projections.df$larvacean_projections_hist)

dataproj_larv <- newdataproj_full[, !colnames(newdataproj_full)%in% c("CHLORSURF", "SALSURF", "TIME_EMIL")]
predicted_abundances_proj_larv <- predict(LarvGam, newdata = dataproj_larv , type = "response")

projections.df$larvacean_projections <- predicted_abundances_proj_larv
mean(projections.df$larvacean_projections)

projections.df<- projections.df %>%
  mutate(larvacean_anom = larvacean_projections -larvacean_projections_hist)

common_minlarv <- min(c(projections.df$larvacean_projections, projections.df$larvacean_projections_hist), na.rm = TRUE)
common_maxlarv <- max(c(projections.df$larvacean_projections, projections.df$larvacean_projections_hist), na.rm = TRUE)
common_limitslarv <- c(common_minlarv, common_maxlarv)

projlarvSDM<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = larvacean_projections), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitslarv, oob = scales::squish) +
  theme_ray_SDM_zoop+  
  guides(fill = guide_colorbar(barheight = unit(1.5, "cm"), barwidth = unit(0.4, "cm")))+
  labs(fill = "Larvacean\nAbundance")
projlarvSDM

projlarvSDM <- ggdraw(projlarvSDM) +
  draw_label("Future Period (2070-2099)", x = 0.25, y = 0.86, size = 4.5, color = "aliceblue")

histlarvSDM <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = larvacean_projections_hist), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitslarv, oob = scales::squish) +
  theme_ray_SDM_zoop+
  guides(fill = guide_colorbar(barheight = unit(1.5, "cm"), barwidth = unit(0.4, "cm")))+
  labs(fill = "Larvacean\nAbundance")
histlarvSDM

histlarvSDM <- ggdraw(histlarvSDM) +
  draw_label("Historical Period (1985-2014)", x = 0.25, y = 0.86, size = 4.5, color = "aliceblue")

larvSDMs<-grid.arrange(histlarvSDM, projlarvSDM,nrow=1)

center_hist_larv <- c(
  weighted.mean(projections.df$x, projections.df$larvacean_projections_hist, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$larvacean_projections_hist, na.rm = TRUE)
)
center_hist_larv

center_proj_larv <- c(
  weighted.mean(projections.df$x, projections.df$larvacean_projections, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$larvacean_projections, na.rm = TRUE)
)
center_proj_larv

distance_meters_larv <- distHaversine(center_hist_larv, center_proj_larv)
print(distance_km <- distance_meters_larv / 1000)
print(initial_bearing <- bearing(center_hist_larv, center_proj_larv))


projections.df$larvacean_percentchange<-
  +   ((projections.df$larvacean_projections-projections.df$larvacean_projections_hist)/projections.df$larvacean_projections_hist)*100

projections.df$larvacean_percentchange


#cyclopoids
CycGam <- gam(Cyclopoids ~ s(TEMPSURF, k=3, m=1) + s(SALSURF, k=3) + s(TEMPMAX, k=3, m=1) +s(CHLORSURF, k=3)+ s(lon, lat),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")

datahist_cyc <- newdatahist_full[, !colnames(newdatahist_full) %in% c("TIME_EMIL")]
predicted_abundances_hist_cyc <- predict(CycGam, newdata = datahist_cyc, type = "response")

projections.df$cyclopoid_projections_hist <- predicted_abundances_hist_cyc
mean(projections.df$cyclopoid_projections_hist)

dataproj_cyc <- newdataproj_full[, !colnames(newdataproj_full)%in% c("TIME_EMIL")]
predicted_abundances_proj_cyc <- predict(CycGam, newdata = dataproj_cyc , type = "response")

projections.df$cyclopoid_projections <- predicted_abundances_proj_cyc
mean(projections.df$cyclopoid_projections)

projections.df<- projections.df %>%
  mutate(cyclopoid_anom = cyclopoid_projections -cyclopoid_projections_hist)

common_mincyc <- min(c(projections.df$cyclopoid_projections, projections.df$cyclopoid_projections_hist), na.rm = TRUE)
common_maxcyc <- max(c(projections.df$cyclopoid_projections, projections.df$cyclopoid_projections_hist), na.rm = TRUE)
common_limitscyc<- c(common_mincyc, common_maxcyc)

projcycSDM<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = cyclopoid_projections), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitscyc, oob = scales::squish) +
  theme_ray_SDM_zoop+
  guides(fill = guide_colorbar(barheight = unit(1.5, "cm"), barwidth = unit(0.5, "cm")))+
  labs(fill = "Cyclopoid\nAbundance")
projcycSDM

projcycSDM <- ggdraw(projcycSDM) +
  draw_label("Future Period (2070-2099)", x = 0.25, y = 0.86, size = 4.5, color = "aliceblue")

histcycSDM <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = cyclopoid_projections_hist), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitscyc, oob = scales::squish) +
  theme_ray_SDM_zoop+
  guides(fill = guide_colorbar(barheight = unit(1.5, "cm"), barwidth = unit(0.4, "cm")))+
  labs(fill = "Cyclopoid\nAbundance")
histcycSDM

histcycSDM <- ggdraw(histcycSDM) +
  draw_label("Historical Period (1985-2014)", x = 0.25, y = 0.86, size = 4.5, color = "aliceblue")

cycSDMs<-grid.arrange(histcycSDM, projcycSDM,nrow=1)

center_hist_cyc <- c(
  weighted.mean(projections.df$x, projections.df$cyclopoid_projections_hist, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$cyclopoid_projections_hist, na.rm = TRUE)
)
center_hist_cyc

center_proj_cyc <- c(
  weighted.mean(projections.df$x, projections.df$cyclopoid_projections, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$cyclopoid_projections, na.rm = TRUE)
)
center_proj_cyc

distance_meters_cyc <- distHaversine(center_hist_cyc, center_proj_cyc)
print(distance_km_cyc<- distance_meters_cyc / 1000)
print(initial_bearing <- bearing(center_hist_cyc, center_proj_cyc))


projections.df$cyclopoid_percentchange<-
  +   ((projections.df$cyclopoid_projections-projections.df$cyclopoid_projections_hist)/projections.df$cyclopoid_projections_hist)*100

projections.df$cyclopoid_percentchange

zoopSDMplots<-grid.arrange(larvSDMs, cycSDMs,calaSDMs,cladSDMs,
                       ncol = 1, nrow = 4)

ggsave("zoopSDMplots.png", zoopSDMplots, dpi = 350, bg = "white",
       width = 1961,
       height = 1400,
       units = "px") 


##### ICHTHYOPLANKTON climate

theme_ray_SDM_fish <- theme(
  legend.position = c(.5, -.48), 
  legend.direction = "horizontal", 
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.background = element_rect(fill = "transparent"),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 8),
  text = element_text(size = 8),
  legend.key.height = unit(1, "lines"),
  panel.background = element_rect(fill = "grey90",color = "black"),
  panel.grid = element_line(color = "grey90", size = 0.8),
  legend.box.margin = margin(12, 0, 0, 0) ,
  plot.margin = unit(c(0.5, 0.5, .5, 0.5), "cm"),
  panel.border = element_rect(color = "black", fill = NA, size = 0.5, linetype = "solid"))

custom_guide_sdm <- guide_colorbar(
  barwidth = unit(8, "lines"),
  barheight = unit(0.5, "lines"),
  title.position = "bottom",
  label.position = "top"
)

#### billfish climate gam
springplankton5 <- springplankton4 %>%
  filter(!(ISTIOPHORIDAE > 0 & TEMPSURF < 25))

ISTGam<- gam(ISTIOPHORIDAE~s(TEMPSURF, k=3, m=1) + s(TEMPMAX, k=3) + s(lon, lat),
              data =springplankton5, tw(),method = "REML", select = TRUE, na.action = "na.fail")

datahist_bill <- newdatahist_full[, !colnames(newdatahist_full) %in% c("TIME_EMIL","SALSURF","CHLORSURF")]
predicted_abundances_hist_bill <- predict(ISTGam, newdata = datahist_bill, type = "response")

projections.df$billfish_projections_hist <- predicted_abundances_hist_bill
mean(projections.df$billfish_projections_hist)

dataproj_bill <- newdataproj_full[, !colnames(newdataproj_full)%in% c("TIME_EMIL","SALSURF","CHLORSURF")]
predicted_abundances_proj_bill <- predict(ISTGam, newdata = dataproj_bill, type = "response")

projections.df$billfish_projections <- predicted_abundances_proj_bill
mean(projections.df$billfish_projections)

projections.df<- projections.df %>%
  mutate(billfish_anom = billfish_projections -billfish_projections_hist)

center_hist_bill <- c(
  weighted.mean(projections.df$x, projections.df$billfish_projections_hist, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$billfish_projections_hist, na.rm = TRUE)
)
center_hist_bill

center_proj_bill <- c(
  weighted.mean(projections.df$x, projections.df$billfish_projections, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$billfish_projections, na.rm = TRUE)
)
center_proj_bill

distance_meters_bill <- distHaversine(center_hist_bill, center_proj_bill)
print(distance_km_bill<- distance_meters_bill / 1000)
print(initial_bearing_bill <- bearing(center_hist_bill, center_proj_bill))

projections.df$billfish_percentchange<-
  +   ((projections.df$billfish_projections-projections.df$billfish_projections_hist)/projections.df$billfish_projections_hist)*100

projections.df$billfish_percentchange



#Thunnus spp gam
THUNGam<- gam(THUNNUS ~ s(TEMPSURF, k=3, m=1) + s(SALSURF, k=3, m=1) + s(TEMPMAX, k=3) +s(CHLORSURF, k=3)+ s(lon, lat),
                  data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(THUNGam)

datahist_thunnus <- newdatahist_full[, !colnames(newdatahist_full) %in% c("TIME_EMIL")]
predicted_abundances_hist_thunnus <- predict(THUNGam, newdata = datahist_thunnus, type = "response")

projections.df$thunnus_projections_hist <- predicted_abundances_hist_thunnus
mean(projections.df$thunnus_projections_hist)

dataproj_thunnus <- newdataproj_full[, !colnames(newdataproj_full)%in% c("TIME_EMIL")]
predicted_abundances_proj_thunnus <- predict(THUNGam, newdata = dataproj_thunnus, type = "response")

projections.df$thunnus_projections <- predicted_abundances_proj_thunnus
mean(projections.df$thunnus_projections)

projections.df<- projections.df %>%
  mutate(thunnus_anom = thunnus_projections -thunnus_projections_hist)

center_hist_thunnus <- c(
  weighted.mean(projections.df$x, projections.df$thunnus_projections_hist, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$thunnus_projections_hist, na.rm = TRUE)
)
center_hist_thunnus

center_proj_thunnus <- c(
  weighted.mean(projections.df$x, projections.df$thunnus_projections, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$thunnus_projections, na.rm = TRUE)
)
center_proj_thunnus

distance_meters_thunnus <- distHaversine(center_hist_thunnus, center_proj_thunnus)
print(distance_km_thunnus<- distance_meters_thunnus / 1000)
print(initial_bearing_thunnus <- bearing(center_hist_thunnus, center_proj_thunnus))


projections.df$thunnus_percentchange<-
  +   ((projections.df$thunnus_projections-projections.df$thunnus_projections_hist)/projections.df$thunnus_projections_hist)*100

projections.df$thunnus_percentchange



#Thunnus thynnus
THYNGam<- gam(THUNNUSTHYNNUS~  s(TEMPSURF, k=8) +s(TEMPMAX, k=3, m=1)+s(CHLORSURF,k=3) +s(lon, lat),
               data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(THYNGam)

datahist_thynnus <- newdatahist_full[, !colnames(newdatahist_full) %in% c("TIME_MIL","SALSURF")]
predicted_abundances_hist_thynnus <- predict(THYNGam, newdata = datahist_thynnus, type = "response")

projections.df$thynnus_projections_hist <- predicted_abundances_hist_thynnus
mean(projections.df$thynnus_projections_hist)

dataproj_thynnus <- newdataproj_full[, !colnames(newdataproj_full)%in% c("TIME_EMIL","SALSURF")]
predicted_abundances_proj_thynnus <- predict(THYNGam, newdata = dataproj_thynnus, type = "response")

projections.df$thynnus_projections <- predicted_abundances_proj_thynnus
mean(projections.df$thynnus_projections)

projections.df<- projections.df %>%
  mutate(thynnus_anom = thynnus_projections -thynnus_projections_hist)

center_hist_thynnus <- c(
  weighted.mean(projections.df$x, projections.df$thynnus_projections_hist, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$thynnus_projections_hist, na.rm = TRUE)
)
center_hist_thynnus

center_proj_thynnus <- c(
  weighted.mean(projections.df$x, projections.df$thynnus_projections, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$thynnus_projections, na.rm = TRUE)
)
center_proj_thynnus

distance_meters_thynnus <- distHaversine(center_hist_thynnus, center_proj_thynnus)
print(distance_km_thynnus<- distance_meters_thynnus / 1000)
print(initial_bearing_thynnus <- bearing(center_hist_thynnus, center_proj_thynnus))


projections.df$thynnus_percentchange<-
  +   ((projections.df$thynnus_projections-projections.df$thynnus_projections_hist)/projections.df$thynnus_projections_hist)*100

projections.df$thynnus_percentchange


### mahi mahi
CORYGam<- gam(CORYPHAENA~s(TEMPSURF, k=3, m=1) + s(SALSURF, k=3, m=1),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(CORYGam)

datahist_cory <- newdatahist_full[, !colnames(newdatahist_full) %in% c("CHLORSURF","TIME_EMIL", "TEMPMAX", "lon", "lat")]
predicted_abundances_hist_cory <- predict(CORYGam, newdata = datahist_cory, type = "response")

projections.df$cory_projections_hist <- predicted_abundances_hist_cory
mean(projections.df$cory_projections_hist)

dataproj_cory<- newdataproj_full[, !colnames(newdataproj_full) %in% c("CHLORSURF","TIME_EMIL", "TEMPMAX", "lon", "lat")]
predicted_abundances_proj_cory <- predict(CORYGam, newdata = dataproj_cory, type = "response")

projections.df$cory_projections <- predicted_abundances_proj_cory
mean(projections.df$cory_projections)

center_hist_cory <- c(
  weighted.mean(projections.df$x, projections.df$cory_projections_hist, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$cory_projections_hist, na.rm = TRUE)
)
center_hist_cory

center_proj_cory <- c(
  weighted.mean(projections.df$x, projections.df$cory_projections, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$cory_projections, na.rm = TRUE)
)
center_proj_cory

distance_meters_cory <- distHaversine(center_hist_cory, center_proj_cory)
print(distance_km_cory<- distance_meters_cory / 1000)
print(initial_bearing_cory <- bearing(center_hist_cory, center_proj_cory))

projections.df$cory_percentchange<-
  +   ((projections.df$cory_projections-projections.df$cory_projections_hist)/projections.df$cory_projections_hist)*100

projections.df$cory_percentchange

###skip jack tuna gam 
KATSUGam<- gam(KATSUWONUSPELAMIS~ s(CHLORSURF, k=3) + s(SALSURF, k=3, m=1) + s(TEMPSURF, k=3, m=1),
               data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(KATSUGam)

datahist_kat <- newdatahist_full[, !colnames(newdatahist_full) %in% c("TEMPMAX", "TIME_EMIL", "lat", "lon")]
predicted_abundances_hist_kat <- predict(KATSUGam, newdata = datahist_kat, type = "response")

projections.df$kat_projections_hist <- predicted_abundances_hist_kat
mean(projections.df$kat_projections_hist)

dataproj_kat<- newdataproj_full[, !colnames(newdataproj_full)%in% c("TEMPMAX", "TIME_EMIL", "lat", "lon")]
predicted_abundances_proj_kat <- predict(KATSUGam, newdata = dataproj_kat, type = "response")

projections.df$kat_projections <- predicted_abundances_proj_kat
mean(projections.df$kat_projections)

center_hist_kat <- c(
  weighted.mean(projections.df$x, projections.df$kat_projections_hist, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$kat_projections_hist, na.rm = TRUE)
)
center_hist_kat

center_proj_kat <- c(
  weighted.mean(projections.df$x, projections.df$kat_projections, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$kat_projections, na.rm = TRUE)
)
center_proj_kat

distance_meters_kat <- distHaversine(center_hist_kat, center_proj_kat)
print(distance_km_kat<- distance_meters_kat / 1000)
print(initial_bearing_kat <- bearing(center_hist_kat, center_proj_kat))

projections.df$kat_percentchange<-
  +   ((projections.df$kat_projections-projections.df$kat_projections_hist)/projections.df$kat_projections_hist)*100

projections.df$kat_percentchange

#### frigate gam
AUXGam<- gam(AUXIS~ s(TEMPSURF, k=3, m=1)+s(TEMPMAX, k=3)+
               s(lon, lat),
             data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(AUXGam)

datahist_aux <- newdatahist_full[, !colnames(newdatahist_full) %in% c("SALSURF", "TIME_EMIL", "CHLORSURF")]
predicted_abundances_hist_aux <- predict(AUXGam, newdata = datahist_aux, type = "response")

projections.df$aux_projections_hist <- predicted_abundances_hist_aux
mean(projections.df$aux_projections_hist)

dataproj_aux<- newdataproj_full[, !colnames(newdataproj_full) %in% c("SALSURF", "TIME_EMIL", "CHLORSURF")]
predicted_abundances_proj_aux <- predict(AUXGam, newdata = dataproj_aux, type = "response")

projections.df$aux_projections <- predicted_abundances_proj_aux
mean(projections.df$aux_projections)

center_hist_aux <- c(
  weighted.mean(projections.df$x, projections.df$aux_projections_hist, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$aux_projections_hist, na.rm = TRUE)
)
center_hist_aux

center_proj_aux <- c(
  weighted.mean(projections.df$x, projections.df$aux_projections, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$aux_projections, na.rm = TRUE)
)
center_proj_aux

distance_meters_aux <- distHaversine(center_hist_aux, center_proj_aux)
print(distance_km_aux<- distance_meters_aux / 1000)
print(initial_bearing_aux <- bearing(center_hist_aux, center_proj_aux))

projections.df$aux_percentchange<-
  +   ((projections.df$aux_projections-projections.df$aux_projections_hist)/projections.df$aux_projections_hist)*100

projections.df$aux_percentchange

#### tunny gam
EUTHGam <- gam(EUTHYNNUSALLETTERATUS ~s(TEMPSURF, k=3, m=1)+s(TEMPMAX, k=3,m=1)+
                 s(lon, lat),
               data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(EUTHGam)

datahist_eut <- newdatahist_full[, !colnames(newdatahist_full) %in% c("SALSURF", "CHLORSURF", "TIME_EMIL")]
predicted_abundances_hist_eut <- predict(EUTHGam, newdata = datahist_eut, type = "response")

projections.df$eut_projections_hist <- predicted_abundances_hist_eut
mean(projections.df$eut_projections_hist)

dataproj_eut<- newdataproj_full[, !colnames(newdataproj_full)%in% c("SALSURF", "CHLORSURF", "TIME_EMIL")]
predicted_abundances_proj_eut <- predict(EUTHGam, newdata = dataproj_eut, type = "response")

projections.df$eut_projections <- predicted_abundances_proj_eut
mean(projections.df$eut_projections)

center_hist_eut<- c(
  weighted.mean(projections.df$x, projections.df$eut_projections_hist, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$eut_projections_hist, na.rm = TRUE)
)
center_hist_eut

center_proj_eut <- c(
  weighted.mean(projections.df$x, projections.df$eut_projections, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$eut_projections, na.rm = TRUE)
)
center_proj_eut

distance_meters_eut <- distHaversine(center_hist_eut, center_proj_eut)
print(distance_km_eut<- distance_meters_eut / 1000)
print(initial_bearing_eut <- bearing(center_hist_eut, center_proj_eut))

projections.df$eut_percentchange<-
  +   ((projections.df$eut_projections-projections.df$eut_projections_hist)/projections.df$eut_projections_hist)*100

projections.df$eut_percentchange


##### ICHTHYOPLANKTON prey

newdataproj_fullprey <- data.frame(
  Calanoids = projections.df$calanoid_projections,
  Cyclopoids = projections.df$cyclopoid_projections,
  Cladocerans = projections.df$cladoceran_projections,
  Larvaceans = projections.df$larvacean_projections,
  lon = projections.df$x,
  lat = projections.df$y)


newdatahist_fullprey <- data.frame(
  Calanoids = projections.df$calanoid_projections_hist,
  Cyclopoids = projections.df$cyclopoid_projections_hist,
  Cladocerans = projections.df$cladoceran_projections_hist,
  Larvaceans = projections.df$larvacean_projections_hist,
  lon = projections.df$x,
  lat = projections.df$y)

#building gams for each of these 7 taxa

#### billfish prey gam
springplankton5 <- springplankton4 %>%
  filter(!(ISTIOPHORIDAE > 0 & TEMPSURF < 25))

ISTGam2<- gam(ISTIOPHORIDAE~ s(Cladocerans, k=3) + s(lon, lat),
             data =springplankton5, tw(),method = "REML", select = TRUE, na.action = "na.fail")

datahist_bill2 <- newdatahist_fullprey[, !colnames(newdatahist_fullprey) %in% c("Calanoids","Cyclopoids", "Larvaceans")]
predicted_abundances_hist_bill2 <- predict(ISTGam2, newdata = datahist_bill2, type = "response")

projections.df$billfish_projections_hist2 <- predicted_abundances_hist_bill2
mean(projections.df$billfish_projections_hist2)

dataproj_bill2 <- newdataproj_fullprey[, !colnames(newdataproj_fullprey)%in% c("Calanoids","Cyclopoids","Larvaceans")]
predicted_abundances_proj_bill2 <- predict(ISTGam2, newdata = dataproj_bill2, type = "response")

projections.df$billfish_projections2 <- predicted_abundances_proj_bill2
mean(projections.df$billfish_projections2)

projections.df<- projections.df %>%
  mutate(billfish_anom2 = billfish_projections2 -billfish_projections_hist2)

center_hist_bill2 <- c(
  weighted.mean(projections.df$x, projections.df$billfish_projections_hist2, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$billfish_projections_hist2, na.rm = TRUE)
)
center_hist_bill2

center_proj_bill2 <- c(
  weighted.mean(projections.df$x, projections.df$billfish_projections2, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$billfish_projections2, na.rm = TRUE)
)
center_proj_bill2

distance_meters_bill2 <- distHaversine(center_hist_bill2, center_proj_bill2)
print(distance_km_bill2<- distance_meters_bill2 / 1000)
print(initial_bearing_bill2 <- bearing(center_hist_bill2, center_proj_bill2))

projections.df$billfish_percentchange2<-
  +   ((projections.df$billfish_projections2-projections.df$billfish_projections_hist2)/projections.df$billfish_projections_hist2)*100

projections.df$billfish_percentchange2

### mahi mahi
CORYGam2<- gam(CORYPHAENA~s(Larvaceans, k=3,m=1) + s(lon, lat),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(CORYGam2)

datahist_cory2 <- newdatahist_fullprey[, !colnames(newdatahist_fullprey) %in% c("Calanoids","Cyclopoids", "Cladocerans")]
predicted_abundances_hist_cory2 <- predict(CORYGam2, newdata = datahist_cory2, type = "response")

projections.df$cory_projections_hist2 <- predicted_abundances_hist_cory2
mean(projections.df$cory_projections_hist2)

dataproj_cory2<- newdataproj_fullprey[, !colnames(newdataproj_fullprey) %in% c("Calanoids","Cyclopoids", "Cladocerans")]
predicted_abundances_proj_cory2 <- predict(CORYGam2, newdata = dataproj_cory2, type = "response")

projections.df$cory_projections2 <- predicted_abundances_proj_cory2
mean(projections.df$cory_projections2)

center_hist_cory2 <- c(
  weighted.mean(projections.df$x, projections.df$cory_projections_hist2, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$cory_projections_hist2, na.rm = TRUE)
)
center_hist_cory2

center_proj_cory2 <- c(
  weighted.mean(projections.df$x, projections.df$cory_projections2, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$cory_projections2, na.rm = TRUE)
)
center_proj_cory2

distance_meters_cory2 <- distHaversine(center_hist_cory2, center_proj_cory2)
print(distance_km_cory2<- distance_meters_cory2 / 1000)
print(initial_bearing_cory2 <- bearing(center_hist_cory2, center_proj_cory2))

projections.df$cory_percentchange2<-
  +   ((projections.df$cory_projections2-projections.df$cory_projections_hist2)/projections.df$cory_projections_hist2)*100

projections.df$cory_percentchange2

#Thunnus spp prey gam
THUNGam2<- gam(THUNNUS ~ s(Larvaceans, k=3, m=1) +s(Cyclopoids, k=3),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(THUNGam2)

datahist_thunnus2 <- newdatahist_fullprey[, !colnames(newdatahist_fullprey) %in% c("Calanoids","Cladocerans", "lat", "lon")]
predicted_abundances_hist_thunnus2 <- predict(THUNGam2, newdata = datahist_thunnus2, type = "response")

projections.df$thunnus_projections_hist2 <- predicted_abundances_hist_thunnus2
mean(projections.df$thunnus_projections_hist2)

dataproj_thunnus2 <- newdataproj_fullprey[, !colnames(newdataproj_fullprey)%in% c("Calanoids","Cladocerans", "lat", "lon")]
predicted_abundances_proj_thunnus2 <- predict(THUNGam2, newdata = dataproj_thunnus2, type = "response")

projections.df$thunnus_projections2 <- predicted_abundances_proj_thunnus2
mean(projections.df$thunnus_projections2)

projections.df<- projections.df %>%
  mutate(thunnus_anom2 = thunnus_projections2 -thunnus_projections_hist2)

center_hist_thunnus2 <- c(
  weighted.mean(projections.df$x, projections.df$thunnus_projections_hist2, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$thunnus_projections_hist2, na.rm = TRUE)
)
center_hist_thunnus2

center_proj_thunnus2 <- c(
  weighted.mean(projections.df$x, projections.df$thunnus_projections2, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$thunnus_projections2, na.rm = TRUE)
)
center_proj_thunnus2

distance_meters_thunnus2 <- distHaversine(center_hist_thunnus2, center_proj_thunnus2)
print(distance_km_thunnus2<- distance_meters_thunnus2 / 1000)
print(initial_bearing_thunnus2 <- bearing(center_hist_thunnus2, center_proj_thunnus2))

projections.df$thunnus_percentchange2<-
  +   ((projections.df$thunnus_projections2-projections.df$thunnus_projections_hist2)/projections.df$thunnus_projections_hist2)*100

projections.df$thunnus_percentchange2

#Thunnus thynnus
THYNGam2<- gam(THUNNUSTHYNNUS~s(Larvaceans, k=3) +s(Cladocerans, k=4, m=1) +s(lon, lat)+s(Cyclopoids,k=3),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(THYNGam2)

datahist_thynnus2 <- newdatahist_fullprey[, !colnames(newdatahist_fullprey) %in% c("Calanoids")]
predicted_abundances_hist_thynnus2 <- predict(THYNGam2, newdata = datahist_thynnus2, type = "response")

projections.df$thynnus_projections_hist2 <- predicted_abundances_hist_thynnus2
mean(projections.df$thynnus_projections_hist2)

dataproj_thynnus2 <- newdataproj_fullprey[, !colnames(newdataproj_fullprey)%in% c("Calanoids")]
predicted_abundances_proj_thynnus2 <- predict(THYNGam2, newdata = dataproj_thynnus2, type = "response")

projections.df$thynnus_projections2 <- predicted_abundances_proj_thynnus2
mean(projections.df$thynnus_projections2)

projections.df<- projections.df %>%
  mutate(thynnus_anom2 = thynnus_projections2 -thynnus_projections_hist2)

center_hist_thynnus2 <- c(
  weighted.mean(projections.df$x, projections.df$thynnus_projections_hist2, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$thynnus_projections_hist2, na.rm = TRUE)
)
center_hist_thynnus2

center_proj_thynnus2 <- c(
  weighted.mean(projections.df$x, projections.df$thynnus_projections2, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$thynnus_projections2, na.rm = TRUE)
)
center_proj_thynnus2

distance_meters_thynnus2 <- distHaversine(center_hist_thynnus2, center_proj_thynnus2)
print(distance_km_thynnus2<- distance_meters_thynnus2 / 1000)
print(initial_bearing_thynnus2 <- bearing(center_hist_thynnus2, center_proj_thynnus2))


projections.df$thynnus_percentchange2<-
  +   ((projections.df$thynnus_projections2-projections.df$thynnus_projections_hist2)/projections.df$thynnus_projections_hist2)*100

projections.df$thynnus_percentchange2


###skip jack tuna gam 
KATSUGam2<- gam(KATSUWONUSPELAMIS~ s(Larvaceans) + s(lon, lat),
               data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(KATSUGam2)

datahist_kat2 <- newdatahist_fullprey[, !colnames(newdatahist_fullprey) %in% c("Calanoids","Cyclopoids","Cladocerans")]
predicted_abundances_hist_kat2 <- predict(KATSUGam2, newdata = datahist_kat2, type = "response")

projections.df$kat_projections_hist2 <- predicted_abundances_hist_kat2
mean(projections.df$kat_projections_hist2)

dataproj_kat2<- newdataproj_fullprey[, !colnames(newdataproj_fullprey)%in% c("Calanoids","Cyclopoids","Cladocerans")]
predicted_abundances_proj_kat2 <- predict(KATSUGam2, newdata = dataproj_kat2, type = "response")

projections.df$kat_projections2 <- predicted_abundances_proj_kat2
mean(projections.df$kat_projections2)

center_hist_kat2 <- c(
  weighted.mean(projections.df$x, projections.df$kat_projections_hist2, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$kat_projections_hist2, na.rm = TRUE)
)
center_hist_kat2

center_proj_kat2 <- c(
  weighted.mean(projections.df$x, projections.df$kat_projections2, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$kat_projections2, na.rm = TRUE)
)
center_proj_kat2

distance_meters_kat2 <- distHaversine(center_hist_kat2, center_proj_kat2)
print(distance_km_kat2<- distance_meters_kat2 / 1000)
print(initial_bearing_kat2 <- bearing(center_hist_kat2, center_proj_kat2))

projections.df$kat_percentchange2<-
  +   ((projections.df$kat_projections2-projections.df$kat_projections_hist2)/projections.df$kat_projections_hist2)*100

projections.df$kat_percentchange2


#### frigate gam
AUXGam2<- gam(AUXIS~ s(Calanoids, k=3,m=1) + s(Larvaceans, k=3)+
               s(lon, lat),
             data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(AUXGam2)

datahist_aux2 <- newdatahist_fullprey[, !colnames(newdatahist_fullprey) %in% c("Cyclopoids", "Cladocerans")]
predicted_abundances_hist_aux2 <- predict(AUXGam2, newdata = datahist_aux2, type = "response")

projections.df$aux_projections_hist2 <- predicted_abundances_hist_aux2
mean(projections.df$aux_projections_hist2)

dataproj_aux2<- newdataproj_fullprey[, !colnames(newdataproj_fullprey)%in% c("Cyclopoids","Cladocerans")]
predicted_abundances_proj_aux2 <- predict(AUXGam2, newdata = dataproj_aux2, type = "response")

projections.df$aux_projections2 <- predicted_abundances_proj_aux2
mean(projections.df$aux_projections2)

center_hist_aux2 <- c(
  weighted.mean(projections.df$x, projections.df$aux_projections_hist2, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$aux_projections_hist2, na.rm = TRUE)
)
center_hist_aux2

center_proj_aux2 <- c(
  weighted.mean(projections.df$x, projections.df$aux_projections2, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$aux_projections2, na.rm = TRUE)
)
center_proj_aux2

distance_meters_aux2 <- distHaversine(center_hist_aux2, center_proj_aux2)
print(distance_km_aux2<- distance_meters_aux2 / 1000)
print(initial_bearing_aux2 <- bearing(center_hist_aux2, center_proj_aux2))

projections.df$aux_percentchange2<-
  +   ((projections.df$aux_projections2-projections.df$aux_projections_hist2)/projections.df$aux_projections_hist2)*100

projections.df$aux_percentchange2

#### tunny gam
EUTHGam2 <- gam(EUTHYNNUSALLETTERATUS ~ s(Calanoids, k=4, m=1) + s(Cyclopoids, k=3, m=1) + s(Larvaceans, k=3) + s(lon, lat),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(EUTHGam2)

datahist_eut2 <- newdatahist_fullprey[, !colnames(newdatahist_fullprey) %in% c("Cladocerans")]
predicted_abundances_hist_eut2 <- predict(EUTHGam2, newdata = datahist_eut2, type = "response")

projections.df$eut_projections_hist2 <- predicted_abundances_hist_eut2
mean(projections.df$eut_projections_hist2)

dataproj_eut2<- newdataproj_fullprey[, !colnames(newdataproj_fullprey)%in% c("Cladocerans")]
predicted_abundances_proj_eut2 <- predict(EUTHGam2, newdata = dataproj_eut2, type = "response")

projections.df$eut_projections2 <- predicted_abundances_proj_eut2
mean(projections.df$eut_projections2)

center_hist_eut2<- c(
  weighted.mean(projections.df$x, projections.df$eut_projections_hist2, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$eut_projections_hist2, na.rm = TRUE)
)
center_hist_eut2

center_proj_eut2 <- c(
  weighted.mean(projections.df$x, projections.df$eut_projections2, na.rm = TRUE),
  weighted.mean(projections.df$y, projections.df$eut_projections2, na.rm = TRUE)
)
center_proj_eut2

distance_meters_eut2 <- distHaversine(center_hist_eut2, center_proj_eut2)
print(distance_km_eut2<- distance_meters_eut2 / 1000)
print(initial_bearing_eut2 <- bearing(center_hist_eut2, center_proj_eut2))

projections.df$eut_percentchange2<-
  +   ((projections.df$eut_projections2-projections.df$eut_projections_hist2)/projections.df$eut_projections_hist2)*100

projections.df$eut_percentchange2

#maps 

theme_ray_SDM_fish <- theme(
  legend.position = c(.5, -.58), 
  legend.direction = "horizontal", 
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.background = element_rect(fill = "transparent"),
  legend.text = element_text(size = 5, margin = margin(t = 2)),  
  legend.title = element_text(size = 6, margin = margin(b = 3)),  
  legend.spacing = unit(0.5, "cm"), 
  text = element_text(size = 5),
  legend.key.height = unit(0.5, "lines"),
  panel.background = element_rect(fill = "grey90",color = "black"),
  panel.grid = element_line(color = "grey90", size = 0.8),
  legend.box.margin = margin(-10, 0, 0, 0) ,
  legend.box.spacing = unit(-2, "cm"),
  axis.ticks.length = unit(0.01, "cm"),
  plot.margin = unit(c(-0.15, 0.15, .75, -0.01), "cm"),
  panel.border = element_rect(color = "black", fill = NA, size = 0.5, linetype = "solid"))

custom_guide_sdm <- guide_colorbar(
  barwidth = unit(7, "lines"),
  barheight = unit(0.3, "lines"),
  legend.text = element_text(size = 5, margin = margin(t = 4)),  
  legend.title = element_text(size = 6, margin = margin(b = 3)),  
  title.position = "bottom",
  label.position = "top")

#mahi mahi plots
#cory abiotic

common_mincory_abiotic <- min(c(projections.df$cory_projections, projections.df$cory_projections_hist), na.rm = TRUE)
common_maxcory_abiotic <- max(c(projections.df$cory_projections, projections.df$cory_projections_hist), na.rm = TRUE)
common_limitscory_abiotic <- c(common_mincory_abiotic, common_maxcory_abiotic)

projcorySDM_abiotic<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = cory_projections), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitscory_abiotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Mahi-Mahi Abundance (Abiotic)*")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
projcorySDM_abiotic

projcorySDM_abiotic <- ggdraw(projcorySDM_abiotic) +
  draw_label("Future Period (2070-2099)", x = 0.32, y = 0.95, size = 4,color = "aliceblue")

histcorySDM_abiotic <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = cory_projections_hist), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"),limits = c(0, .55),oob = scales::squish) +  
  theme_ray_SDM_fish+
  labs(fill = "Mahi-Mahi Abundance (Abiotic)*")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
histcorySDM_abiotic

histcorySDM_abiotic<- ggdraw(histcorySDM_abiotic) +
  draw_label("Historical Period (1985-2014)", x = 0.32, y = 0.95, size = 4.0, color = "aliceblue")

#cory biotic
common_mincory_biotic <- min(c(projections.df$cory_projections2, projections.df$cory_projections_hist2), na.rm = TRUE)
common_maxcory_biotic <- max(c(projections.df$cory_projections2, projections.df$cory_projections_hist2), na.rm = TRUE)
common_limitscory_biotic <- c(common_mincory_biotic, common_maxcory_biotic)

projcorySDM_biotic<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = cory_projections2), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitscory_biotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Mahi-Mahi Abundance (Biotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
projcorySDM_biotic

projcorySDM_biotic <- ggdraw(projcorySDM_biotic) +
  draw_label("Future Period (2070-2099)", x = 0.32, y = 0.95, size = 4,color = "aliceblue")

histcorySDM_biotic <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = cory_projections_hist2), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitscory_biotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Mahi-Mahi Abundance (Biotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
histcorySDM_biotic

histcorySDM_biotic<- ggdraw(histcorySDM_biotic) +
  draw_label("Historical Period (1985-2014)", x = 0.32, y = 0.95, size = 4.0, color = "aliceblue")


#ABFT plots 
#abft abiotic
common_minthynnus_abiotic <- min(c(projections.df$thynnus_projections, projections.df$thynnus_projections_hist), na.rm = TRUE)
common_maxthynnus_abiotic <- max(c(projections.df$thynnus_projections, projections.df$thynnus_projections_hist), na.rm = TRUE)
common_limitsthynnus_abiotic <- c(common_minthynnus_abiotic, common_maxthynnus_abiotic)

projthynnusSDM_abiotic<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = thynnus_projections), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsthynnus_abiotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Alt. Bluefin Tuna Abundance (Abiotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
projthynnusSDM_abiotic

projthynnusSDM_abiotic <- ggdraw(projthynnusSDM_abiotic) +
  draw_label("Future Period (2070-2099)", x = 0.32, y = 0.95, size = 4,color = "aliceblue")

histthynnusSDM_abiotic <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = thynnus_projections_hist), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsthynnus_abiotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Alt. Bluefin Tuna Abundance (Abiotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
histthynnusSDM_abiotic

histthynnusSDM_abiotic <- ggdraw(histthynnusSDM_abiotic) +
  draw_label("Historical Period (1985-2014)", x = 0.32, y = 0.95, size = 4.0, color = "aliceblue")


#abft biotic
common_minthynnus_biotic <- min(c(projections.df$thynnus_projections2, projections.df$thynnus_projections_hist2), na.rm = TRUE)
common_maxthynnus_biotic <- max(c(projections.df$thynnus_projections2, projections.df$thynnus_projections_hist2), na.rm = TRUE)
common_limitsthynnus_biotic <- c(common_minthynnus_biotic, common_maxthynnus_biotic)

projthynnusSDM_biotic<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = thynnus_projections2), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsthynnus_biotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Alt. Bluefin Tuna Abundance (Biotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
projthynnusSDM_biotic

projthynnusSDM_biotic <- ggdraw(projthynnusSDM_biotic) +
  draw_label("Future Period (2070-2099)", x = 0.32, y = 0.95, size = 4,color = "aliceblue")

histthynnusSDM_biotic <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = thynnus_projections_hist2), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsthynnus_biotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Alt. Bluefin Tuna Abundance (Biotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
histthynnusSDM_biotic

histthynnusSDM_biotic <- ggdraw(histthynnusSDM_biotic) +
  draw_label("Historical Period (1985-2014)", x = 0.32, y = 0.95, size = 4.0, color = "aliceblue")


#other tuna plots
#other tuna abiotic
common_minthunnus_abiotic <- min(c(projections.df$thunnus_projections, projections.df$thunnus_projections_hist), na.rm = TRUE)
common_maxthunnus_abiotic <- max(c(projections.df$thunnus_projections, projections.df$thunnus_projections_hist), na.rm = TRUE)
common_limitsthunnus_abiotic <- c(common_minthunnus_abiotic, common_maxthunnus_abiotic)

projthunnusSDM_abiotic<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = thunnus_projections), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsthunnus_abiotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Other True Tuna Abundance (Abiotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
projthunnusSDM_abiotic

projthunnusSDM_abiotic <- ggdraw(projthunnusSDM_abiotic) +
  draw_label("Future Period (2070-2099)", x = 0.32, y = 0.95, size = 4,color = "aliceblue")

histthunnusSDM_abiotic <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = thunnus_projections_hist), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsthunnus_abiotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Other True Tuna Abundance (Abiotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
histthunnusSDM_abiotic

histthunnusSDM_abiotic <- ggdraw(histthunnusSDM_abiotic) +
  draw_label("Historical Period (1985-2014)", x = 0.32, y = 0.95, size = 4.0, color = "aliceblue")

#other tuna biotic
common_minthunnus_biotic <- min(c(projections.df$thunnus_projections2, projections.df$thunnus_projections_hist2), na.rm = TRUE)
common_maxthunnus_biotic <- max(c(projections.df$thunnus_projections2, projections.df$thunnus_projections_hist2), na.rm = TRUE)
common_limitsthunnus_biotic <- c(common_minthunnus_biotic, common_maxthunnus_biotic)

projthunnusSDM_biotic<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = thunnus_projections2), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsthunnus_biotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Other True Tuna Abundance (Biotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
projthunnusSDM_biotic

projthunnusSDM_biotic <- ggdraw(projthunnusSDM_biotic) +
  draw_label("Future Period (2070-2099)", x = 0.32, y = 0.95, size = 4,color = "aliceblue")

histthunnusSDM_biotic <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = thunnus_projections_hist2), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsthunnus_biotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Other True Tuna Abundance (Biotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
histthunnusSDM_biotic

histthunnusSDM_biotic <- ggdraw(histthunnusSDM_biotic) +
  draw_label("Historical Period (1985-2014)", x = 0.32, y = 0.95, size = 4.0, color = "aliceblue")

#billfish 
#billfish abiotic
common_minbillfish_abiotic <- min(c(projections.df$billfish_projections, projections.df$billfish_projections_hist), na.rm = TRUE)
common_maxbillfish_abiotic <- max(c(projections.df$billfish_projections, projections.df$billfish_projections_hist), na.rm = TRUE)
common_limitsbillfish_abiotic <- c(common_minbillfish_abiotic, common_maxbillfish_abiotic)

projbillSDM_abiotic<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = billfish_projections), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsbillfish_abiotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Billfish Abundance (Abiotic)*")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
projbillSDM_abiotic

projbillSDM_abiotic <- ggdraw(projbillSDM_abiotic) +
  draw_label("Future Period (2070-2099)", x = 0.32, y = 0.95, size = 4,color = "aliceblue")

histbillSDM_abiotic <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = billfish_projections_hist), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"),limits = c(0, .4),oob = scales::squish) +  
  theme_ray_SDM_fish+
  labs(fill = "Billfish Abundance (Abiotic)*")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
histbillSDM_abiotic

histbillSDM_abiotic <- ggdraw(histbillSDM_abiotic) +
  draw_label("Historical Period (1985-2014)", x = 0.32, y = 0.95, size = 4.0, color = "aliceblue")

#billfish biotic
common_minbillfish_biotic <- min(c(projections.df$billfish_projections2, projections.df$billfish_projections_hist2), na.rm = TRUE)
common_maxbillfish_biotic <- max(c(projections.df$billfish_projections2, projections.df$billfish_projections_hist2), na.rm = TRUE)
common_limitsbillfish_biotic <- c(common_minbillfish_biotic, common_maxbillfish_biotic)

projbillSDM_biotic<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = billfish_projections2), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsbillfish_biotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Billfish Abundance (Biotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
projbillSDM_biotic

projbillSDM_biotic <- ggdraw(projbillSDM_biotic) +
  draw_label("Future Period (2070-2099)", x = 0.32, y = 0.95, size = 4,color = "aliceblue")

histbillSDM_biotic <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = billfish_projections_hist2), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsbillfish_biotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Billfish Abundance (Biotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
histbillSDM_biotic

histbillSDM_biotic <- ggdraw(histbillSDM_biotic) +
  draw_label("Historical Period (1985-2014)", x = 0.32, y = 0.95, size = 4.0, color = "aliceblue")


#tunny
#tunny abiotic
common_mineut_abiotic <- min(c(projections.df$eut_projections, projections.df$eut_projections_hist), na.rm = TRUE)
common_maxeut_abiotic <- max(c(projections.df$eut_projections, projections.df$eut_projections_hist), na.rm = TRUE)
common_limitseut_abiotic <- c(common_mineut_abiotic, common_maxeut_abiotic)

projeutSDM_abiotic<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = eut_projections), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitseut_abiotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Little Tunny Abundance (Abiotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
projeutSDM_abiotic

projeutSDM_abiotic <- ggdraw(projeutSDM_abiotic) +
  draw_label("Future Period (2070-2099)", x = 0.32, y = 0.95, size = 4,color = "aliceblue")

histeutSDM_abiotic <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = eut_projections_hist), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitseut_abiotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Little Tunny Abundance (Abiotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
histeutSDM_abiotic

histeutSDM_abiotic<- ggdraw(histeutSDM_abiotic) +
  draw_label("Historical Period (1985-2014)", x = 0.32, y = 0.95, size = 4.0, color = "aliceblue")

#tunny biotic
common_mineut_biotic <- min(c(projections.df$eut_projections2, projections.df$eut_projections_hist2), na.rm = TRUE)
common_maxeut_biotic <- max(c(projections.df$eut_projections2, projections.df$eut_projections_hist2), na.rm = TRUE)
common_limitseut_biotic <- c(common_mineut_biotic, common_maxeut_biotic)

projeutSDM_biotic<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = eut_projections2), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitseut_biotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Little Tunny Abundance (Biotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
projeutSDM_biotic

projeutSDM_biotic <- ggdraw(projeutSDM_biotic) +
  draw_label("Future Period (2070-2099)", x = 0.32, y = 0.95, size = 4,color = "aliceblue")

histeutSDM_biotic <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = eut_projections_hist2), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitseut_biotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Little Tunny Abundance (Biotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
histeutSDM_biotic

histeutSDM_biotic<- ggdraw(histeutSDM_biotic) +
  draw_label("Historical Period (1985-2014)", x = 0.32, y = 0.95, size = 4.0, color = "aliceblue")



#skipjack
#skipjack abiotic
common_minkat_abiotic <- min(c(projections.df$kat_projections, projections.df$kat_projections_hist), na.rm = TRUE)
common_maxkat_abiotic <- max(c(projections.df$kat_projections, projections.df$kat_projections_hist), na.rm = TRUE)
common_limitskat_abiotic <- c(common_minkat_abiotic, common_maxkat_abiotic)

projkatSDM_abiotic<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = kat_projections), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitskat_abiotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Skipjack Tuna Abundance (Abiotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
projkatSDM_abiotic

projkatSDM_abiotic <- ggdraw(projkatSDM_abiotic) +
  draw_label("Future Period (2070-2099)", x = 0.32, y = 0.95, size = 4,color = "aliceblue")

histkatSDM_abiotic <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = kat_projections_hist), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitskat_abiotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Skipjack Tuna Abundance (Abiotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
histkatSDM_abiotic

histkatSDM_abiotic <- ggdraw(histkatSDM_abiotic) +
  draw_label("Historical Period (1985-2014)", x = 0.32, y = 0.95, size = 4.0, color = "aliceblue")

#skipjack biotic
common_minkat_biotic <- min(c(projections.df$kat_projections2, projections.df$kat_projections_hist2), na.rm = TRUE)
common_maxkat_biotic <- max(c(projections.df$kat_projections2, projections.df$kat_projections_hist2), na.rm = TRUE)
common_limitskat_biotic <- c(common_minkat_biotic, common_maxkat_biotic)

projkatSDM_biotic<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = kat_projections2), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitskat_biotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Skipjack Tuna Abundance (Biotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
projkatSDM_biotic

projkatSDM_biotic <- ggdraw(projkatSDM_biotic) +
  draw_label("Future Period (2070-2099)", x = 0.32, y = 0.95, size = 4,color = "aliceblue")

histkatSDM_biotic <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = kat_projections_hist2), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitskat_biotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Skipjack Tuna Abundance (Biotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
histkatSDM_biotic

histkatSDM_biotic <- ggdraw(histkatSDM_biotic) +
  draw_label("Historical Period (1985-2014)", x = 0.32, y = 0.95, size = 4.0, color = "aliceblue")


#frigate
#frigate abiotic
common_minaux_abiotic <- min(c(projections.df$aux_projections, projections.df$aux_projections_hist), na.rm = TRUE)
common_maxaux_abiotic <- max(c(projections.df$aux_projections, projections.df$aux_projections_hist), na.rm = TRUE)
common_limitsaux_abiotic <- c(common_minaux_abiotic, common_maxaux_abiotic)

projauxSDM_abiotic<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = aux_projections), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsaux_abiotic, oob = scales::squish, guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Frigate Tuna Abundance (Abiotic)")+
guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)), label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
projauxSDM_abiotic

projauxSDM_abiotic <- ggdraw(projauxSDM_abiotic) +
  draw_label("Future Period (2070-2099)", x = 0.32, y = 0.95, size = 4,color = "aliceblue")

histauxSDM_abiotic <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = aux_projections_hist), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsaux_abiotic, oob = scales::squish,guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Frigate Tuna Abundance (Abiotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)),
  label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
histauxSDM_abiotic

histauxSDM_abiotic <- ggdraw(histauxSDM_abiotic) +
  draw_label("Historical Period (1985-2014)", x = 0.32, y = 0.95, size = 4.0, color = "aliceblue")

#frigate biotic
common_minaux_biotic <- min(c(projections.df$aux_projections2, projections.df$aux_projections_hist2), na.rm = TRUE)
common_maxaux_biotic <- max(c(projections.df$aux_projections2, projections.df$aux_projections_hist2), na.rm = TRUE)
common_limitsaux_biotic <- c(common_minaux_biotic, common_maxaux_biotic)

projauxSDM_biotic<-ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = aux_projections2), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsaux_biotic, oob = scales::squish,guide = custom_guide_sdm) +
  theme_ray_SDM_fish+
  labs(fill = "Frigate Tuna Abundance (Biotic)")+
 guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)),  
  label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines")))
projauxSDM_biotic

projauxSDM_biotic <- ggdraw(projauxSDM_biotic) +
  draw_label("Future Period (2070-2099)", x = 0.32, y = 0.95, size = 4,color = "aliceblue")

histauxSDM_biotic <- ggplot() + 
  geom_raster(data = projections.df, aes(x = x, y = y, fill = aux_projections_hist2), interpolate = TRUE, na.rm = TRUE) +
  geom_sf(data = country, fill = "gray33") +
  coord_sf(xlim = c(-100.15, -79.12), ylim = c(30.85, 23.37)) +
  theme_ray_SDM_fish+
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = common_limitsaux_biotic, oob = scales::squish) +
  labs(fill = "Frigate Tuna Abundance (Biotic)")+
  guides(fill = guide_colorbar(title.position = "top",label.position = "bottom",  title.theme = element_text(margin = margin(b = 1)),  
    label.theme = element_text(margin = margin(t = 1)), barwidth = unit(7, "lines"), barheight = unit(0.3, "lines"),))
histauxSDM_biotic

histauxSDM_biotic <- ggdraw(histauxSDM_biotic) +
  draw_label("Historical Period (1985-2014)", x = 0.32, y = 0.95, size = 4.0, color = "aliceblue")


### final fish prey plot

fishSDMS<-grid.arrange(
  histthunnusSDM_biotic, projthunnusSDM_biotic, histthunnusSDM_abiotic, projthunnusSDM_abiotic,
  histthynnusSDM_biotic, projthynnusSDM_biotic, histthynnusSDM_abiotic, projthynnusSDM_abiotic,
  histkatSDM_biotic, projkatSDM_biotic, histkatSDM_abiotic, projkatSDM_abiotic,
  histbillSDM_biotic, projbillSDM_biotic, histbillSDM_abiotic, projbillSDM_abiotic,
  histeutSDM_biotic, projeutSDM_biotic, histeutSDM_abiotic, projeutSDM_abiotic,
  histcorySDM_biotic, projcorySDM_biotic, histcorySDM_abiotic, projcorySDM_abiotic,
  histauxSDM_biotic, projauxSDM_biotic, histauxSDM_abiotic, projauxSDM_abiotic,
  nrow=7, ncol=4)


ggsave("fishSDMS.png", fishSDMS, dpi = 300, bg = "white",
       width = 1961,
       height = 2200,
       units = "px") 

