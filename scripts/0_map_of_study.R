#install.packages(c("ggplot2", "sf", "rnaturalearth", "rnaturalearthdata", "ggspatial"))
#install.packages(c("ggplot2","elevatr", "sf", "terra", "geodata","ggspatial", "dplyr" )
library(ggplot2)
library(sf)
library(ggspatial)
library(dplyr)
library(elevatr)
library(geodata)
library(terra)

# load locality data
locals <- read.csv("formatted_data/locals.csv") %>% filter(site !="Garrapatero")

# Convert to sf object
locals <- st_as_sf(locals, coords = c("longitude", "latitude"), crs = 4326)

# Define a bounding box around the Santa Cruz Island
santa_cruz_bbox <- st_as_sf(st_sfc(
  st_polygon(list(rbind(
    c(-90.6, -0.788),  # SW
    c(-90.6, -0.403),  # NW
    c(-90.11, -0.403),  # NE
    c(-90.11, -0.788),  # SE
    c(-90.6, -0.788)   # close
  ))),
  crs = 4326
))

# Download galapagos shape file, with administrative boundaries that will
# include the island of Santa Cruz
ecu_admin1 <- gadm(country = "ECU", level = 2, path = "map_files/")
galapagos <- ecu_admin1[ecu_admin1$NAME_2 == "Santa Cruz", ]
galapagos <- crop(galapagos, santa_cruz_bbox)
plot(galapagos)
galapagos2 <- sf::st_as_sf(galapagos)

# Get elevation raster (SRTM-based), resolution depends on zoom level (z)
# 2. Get elevation data (z = 10 for ~30 m resolution)
santa_cruz_elev <- get_elev_raster(locations = santa_cruz_bbox, z = 10, clip = "locations") %>%
                    mask(., galapagos2)

# Transform raster into dataframe for plotting in ggplot
santa_cruz_elev_df <- as.data.frame(santa_cruz_elev, xy = T)
colnames(santa_cruz_elev_df)[3] <- "elev"

# Turn negative elevations into NAs for cleaner plotting
santa_cruz_elev_df <- santa_cruz_elev_df %>% mutate(elev = ifelse(elev < 0, NA_real_, elev))



# Make a map of elevation --------------------------------------------------------------



# pdf("output_plots/sampling_map.pdf", width = 10, height = 7)
# ggplot() +
#   geom_raster(data = santa_cruz_elev_df, aes(x = x, y = y, fill = elev))+
#   scale_fill_viridis_b(limits = c(0, 1000), breaks = seq(from = 0, to = 1000, by = 200),
#                        option = "magma", direction = -1, end = 1,
#                        na.value = NA) +
#   geom_sf(data = locals, aes(shape = habitat, size = 2, lwd = 1.5)) +
#   scale_shape_manual(values = c("arid" = 1, "natural" = 0, "fruit_vegetable" = 8,
#                                 "coffee" = 10, "pasture" = 7 ),
#                      labels = c("natural" = "Forest","arid" = "Coast/Urban",  "coffee" = "Coffee",
#                                 "fruit_vegetable" = "Fruit/Veg", "pasture" = "Pasture")) +
#   geom_sf(data = galapagos2, fill = NA) +
#   labs( x = "Longitude", y = "Latitude", fill = "Elevation (m)", shape = "Habitat") +
#   guides(size = "none", linewidth = "none",
#          shape = guide_legend(override.aes = list(size = 3.5))) +
#   theme(text = element_text(size = 18)) +
#   annotation_scale(location = "bl", width_hint = 0.3) +
#   annotation_north_arrow(location = "tl", which_north = "true")
# dev.off()



# would precip be better than elevation?  Yes ---------------------------------

#Download WorldClim bioclimatic variables (e.g. resolution = 10 minutes)
bioclim <- worldclim_global(var = "bio", res = 0.5, path = "map_files/", country = "ECU")

#Crop the raster stack to the Santa Cruz bounding box
#bbox_ext <- ext(bbox)
bioclim_crop <- crop(bioclim, santa_cruz_bbox) %>% mask(., galapagos2)
bioclim_crop2 <- mask(bioclim_crop, galapagos)

# Step 4: Extract annual precipitation (bio 12 )
precip_annual <- bioclim_crop$wc2.1_30s_bio_12 %>%
  crop(., galapagos) %>%
  mask(., galapagos, touches = TRUE)
precip_annual_df <-  as.data.frame(precip_annual, xy = T)
colnames(precip_annual_df)[3] <- "precip"


# Rasterize the vector to match raster resolution
vector_rast <- rasterize(galapagos, precip_annual)

# Multiply the raster and mask layer to get better overlay
r_clean <- precip_annual * vector_rast


pdf("output_plots/precip_map.pdf", width = 10)
ggplot() +
  # set plot bounds
  geom_sf(data = galapagos2, fill = NA) +

  # add precipitation layer
  geom_raster(data = precip_annual_df, aes(x = x, y = y, fill = precip))+
  scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.1) +

  # Add elevation contour lines
  geom_contour(data = santa_cruz_elev_df, aes(x = x, y = y, z = elev),
               breaks = seq(0, 1000, by = 200), color = "grey30") +
  # Add contour labels
  annotate("text", x = -90.270, y = -.587, label = "200", size = 4, color = "grey30", angle = 300) +
  annotate("text", x = -90.290, y = -.618, label = "400", size = 4, color = "grey30", angle = 300) +
  annotate("text", x = -90.310, y = -.635, label = "600", size = 4, color = "grey30", angle = 300) +
  annotate("text", x = -90.320, y = -.641, label = "800", size = 4, color = "grey30", angle = 300) +

  # Add sampling locations
  geom_sf(data = locals, aes(shape = habitat, size = 2, lwd = 1.5), color = "#f39c12") +
  scale_shape_manual(values = c("arid" = 16,
                                "natural" = 18,
                                "fruit_vegetable" = 8,
                                "coffee" = 10,
                                "pasture" = 0 ),
                     labels = c("natural" = "Forest",
                                "arid" = "Coast/Urban",
                                "coffee" = "Coffee",
                                "fruit_vegetable" = "Fruit/Veg",
                                "pasture" = "Pasture")) +
  # Add border around island
  geom_sf(data = galapagos2, fill = NA, linewidth = .8) +

  # Clean up axes and legend labels
  labs( x = "Longitude", y = "Latitude", fill = "Annual\nprecipitation (mm)", shape = "Habitat") +
  guides(size = "none", linewidth = "none",
         shape = guide_legend(override.aes = list(size = 3.5))) +
  theme(text = element_text(size = 18)) +

  # add arrow and scale
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true")

dev.off()




# For picking coordinates of map
#install.packages("plotly")
#library(plotly)


#p <- ggplot() +
#  geom_sf(data = galapagos2, fill = NA) +
#  geom_raster(data = precip_annual_df, aes(x = x, y = y, fill = precip))+
#  scale_fill_viridis_c(option = "rocket", na.value = NA, direction = -1, begin = .1) +
#  geom_contour(data = santa_cruz_elev_df, aes(x = x, y = y, z = elev),
#               breaks = seq(0, 1000, by = 200), color = "grey30")
#ggplotly(p)
