###########################################################
####                                                   ####
#### Script to clean and explore census commuting data ####
####                                                   ####
###########################################################


#### Load packages ####
pacman::p_load(tidyverse, data.table, mgcv, sf, spdep, spam, igraph, 
               sfnetworks, lwgeom)

## Create function opposite of %in%
'%!in%' <- function(x,y)!('%in%'(x,y))



#### Load data ####
## Brazil shapefile 
shp <- read_sf("Data/BRMUE250GC_SIR.shp") %>% 
  mutate(municip_code_ibge = as.numeric(CD_GEOCMU))

## Dengue data
df_deng <- fread("Data/dengue_year.csv")


## Census commuting data
df_move <- fread("Data/all_FUs-redistributed_mobility_matrix.csv") %>% 
  # replace spaces from variable names with _
  rename_with(~tolower(str_replace_all(., "\\s+", "_"))) %>% 
  # Remove international + missing links
  filter(origin_country == "BRASIL" & destination_country == "BRASIL" &
           substr(destination_geocode, 3, 7) != 99999) %>% 
  mutate(destination_geocode = as.numeric(destination_geocode))


## Load file to convert 2010 municipalities into 'parent' municipalities
parent_municip <- fread("Data/parent_municip.csv") %>% 
  dplyr::select(municip_code_ibge, municip_parent_code)


## Convert shapefile to 'parent' municipalities
shp_parent <- right_join(shp, parent_municip, 
                         by = "municip_code_ibge") %>% 
  group_by(municip_parent_code) %>% 
  summarise(geometry = st_union(geometry)) %>% 
  ungroup() %>% 
  rename(municip_code_ibge = municip_parent_code)


# Save converted shapefile for analysis
write_rds(shp_parent, file = "Output/shp_parent.RDS")


#### Clean movement data ####
## Combine population in origin cities between 'parent' municipalities
population_conv <- df_move %>% 
  dplyr::select(origin_geocode, population) %>% 
  # Aggregate to 'parent' municipalities
  left_join(., parent_municip, by = c("origin_geocode" = "municip_code_ibge")) %>% 
  # Avoid duplicates in sums
  group_by(origin_geocode) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(municip_parent_code) %>% 
  # Combine population counts
  summarise(pop_parent = sum(population),
            n = n()) %>% 
ungroup() 


## Aggregate numbers moving between 'parent' municipalities
df_move_conv <- df_move %>% 
  # Aggregate origin cities
  left_join(., parent_municip, by = c("origin_geocode" = "municip_code_ibge")) %>% 
  rename(origin_parent_code = municip_parent_code) %>% 
  group_by(origin_parent_code, destination_geocode) %>% 
  summarise(total_parent = sum(total)) %>% 
  ungroup() %>% 
  # Aggregate destination cities
  left_join(., parent_municip, by = c("destination_geocode" = "municip_code_ibge")) %>% 
  rename(dest_parent_code = municip_parent_code) %>% 
  group_by(origin_parent_code, dest_parent_code) %>% 
  summarise(total_parent = sum(total_parent)) %>% 
  ungroup() %>% 
  # Remove duplicates
  filter(origin_parent_code != dest_parent_code) %>% 
  # Add population
  full_join(., population_conv, by = c("origin_parent_code" = "municip_parent_code")) %>%
  # Estimate the rate of movement (proportion of population moving)
  mutate(density = (total_parent/pop_parent),
         # Convert to 'distance' i.e., the proportion not moving
         density_1 = 1 - (total_parent/pop_parent)) %>% 
  rename(origin_code = origin_parent_code,
         dest_code = dest_parent_code)



#### Plot links across Brazil ####
## Add coordinates to movement data 
df_coord <- group_by(df_deng, municip_code_ibge, lon, lat) %>% 
  # Only keep first row of each municipality
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  dplyr::select(municip_code_ibge, lon, lat)


df_move_coord <-  df_move_conv %>% 
  full_join(., df_coord, by = c("origin_code" = "municip_code_ibge")) %>% 
  full_join(., df_coord, by = c("dest_code" = "municip_code_ibge"),
            suffix = c("_origin", "_dest"))


# Plot connections to/from Brasilia
df_brasilia <- df_move_coord %>% 
  filter(origin_code == 5300108 | dest_code == 5300108)

brasilia_connect <- ggplot( ) +
  geom_sf(data = shp, lwd = .05, fill = "White") + 
  geom_segment(data = df_move_coord %>% 
                 filter(origin_code == 5300108 | dest_code == 5300108),
               aes(x = lon_origin, y = lat_origin,
                   xend = lon_dest, yend = lat_dest, colour = density*10),
               lwd = .5) + 
  scale_colour_viridis_c(name = "% moving", direction = -1) +
  expand_limits(colour = c(0, 3)) +
  theme_void() 

ggsave(brasilia_connect, filename = "Output/brasilia_connect.png")


# Plot connections with > 1% population moving to/from Brasilia
brasilia_connect_zoom <- ggplot( ) +
  geom_sf(data = shp, lwd = .05, fill = "White") + 
  geom_segment(data = df_move_coord %>% 
                 filter(origin_code == 5300108 | dest_code == 5300108,
                        density > .01),
               aes(x = lon_origin, y = lat_origin,
                   xend = lon_dest, yend = lat_dest, colour = density*10),
               lwd = .5) + 
  scale_colour_viridis_c(name = "% moving", direction = -1) +
  expand_limits(colour = c(0, 3)) +
  # Only contains close municipalities so zoom in 
  coord_sf(x = c(-50, -47), y = c(-17, -15)) +
  theme_void() 

ggsave(brasilia_connect_zoom, filename = "Output/brasilia_connect_zoom.png")



# Plot connections to/from Rio Branco
riobranco_connect <- ggplot( ) +
  geom_sf(data = shp, lwd = .05, fill = "white") + 
  geom_segment(data = df_move_coord %>% 
                 filter(origin_code == 1200401 | dest_code == 1200401),
               aes(x = lon_origin, y = lat_origin, 
                   xend = lon_dest, yend = lat_dest, colour = density*10),
               lwd = .5) + 
  scale_colour_viridis_c(name = "% moving", direction = -1) +
  expand_limits(colour = c(0, 3)) +
  theme_void()

ggsave(riobranco_connect, file = "Output/riobranco_connect.png")


# Plot connections with > 1% population moving to/from Rio Branco
riobranco_connect_zoom <- ggplot( ) +
  geom_sf(data = shp, lwd = .05, fill = "white") + 
  geom_segment(data = df_move_coord %>% 
                 filter(origin_code == 1200401 | dest_code == 1200401, 
                        density > .01),
               aes(x = lon_origin, y = lat_origin, 
                   xend = lon_dest, yend = lat_dest, colour = density*10),
               lwd = .5) + 
  scale_colour_viridis_c(name = "% moving", direction = -1) +
  expand_limits(colour = c(0, 3)) +
  coord_sf(xlim = c(-70, -67), ylim = c(-11, -9)) +
  theme_void()

ggsave(riobranco_connect_zoom, file = "Output/riobranco_connect_zoom.png")


# Plot connections to/from Porto Alegre
portoalegre_connect <- ggplot( ) +
  geom_sf(data = shp, lwd = .05, fill = "white") + 
  geom_segment(data = df_move_coord %>% 
                 filter(origin_code == 4314902 | dest_code == 4314902),
               aes(x = lon_origin, y = lat_origin, 
                   xend = lon_dest, yend = lat_dest, colour = density*10),
               lwd = .5) + 
  scale_colour_viridis_c(name = "% moving", direction = -1) +
  expand_limits(colour = c(0, 3)) +
  theme_void()

ggsave(portoalegre_connect, filename = "Output/portoalegre_connect.png")


# Plot connections with > 1% population moving to/from Porto Alegre
portoalegre_connect_zoom <- ggplot( ) +
  geom_sf(data = shp, lwd = .05, fill = "white") + 
  geom_segment(data = df_move_coord %>% 
                 filter((origin_code == 4314902 | dest_code == 4314902),
                        density > .01),
               aes(x = lon_origin, y = lat_origin, 
                   xend = lon_dest, yend = lat_dest, colour = density*10),
               lwd = .5) + 
  scale_colour_viridis_c(name = "% moving", direction = -1) +
  expand_limits(fill = c(0, 3)) +
  coord_sf(xlim = c(-53, -50), ylim = c(-31, -29)) +
  theme_void()

ggsave(portoalegre_connect_zoom, filename = "Output/portoalegre_connect_zoom.png")


#### Create movement coordinates based on density of movement  ####
## Convert pairs of municipalities to matrix with codes as column and row names
move_mat <- df_move_conv %>% 
  dplyr::select(origin_code, dest_code, density_1) %>% 
  arrange(origin_code, dest_code) %>% 
  pivot_wider(names_from = origin_code, values_from = density_1,
              values_fill = 0) %>% 
  column_to_rownames(var = "dest_code") 



## Apply multidimensional scaling to convert movement matrix into 2-d coordinates
msd_scale <- cmdscale(move_mat, eig = T, k = 2)


# Convert into a df with 'connectivity coordinates' per city
connect_coords <-  as.data.frame(msd_scale$points)
names(connect_coords) <- c("connect_coord1", "connect_coord2")
connect_coords$municip_code_ibge <- as.numeric(rownames(connect_coords))

# Plot to view coordinate system
human_coord_plot <- ggplot(data = connect_coords) +
  geom_point(aes(x = connect_coord1, connect_coord2)) +
  labs(x = "Connectivity coordinate 1", y = "Connectivity coordinate 2") +
  expand_limits(x = c(-.6, .6), y = c(-.6, .6)) +
  theme_bw()

ggsave(human_coord_plot, filename = "Output/human_coord_plot.png")


#### Add human movement coordinates to other data and save ####
df <- full_join(df_deng, connect_coords, by = "municip_code_ibge")

## Save full data for analysis
fwrite(df, file = "Data/df_full.csv")


#### Condense data to single year for models ####
df_binom <- df %>% 
  # Convert dengue cases to incidence and outbreak indicator (DIR > 300)
  mutate(DIR = (dengue_year/population) * 10^5,
         outbreak = ifelse(DIR >= 300, 1, 0)) %>% 
  group_by(municip_code_ibge, region_code, region_name, urban10,
           level18_num, lon, lat, connect_coord1, connect_coord2) %>%
  # Return the number of outbreaks over the time period
  summarise(n_outbreaks = sum(outbreak),
            median_suitable = median(months_suitable.era),
            n_yrs = n()) %>% 
  ungroup() 

# Save data for analysis
fwrite(df_binom, file = "data/df_model.csv")

