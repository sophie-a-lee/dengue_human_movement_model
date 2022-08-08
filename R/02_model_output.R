#############################################
####                                     ####
####    Script to explore model output   ####
####                                     ####
#############################################


#### Load packages ####
pacman::p_load(tidyverse, data.table, mgcv, sf, nimble, coda, mgcViz, cowplot,
               pROC)
sf_use_s2(F)

## Set up colour palette for regions 
region_col <- c("#ef476f", "#FFBA08", "#06d6a0", "#118ab2", "#073b4c")


#### Load data ####
## NOTE: This code will only work if 00_human_movement_collate.R has been run
# Dengue outbreak data
df_model <- fread("data/df_model.csv")

# Brazil shapefile
shp <- readRDS("output/shp_parent.RDS")


#### Load model results ####
## NOTE: This code will only work if 01_run_model.R has been run
model_results <- read_rds("Output/mcmc_results.rds")


#### Return proportion of variance per spatial term ####
n <- nrow(df_model)

# Return column numbers with distance-based random effect simulations
u_dist_min <- which(colnames(model_results[[1]]) == "u_dist[1]")
u_dist_max <- which(colnames(model_results[[1]]) == paste0("u_dist[", n, "]"))

# Return column numbers with human movement-based random effect simulations
u_human_min <- which(colnames(model_results[[1]]) == "u_human[1]")
u_human_max <- which(colnames(model_results[[1]]) == paste0("u_human[", n, "]"))

# Return column numbers with unstructured random effect simulations
vmin <- which(colnames(model_results[[1]]) == "v[1]")
vmax <- which(colnames(model_results[[1]]) == paste0("v[", n, "]"))


# Extract simulations of distance-based structured random effect
u_dist_mat <- do.call(rbind, model_results)[,u_dist_min:u_dist_max]
# Estimate the variance of each simulation
u_dist_var <- apply(u_dist_mat, 1, var)

# Extract simulations of human movement-based structured random effect
u_human_mat <- do.call(rbind, model_results)[,u_human_min:u_human_max]
# Estimate the variance of each simulation
u_human_var <- apply(u_human_mat, 1, var)

# Extract simulations of unstructured random effect
v_mat <- do.call(rbind, model_results)[,vmin:vmax]
# Estimate the variance of each simulation
v_var <- apply(v_mat, 1, var)

# Calculate the proportion of variance explained by each term
propn_spat_dist <- u_dist_var/(u_dist_var +  u_human_var + v_var)
propn_spat_human <- u_human_var/(u_dist_var +  u_human_var + v_var)
propn_spat_iid <- v_var/(u_dist_var +  u_human_var + v_var)

re_vars <- data.table(spat_dist_var = u_dist_var,
                      spat_human_var = u_human_var,
                      iid_var = v_var,
                      propn_var_dist = propn_spat_dist,
                      propn_var_human = propn_spat_human,
                      propn_var_iid = propn_spat_iid)


## Plot proportion of variance explained by each spatial term
propn_var_plot <- ggplot(data = re_vars) +
  geom_density(aes(x = propn_var_dist, ..scaled.., fill = "Distance"),
               alpha = .7) +
  geom_density(aes(x = propn_var_human, ..scaled.., 
                   fill = "Human movement"), alpha = .7) +
  geom_density(aes(x = propn_var_iid, ..scaled.., 
                   fill = "Unstructured"), alpha = .7) +
  labs(x = "Proportion of variance", y = "Density") +
  scale_fill_manual("", values = c("#F26430", "#6761A8", "#009B72"),
                    breaks = c("Distance", "Human movement",
                               "Unstructured")) +
  theme_bw()


ggsave(propn_var_plot, filename = "Output/propn_var_brazil.png")


## Return estimates and 95% CI for each proportion
propn_var_table <- data.table(distance = paste0(round(median(re_vars$propn_var_dist), 3), 
                                                " (",
                                                round(quantile(re_vars$propn_var_dist, .025), 3), 
                                                ", ",
                                                round(quantile(re_vars$propn_var_dist, .975), 3), 
                                                ")"),           
                              human = paste0(round(median(re_vars$propn_var_human), 3), 
                                             " (",
                                             round(quantile(re_vars$propn_var_human, .025), 3), 
                                             ", ",
                                             round(quantile(re_vars$propn_var_human, .975), 3), 
                                             ")"),
                              iid = paste0(round(median(re_vars$propn_var_iid), 3), 
                                           " (",
                                           round(quantile(re_vars$propn_var_iid, .025), 3), 
                                           ", ",
                                           round(quantile(re_vars$propn_var_iid, .975), 3), ")"))

fwrite(propn_var_table, file = "Output/propn_var_table.csv")


#### Compare proportion of variance explained by spatial terms between regions ####
## Add municipality names to columns
colnames(u_dist_mat) <- df_model$municip_code_ibge

## Extract variances for each region separately
u_dist_reg_mat <- list()
u_human_reg_mat <- list()
v_reg_mat <- list()


propn_reg_dist <- list()
propn_reg_human <- list()
propn_reg_iid <- list()


for(i in 1:5) {
  u_dist_reg_num <-  which(substr(colnames(u_dist_mat), 1, 1) == i)
  
  u_dist_reg_mat[[i]] <- u_dist_mat[ , u_dist_reg_num]
  u_human_reg_mat[[i]] <- u_human_mat[ , u_dist_reg_num]
  v_reg_mat[[i]] <- v_mat[ , u_dist_reg_num]
  
  u_dist_var <- apply(u_dist_reg_mat[[i]], 1, var)
  u_human_var <- apply(u_human_reg_mat[[i]], 1, var)
  v_regvar <- apply(v_reg_mat[[i]], 1, var)
  
  propn_reg_dist[[i]] <- u_dist_var/(u_dist_var +  u_human_var + v_var)
  propn_reg_human[[i]] <- u_human_var/(u_dist_var +  u_human_var + v_var)
  propn_reg_iid[[i]] <- v_var/(u_dist_var +  u_human_var + v_var)
  
}


# Combine all sources into a single data table
re_vars_reg <- data.table(region = factor(rep(1:5, each = 12000), levels = 1:5,
                                          labels = c("North", "Northeast", "Southeast",
                                                     "South", "Centre-West")),
                          distance = reduce(propn_reg_dist,  c),
                          human = reduce(propn_reg_human,  c),
                          independent = reduce(propn_reg_iid,  c))



## Convert regional proportions into long format to include on same plot
re_vars_reg_long <- re_vars_reg %>% 
  pivot_longer(distance:independent,
               values_to = "propn",
               names_to = "source")


## Plot proportion of spatial variance explained by each term, faceted by source
propn_var_regplot <- ggplot(data = re_vars_reg_long) +
  geom_density(aes(x = propn, ..scaled.., fill = region, group = region),
               alpha = .5) +
  labs(x = "Proportion of variance", y = "Density") +
  scale_fill_manual(name = "Region", values = region_col) +
  facet_grid(rows = vars(source)) +
  theme_bw()


ggsave(propn_var_regplot, filename = "Output/propn_var_regplot.png",
       width = 15, height = 10)


## Return estimates and 95% CI for each proportion
propn_summary <- function(variable) {
  paste0(round(median(variable), 3), " (",
         round(quantile(variable, .025), 3), ", ",
         round(quantile(variable, .975), 3), ")")
}


distance_reg_propn <- c("Distance", 
                        tapply(re_vars_reg$distance, re_vars_reg$region, propn_summary))
human_reg_propn <- c("Human", 
                     tapply(re_vars_reg$human, re_vars_reg$region, propn_summary))
iid_reg_propn <- c("IID", 
                   tapply(re_vars_reg$independent, re_vars_reg$region, propn_summary))

propn_var_regtab  <- data.table(reduce(list(distance_reg_propn, 
                                            human_reg_propn, iid_reg_propn),
                                       rbind))  

fwrite(propn_var_regtab, file = "Output/propn_var_reg_table.csv")


#### Extract estimates of spatial random terms ####
u_dist_mean <- apply(u_dist_mat, 2, mean)
u_human_mean <- apply(u_human_mat, 2, mean)
v_mean <- apply(v_mat, 2, mean)

# Estimate total random term 
re_mean <- apply((u_dist_mat + u_human_mat + v_mat), 2, mean)

# Add to shapefile to plot
shp <- shp %>% 
  mutate(tot_re = re_mean,
         spat_dist_re = u_dist_mean, 
         spat_human_re = u_human_mean,
         iid_re = v_mean)


# Plot estimates of each random term and the sum
tot_re_map <- ggplot(data = filter(shp, municip_code_ibge != 2605459)) +
  geom_sf(aes(fill = tot_re), lwd = .05) +
  scale_fill_gradient2(name = "Total", low = "#4d9221", 
                       mid = "white", high = "#c51b7d", midpoint = 0) +
  expand_limits(fill = c(-7, 7)) +
  theme_void()

dist_re_map <- ggplot(data = filter(shp, municip_code_ibge != 2605459)) +
  geom_sf(aes(fill = spat_dist_re), lwd = .05) +
  scale_fill_gradient2(name = "Distance", low = "#4d9221", 
                       mid = "white", high = "#c51b7d", midpoint = 0) +
  expand_limits(fill = c(-7, 7)) +
  theme_void()

hum_re_map <- ggplot(data = filter(shp, municip_code_ibge != 2605459)) +
  geom_sf(aes(fill = spat_human_re), lwd = .05) +
  scale_fill_gradient2(name = "Human", low = "#4d9221", 
                       mid = "white", high = "#c51b7d", midpoint = 0) +
  expand_limits(fill = c(-2, 2)) +
  theme_void()

iid_re_map <- ggplot(data = filter(shp, municip_code_ibge != 2605459)) +
  geom_sf(aes(fill = iid_re), lwd = .05) +
  scale_fill_gradient2(name = "IID", low = "#4d9221", 
                       mid = "white", high = "#c51b7d", midpoint = 0) +
  expand_limits(fill = c(-3, 3)) +
  theme_void()


## Save plots
ggsave(tot_re_map, filename = "Output/tot_re_map.png")
ggsave(dist_re_map, filename = "Output/dist_re_map.png")
ggsave(hum_re_map, filename = "Output/hum_re_map.png")
ggsave(iid_re_map, filename = "Output/iid_re_map.png")
