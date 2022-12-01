#######################################################################
####                                                               ####
#### Script to fit baseline binomial model to dengue outbreak data ####
####                                                               ####
#######################################################################


#### Load packages ####
pacman::p_load(tidyverse, data.table, mgcv, sf, nimble, coda, mgcViz, cowplot,
               ggthemes, pROC)
sf_use_s2(F)


#### Load data ####
## NOTE: This code will only work if 00_human_movement_collate.R has been run
df <- fread("Data/df_model.csv")



#### Plot number of outbreaks on a map ####
## NOTE: This code will only work if 00_human_movement_collate.R has been run
shp <- read_rds("Output/shp_parent.RDS")


n_outbreak_map <- shp %>% 
  full_join(., df, by = "municip_code_ibge") %>% 
  mutate(n_outbreak_na = ifelse(n_outbreaks == 0, NA, n_outbreaks)) %>% 
  ggplot( ) +
  geom_sf(aes(fill = n_outbreak_na, colour = ""), lwd = 0) +
  scale_fill_gradient2(name = "Number of \noutbreaks", 
                       low = "#f20089", mid = "#240046", high = "#2d00f7",
                       midpoint = 10,
                       # Change colour for municipalities with no outbreaks
                       na.value = "grey") +
  scale_colour_manual(values = NA) +
  guides(colour = guide_legend("No outbreaks")) +
  theme_void()


# Save figure
ggsave(n_outbreak_map, file = "Output/n_outbreak_map.png")


## Count the number of municipalities with/without outbreaks by region
addmargins(table(df$region_name, df$n_outbreaks))


#### Fit baseline binomial model ####
## Formulate spatial smooth splines
jagam_output <- jagam(n_outbreaks/n_yrs ~ s(lon, lat, k = 10, bs = "tp") +
                        s(connect_coord1, connect_coord2, k = 10, bs = "tp"),
                      weights = df$n_yrs,
                      file = "jagam.txt", 
                      # Set prior for smoothing function (in this case, the lon/lat spatial smooth)
                      sp.prior = "gamma", 
                      family = binomial,
                      diagonalize = F, 
                      data = df)


# Extract basis functions to use as linear predictors in the model code
X_dist <- jagam_output$jags.data$X[,2:10]
X_human <- jagam_output$jags.data$X[,11:19]

# Set constants for model (number obs +  number coefficients)
Consts <- list(n = nrow(df), n_yrs = df$n_yrs[1], 
               m_dist = ncol(X_dist), 
               m_human = ncol(X_human), 
               m_tot = ncol(jagam_output$jags.data$X))


## Write model formula
Model <- nimbleCode({ 
  
  # u_dist = spatial smooth term based on distance
  u_dist[1:n] <- X_dist[1:n, ] %*% b[2:(m_tot - m_human)] 
  
  # u_human = spatial smooth term based on human movement
  u_human[1:n] <- X_human[1:n, ] %*% b[(m_tot - m_human + 1):(m_tot)] 
  
  for (i in 1:n) { 
    # y = number of cases
    y[i] ~ dbin(p[i], n_yrs) 
    
    logit(p[i]) <- b[1] + u_dist[i] + u_human[i] + v[i]
    
    # v = iid random effect
    v[i] ~ dnorm(0, sd = sig_re)
    
  } 
  
  # Priors
  # Random effect SD
  sig_re ~ dexp(.1)
  
  # Intercept
  b[1] ~ dnorm(0, sd = 5) 
  
  ## prior for sd(s(lon,lat))
  K1[1:(m_dist),1:(m_dist)] <- 
    S1[1:(m_dist),1:(m_dist)] * lambda[1] + 
    S1[1:(m_dist), (m_dist + 1):(2*m_dist)] * lambda[2]
  
  
  ## prior for sd(s(connect_coord1, connect_coord2))
  K2[1:m_human,1:m_human] <- 
    S2[1:m_human, 1:m_human] * lambda[3] + 
    S2[1:m_human, (m_human + 1):(2*m_human)] * lambda[4]
  
  # Prior for smooth coefficient
  b[2:(m_tot - m_human)] ~ dmnorm(zero[1:m_dist], 
                                  K1[1:m_dist, 1:m_dist]) 
  
  b[(m_tot - m_human + 1):(m_tot)] ~ dmnorm(zero[1:m_human], 
                                            K2[1:m_human, 1:m_human]) 
  
  ## smoothing parameter priors 
  for (i in 1:4) {
    
    lambda[i] ~ T(dgamma(.05,.005), 0, 10)
    
  }
  
} )


# Convert jagam data into data suitable for nimble
nimbleData <- list(y = df$n_outbreaks, 
                   X_dist = X_dist, 
                   X_human = X_human,
                   zero = jagam_output$jags.data$zero,
                   S1 = jagam_output$jags.data$S1,
                   S2 = jagam_output$jags.data$S2)



# Set initial values for MCMC
inits <- list(b = rnorm(ncol(jagam_output$jags.data$X), sd=0.1), 
              v =  rnorm(nrow(df), 1), 
              lambda = rep(1, 4), sig_re = .5)


# Set up model in nimble code
nimbleModel <- nimbleModel(code = Model, name = 'nimbleModel', 
                           constants = Consts, data = nimbleData, 
                           inits = inits)

# Tell model which parameter to estimate and return
MCMCconfig <- configureMCMC(nimbleModel,monitors=c("b", "lambda", "u_dist", 
                                                   "u_human", "v", 
                                                   "sig_re", "p", "y"))

# Configure and run MCMC
modelMCMC <- buildMCMC(MCMCconfig)
compiled_model <- compileNimble(nimbleModel)
compiled_model_MCMC <- compileNimble(modelMCMC, project = nimbleModel) 


results <- runMCMC(compiled_model_MCMC, niter = 600000, nburnin = 200000, 
                   nchains = 3, thin=100, inits=inits, progressBar = T, 
                   samplesAsCodaMCMC = T)



#### Save model results as RData ####
write_rds(model_results, "Output/mcmc_results_base.rds")

