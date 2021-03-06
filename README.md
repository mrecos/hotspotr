
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hotspotr

<!-- badges: start -->

<!-- badges: end -->

The goal of hotspotr is to mimic the ESRI ArcGIS vesrion of “HotSpot
Analysis tool”. It is not a duplication of that tool, but replicates it
in spirit…

## Installation

The development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mrecos/hotspotr")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(hotspotr)
library(gstat)
library(sp)
library(sf)
library(raster)
library(ggplot2)
library(tidyverse)

library(janitor)
library(spdep)
library(velox)
library(dismo)
```

``` r
## basic example code
# http://santiago.begueria.es/2010/10/generating-spatially-correlated-random-fields-with-r/
# unconditional simulations on a 100 x 100 grid using gstat
set.seed(717)

# create structure
xy <- expand.grid(1:100, 1:100)
names(xy) <- c("x","y")

# define the gstat object (spatial model)
g_dummy <- gstat::gstat(formula   = z~1+x+y, 
                 locations = ~x+y, 
                 dummy     = TRUE, 
                 beta      = c(1,0.01,0.005),
                 model     = vgm(psill=0.025, range=15, model='Exp'),
                 nmax      = 20)

# make four simulations based on the stat object
g_pred <- predict(g_dummy, newdata = xy, nsim = 1)
#> [using unconditional Gaussian simulation]

#### Create points
points_xy <- data.frame(x = rnorm(150, 50, 12),
                     y = rnorm(150, 50, 12)) %>% 
  filter(x <= 100 & y <= 100) %>% 
  filter(x >= 0 & y >= 0) %>% 
  mutate(pnt_id = seq(1,n())) %>% 
   st_as_sf(coords = c("x","y"),
            remove = FALSE)

points_xy$sim1 <- raster::extract(x = raster::rasterFromXYZ(g_pred),
                               y = points_xy)

survey_area <- st_buffer(points_xy,
                         dist = 12) %>% 
  st_union() %>% 
  st_convex_hull()

ggplot2::ggplot() +
  geom_raster(data = g_pred, 
              aes(x = x, y = y, fill = sim1),
              interpolate = FALSE) +
  geom_point(data = points_xy, 
             aes(x = x, y = y),
             size = 4,
             color = "black") +
  geom_sf(data = survey_area, 
          fill = NA,
          size = 1.5,
          color = "black") +
  scale_fill_viridis_c() +
  theme_void()
```

<img src="man/figures/README-create_sim_data-1.png" width="100%" />

``` r

# fishnet
analysis_fishnet <- st_make_grid(st_union(survey_area),
                                 cellsize = 10, square = FALSE) %>%
  st_sf() 
analysis_fishnet <- analysis_fishnet[survey_area, , op = st_intersects] %>% 
  mutate(net_id = 1:n())

ggplot2::ggplot() +
  # geom_raster(data = g_pred, 
  #             aes(x = x, y = y, fill = sim1),
  #             interpolate = FALSE) +
  geom_point(data = points_xy, 
             aes(x = x, y = y),
             size = 4,
             color = "black") +
  geom_sf(data = analysis_fishnet, 
          fill = NA,
          size = 1.5,
          color = "purple") +
    geom_sf(data = survey_area, 
          fill = NA,
          size = 1.5,
          color = "black") +
  scale_fill_viridis_c() +
  theme_void()
```

<img src="man/figures/README-create_sim_data-2.png" width="100%" />

``` r
  small_net <- analysis_fishnet %>% 
    mutate(net_id = seq(1:n())) %>% 
    st_join(points_xy) %>% 
    group_by(net_id) %>% 
    summarise(mean_value = mean(sim1, na.rm=TRUE)) %>% 
    mutate(mean_value = ifelse(is.nan(mean_value), NA, mean_value),
           mean_value = ifelse(is.infinite(mean_value), NA, mean_value)) 

  small_net_na_omit <- small_net %>% 
    na.omit() %>% 
    mutate(net_id2 = seq(1:n()))
```

``` r
# loop over bands for flakes FOR EACH SITE 
  # results in list of long-dfs that are 10xlonger than obs for each site
  # each total site obs x 10 (number of different distance bands)
  # total lenght of that is 58260 (at the moment), checks out

  # Gi_results_j <- vector(mode = "list", length = length(unique(dat_join$site_group)))
bands <- seq(5,25,5)
net_coords <- st_coordinates(points_xy)[,c("X","Y")]
results_gistar <- data.frame() #results holder
for(i in seq_along(bands)){
  cat("Looping over band:",bands[i],"\n")
  nb_i <- include.self(dnearneigh(net_coords, 0, bands[i]))
  sac_nb_i <- nb2listw(nb_i, style="B", zero.policy = TRUE) #NULL
  localg_i <- localG(points_xy$sim1, sac_nb_i)
  result_i <- data.frame(pnt_id  = points_xy$pnt_id,
                         gi_star = as.numeric(localg_i),
                         m       = unlist(lapply(sac_nb_i$weights, sum)),
                         nb      = bands[i])
  results_gistar <- rbind(results_gistar, result_i) #long DF
}
#> Looping over band: 5 
#> Looping over band: 10 
#> Looping over band: 15 
#> Looping over band: 20 
#> Looping over band: 25
```

``` r
## z -> p -> bon correct alpha_star -> z-threshold, compare z to threshold
gistar_corrected1 <- results_gistar %>%
 mutate(p_val   = 2*pnorm(-abs(gi_star)),
        sign    = sign(gi_star),
        a_star_95  = 0.05/m,
        a_star_99  = 0.01/m,
        a_star_999 = 0.001/m,
        # no a_star/2 used. without best approximates table 3 in ord/getis
        z_95    = qnorm(a_star_95,  lower.tail = FALSE)*sign(gi_star),
        z_99    = qnorm(a_star_99,  lower.tail = FALSE)*sign(gi_star),
        z_999   = qnorm(a_star_999, lower.tail = FALSE)*sign(gi_star),
        sig_95  = ifelse(abs(gi_star) >= abs(z_95),  1, 0),
        sig_99  = ifelse(abs(gi_star) >= abs(z_99),  1, 0),
        sig_999 = ifelse(abs(gi_star) >= abs(z_999), 1, 0)) %>% 
 mutate_at(vars(contains("sig_")),
           ~ case_when(.x == 1 & sign == 1  ~ "Hot Spot",
                       .x == 1 & sign == -1 ~ "Cold Spot",
                       .x == 0              ~ "Neutral")) %>% 
 mutate_at(vars(contains("sig_")), ~factor(.x,
                                           levels = c("Cold Spot",
                                                      "Neutral",
                                                      "Hot Spot")))

# back to same length as original site obs. a wide DF
gistar_corrected_wide <- gistar_corrected1 %>% 
 dplyr::select(pnt_id, nb, sig_95, sig_99, sig_999) %>% 
 pivot_wider(names_from = nb,
             values_from = c(sig_95, sig_99, sig_999),
                  names_prefix = "nb_")

Gi_corrected_results <- left_join(gistar_corrected_wide, points_xy, by = "pnt_id") %>% 
  st_sf
```

``` r
 hot_cold_select <- function(.x){
    # cat(.x,"\n")
    # case of a single or uniquly NA string
    if(length(unique(.x)) == 1){
      if(is.na(unique(.x))){
        return("NA")
      }
    }
    # get rid of NA
    .y <- unique(as.character(na.omit(.x)))
    # if just Nuetral, pick it, otherwise Hot or Cold or Error if both
    # convert to numbers for rasterizing
    .y <- setdiff(.y,"Neutral")
    if(length(.y) == 0){
      .z <- "0"
    } else if(length(.y) == 1){
      if(.y == "Hot Spot"){
        .z <- "1" 
      } else if(.y == "Cold Spot"){
        .z <- "-1"
      } 
    } else if(length(.y) > 1){
      .z <- "ERROR"
    }
    return(.z)
 }


  ### aggregate hot/cold spot points to cell
  ### Pick hot or cold over nuetral, but "ERROR" if both hot and cold present
net_gi_star <- small_net %>% 
  st_sf() %>% 
  st_join(Gi_corrected_results) %>% 
  mutate_if(is.factor, as.character) %>% 
  group_by(net_id) %>% 
  summarise_at(.vars = vars(sig_95_nb_5:sig_999_nb_25), 
               list(~hot_cold_select(.))) %>% 
  mutate_at(.vars = vars(sig_95_nb_5:sig_999_nb_25), 
            list(~ifelse(. == "NA", NA, .))) %>% 
  mutate_at(.vars = vars(sig_95_nb_5:sig_999_nb_25), 
            list(~as.numeric(.)))

net_gi_star_NA_OMIT <- net_gi_star %>% 
  filter(!is.na(sig_95_nb_5))
  
```

``` r
 # ## Interpolation of z-score...
  # would need to loop across all z-score cols instead of summarise_at
  # in order to deal with nan in z-scores
  v <- dismo::voronoi(as(st_centroid(net_gi_star_NA_OMIT), "Spatial"))
#> Warning in st_centroid.sf(net_gi_star_NA_OMIT): st_centroid assumes attributes
#> are constant over geometries of x
  ca <- aggregate(as(st_buffer(survey_area,50), "Spatial"))
  vca <- raster::intersect(v, ca)
  # spplot(vca, 'sig_95_nb_5', col.regions=rev(get_col_regions()))
  rast_bounds  <- raster(as(st_buffer(survey_area,5), "Spatial"), res=6)
  
  # mapview(net_gi_star, zcol = "sig_95_nb_5") + 
  #   mapview(net_gi_star %>% filter(net_id == 1))
```

``` r
  ## aggregate mean value to all cells of just those that are NA
  ## the latter gives the data for tested cells more original data
  ## the former leaves a smoothed over and consistent method
  ## I think I like the latter. just rasterized values to NA cells.
  ## but need to do it for all different distance bands?
  
  sigs <- c(95,99,999)
  ## Reassign for safe keeping
  net_gi_star_new <- net_gi_star
  # list for results sigs * bands long goes here; name the list
  NA_cell_list <- vector(mode="list",length=length(sigs)*length(bands))
  #names here of in the loop
  list_names <- NULL
  loop_iter <- 1
  ## reminder: all of this here is just to fill in NA cells with no excvations
  for(i in seq_along(bands)){
    for(j in seq_along(sigs)){
      var_j <- paste0("sig_",sigs[j],"_nb_",bands[i])
      list_names <- c(list_names,var_j)
      cat("Looping on", var_j, ". Loop", loop_iter,"\n")
      net_gi_star_N_j <- net_gi_star_new %>%
        filter(is.na(get(var_j)))
      cat(nrow(net_gi_star_N_j),"\n")
      ## if >0 rows above, do raster stuff (slow'ish)
      if(nrow(net_gi_star_N_j)>0){
        cat("rasterize for",var_j, ". Loop", loop_iter,"\n")
        vr_j <- rasterize(vca, rast_bounds, var_j)
        # exract raster to the NA cells and 
        rv_j <- velox(vr_j)
        ## NEED TO RESCALE TO 1 or ZERO becuase of Mean
        mean_z <- as.numeric(rv_j$extract(net_gi_star_N_j, 
                                          fun = mean, small = TRUE))
        mean_z <- ifelse(mean_z > 0, 1, 
                         ifelse(mean_z < 0, -1, 0))
        # turn NA into zeros
        # NAs are b/c of no artefacts or NB size too large
        mean_z <- ifelse(is.na(mean_z), 0, mean_z)
        
        # update data
        net_gi_star_N_j[,var_j] <- mean_z
        net_gi_star_new <- net_gi_star_new %>% 
          left_join(., st_drop_geometry(net_gi_star_N_j[,c("net_id",var_j)]),
                    by = "net_id")
        # fill in NA spots for cells that did not have a hot/cold/neutral
        net_gi_star_new[[paste0(var_j,".x")]] <- ifelse(is.na(net_gi_star_new[[paste0(var_j,".x")]]),
                                                        net_gi_star_new[[paste0(var_j,".y")]],
                                                        net_gi_star_new[[paste0(var_j,".x")]])
        
        # add to results list
        NA_cell_list[[loop_iter]] <- mean_z
      } else {
        cat("Adding NULL for",var_j,". Loop", loop_iter,"\n")
        # add NULL to results list
        NA_cell_list[[loop_iter]] <- NULL
      }
      loop_iter <- loop_iter + 1
    } # end j
  }# end i
#> Looping on sig_95_nb_5 . Loop 1 
#> 45 
#> rasterize for sig_95_nb_5 . Loop 1 
#> Looping on sig_99_nb_5 . Loop 2 
#> 45 
#> rasterize for sig_99_nb_5 . Loop 2 
#> Looping on sig_999_nb_5 . Loop 3 
#> 45 
#> rasterize for sig_999_nb_5 . Loop 3 
#> Looping on sig_95_nb_10 . Loop 4 
#> 45 
#> rasterize for sig_95_nb_10 . Loop 4 
#> Looping on sig_99_nb_10 . Loop 5 
#> 45 
#> rasterize for sig_99_nb_10 . Loop 5 
#> Looping on sig_999_nb_10 . Loop 6 
#> 45 
#> rasterize for sig_999_nb_10 . Loop 6 
#> Looping on sig_95_nb_15 . Loop 7 
#> 45 
#> rasterize for sig_95_nb_15 . Loop 7 
#> Looping on sig_99_nb_15 . Loop 8 
#> 45 
#> rasterize for sig_99_nb_15 . Loop 8 
#> Looping on sig_999_nb_15 . Loop 9 
#> 45 
#> rasterize for sig_999_nb_15 . Loop 9 
#> Looping on sig_95_nb_20 . Loop 10 
#> 45 
#> rasterize for sig_95_nb_20 . Loop 10 
#> Looping on sig_99_nb_20 . Loop 11 
#> 45 
#> rasterize for sig_99_nb_20 . Loop 11 
#> Looping on sig_999_nb_20 . Loop 12 
#> 45 
#> rasterize for sig_999_nb_20 . Loop 12 
#> Looping on sig_95_nb_25 . Loop 13 
#> 45 
#> rasterize for sig_95_nb_25 . Loop 13 
#> Looping on sig_99_nb_25 . Loop 14 
#> 45 
#> rasterize for sig_99_nb_25 . Loop 14 
#> Looping on sig_999_nb_25 . Loop 15 
#> 45 
#> rasterize for sig_999_nb_25 . Loop 15
  names(NA_cell_list) <- list_names
  
  net_gi_star_new2 <- net_gi_star_new %>% 
    dplyr::select(contains(".x"), net_id) %>% 
    janitor::clean_names()
```
