
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
library(raster)
library(ggplot2)
library(dplyr)
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
points <- data.frame(x = rnorm(200, 50, 15),
                     y = rnorm(200, 50, 15)) %>% 
  filter(x <= 100 & y <= 100) %>% 
  filter(x >= 0 & y >= 0)

points$sim1 <- raster::extract(x = raster::rasterFromXYZ(g_pred),
                               y = points)

ggplot2::ggplot() +
  geom_raster(data = g_pred, 
              aes(x = x, y = y, fill = sim1),
              interpolate = FALSE) +
  geom_point(data = points, 
             aes(x = x, y = y),
             size = 4,
             color = "black") +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_void()
```

<img src="man/figures/README-create_sim_data-1.png" width="100%" />
