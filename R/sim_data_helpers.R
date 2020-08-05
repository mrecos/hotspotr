
#' rescale_sim_raster
#'
#' Helper function to rescale a simulated landscape to match the mean and standard deviation of the simulated data used to fit a model.
#' `rescale_sim_raster()` - Function that rescales simulated rasters from `NLMR::nlm_gaussianfield` (Sciaini, Fritsch, and Simpkins 2017) or whatever you want to use, to the mean and standard deviation of the simulated data used to fit the klr model. You will have to add the mean and sd arguments manually based on what you put into the get_sim_data function. The example in the code above inputs the default mean and sd values from the defualts of the `get_sim_data()` function. Returned is a raster scaled too your simualted training data.
#'
#' @param rast [raster] input raster object
#' @param mean [numeric] value of mean used in `get_sim_data()` function.
#' @param sd   [numeric] value of standard deviation used in `get_sim_data()` function.
#'
#' @return [raster] a raster object scaled to the input data
#' @export
#'
#' @examples
#' \dontrun{
#' ### width and hieght of roving focal window (required)
#' ngb = 5
#' ### Number of rows and columns in prediction rasters
#' ## needed for making simulated rasters, as well as for predicting real-world rasters
#' cols = 100
#' rows = 100
#'
#' ### Create simulated environmental rasters  (sim data only) ####
#' s_var1r <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = 20)
#' s_var1 <- rescale_sim_raster(s_var1r, 50, 10)
#' s_var2 <- rescale_sim_raster(s_var1r, 3, 2)
#' b_var1r <- NLMR::nlm_gaussianfield(cols,rows,autocorr_range = 20)
#' b_var1 <- rescale_sim_raster(b_var1r, 100, 20)
#' b_var2 <- rescale_sim_raster(b_var1r, 6, 3)
#' ### Create a site-present trend surface  (sim data only)
#' trend_coords <- sim_trend(cols, rows, n = 3)
#' coords <- trend_coords$coords
#' trend <- trend_coords$trend
#' inv_trend <- abs(1-trend)
#' var1 <- (s_var1 * trend) + (b_var1 * inv_trend)
#' var2 <- (s_var2 * trend) + (b_var2 * inv_trend)
#' #### end simulated data creation ####
#'
#' ### Create raster stack of predictor variables
#' pred_var_stack <- raster::stack(var1, var2)
#' names(pred_var_stack) <- c("var1","var2")
#' ### scale rasters to training data
#' pred_var_stack_scaled <- scale_prediction_rasters(pred_var_stack, params, verbose = 0)
#' ### Predict raster (single chunk, not in parallel)
#' pred_rast <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params, split = FALSE, ppside = NULL,
#'                                 progress = FALSE, parallel = FALSE)
#' }
#'
rescale_sim_raster <- function(rast,mean,sd){
  r_vals <- rast@data@values
  rast[] <- mean + (r_vals - mean(r_vals)) * (sd/sd(r_vals))
  return(rast)
}


#' sim_trend
#'
#' A helper function to create a trend surface that when applied to a simulated landscape raster correlates it with simulated site locations
#'
#' `sim_trend()` - Function is used to take create `n` number of simulated site locations of `size` cell dimensions on a rows by cols raster. The latter two arguments should match the size of your simulated rasters. The function randomly locates the sites and then creates a distance gradient (trend) from the site locations outward. The trend is a value 1 at the sites and reduces to 0 at the maximum combined distance from all sites. The output of this function is a list of a matrix of simulated site x/y coordinates (centers) and a raster of the trend surface. The point of the trend is to then combine it with the simulated rasters (as down in the code above) such that the raster depicting site-likely conditions is multiplied by the trend to output a raster retaining site-likely conditions near simulated site locations. Conversely, the site-unlikely simulated raster is multiplied by the inverse of the trend to result in a raster retaining site-unlikely characteristics away from the site locations. When those two rasters are added you get a simulated environment that is more preferable to site locations near site locations. It is a bit involved for something that had nothing to do with the actual KLRfome model, but it is needed to produce actual correlated environments for model testing.
#'
#' @param cols [integer] the number of columns in the simulated landscape raster
#' @param rows [integer] the number of rows in the simulated landscape raster
#' @param n    [integer] the number of columns in the simulated landscape raster
#' @param size [integer] the pixel dimensions for simualted sites
#'
#' @importFrom NLMR nlm_distancegradient
#' @return [list] a list with two elements, `coords` are the center coordiantes of the simulated sites, and `trend` is the raster trend surface to be applied to the simulated landscape.
#' @export
#'
#' @examples
#' \dontrun{
#' ### width and hieght of roving focal window (required)
#' ngb = 5
#' ### Number of rows and columns in prediction rasters
#' ## needed for making simulated rasters, as well as for predicting real-world rasters
#' cols = 100
#' rows = 100
#'
#' ### Create simulated environmental rasters  (sim data only) ####
#' s_var1r <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = 20)
#' s_var1 <- rescale_sim_raster(s_var1r, 50, 10)
#' s_var2 <- rescale_sim_raster(s_var1r, 3, 2)
#' b_var1r <- NLMR::nlm_gaussianfield(cols,rows,autocorr_range = 20)
#' b_var1 <- rescale_sim_raster(b_var1r, 100, 20)
#' b_var2 <- rescale_sim_raster(b_var1r, 6, 3)
#' ### Create a site-present trend surface  (sim data only)
#' trend_coords <- sim_trend(cols, rows, n = 3)
#' coords <- trend_coords$coords
#' trend <- trend_coords$trend
#' inv_trend <- abs(1-trend)
#' var1 <- (s_var1 * trend) + (b_var1 * inv_trend)
#' var2 <- (s_var2 * trend) + (b_var2 * inv_trend)
#' #### end simulated data creation ####
#'
#' ### Create raster stack of predictor variables
#' pred_var_stack <- raster::stack(var1, var2)
#' names(pred_var_stack) <- c("var1","var2")
#' ### scale rasters to training data
#' pred_var_stack_scaled <- scale_prediction_rasters(pred_var_stack, params, verbose = 0)
#' ### Predict raster (single chunk, not in parallel)
#' pred_rast <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params, split = FALSE, ppside = NULL,
#'                                 progress = FALSE, parallel = FALSE)
#' }
#'
sim_trend <- function(cols, rows, n = 2, size = 2){
  trend <- raster(ncols = cols,nrows = rows,res = 1,
                  xmn=0, xmx = cols, ymn=0, ymx=rows, vals = 0)
  coords <- matrix(ncol=2,nrow=n)
  colnames(coords) <- c("x","y")
  for(i in seq_len(n)){
    crd <- c(sample(1:(rows-size),1),sample(1:(cols-size),1))
    coords[i,1] <- crd[1]
    coords[i,2] <- crd[2]
    # had to adjust origin to match expected from crds.
    site_i <- NLMR::nlm_distancegradient(cols,rows,
                                         origin = c(rows-crd[2],rows-crd[2]+size,
                                                    crd[1],crd[1]-size))^0.5
    trend <- trend + site_i
  }
  mn <- raster::cellStats(trend, min)
  mx <- raster::cellStats(trend, max)
  trend <- abs(1-((trend-mn)/(mx-mn))) # standardize and invert
  return(list(trend = trend, coords = coords))
}

#' scale_prediction_rasters
#'
#' Center and scale the prediction raster stack to the parameters of the training data used to fit a model
#'
#'  `scale_prediction_rasters()` - Function scales your predictor rater stack based on the params list created in the model fitting process. This script simply loops over the rasters in the stack and centers and scales based on mean and sd of the training data used to fit the klr model. The function outputs a raster stack.
#'
#' @param pred_var_stack [raster stack] stack of input landscape rasters for prediction
#' @param params [list] list of mean and standard deviation parameters returned from formatting data. (See example below)
#' @param verbose [integer] a value of `1` to print out the landscape variable name as it is being normalized.
#'
#' @importFrom raster scale
#' @return [raster stack] a stack of the prediction rasters scaled to the the input training parameters.
#' @export
#'
#' @examples
#' \dontrun{
#' ### Create param list
#' params <- list(train_data = train_data,
#' alphas_pred = train_log_pred[["alphas"]],
#' sigma = sigma,
#' lambda = lambda,
#' means = formatted_data$means,
#' sds = formatted_data$sds)
#'
#' ### width and hieght of roving focal window (required)
#' ngb = 5
#' ### Number of rows and columns in prediction rasters
#' ## needed for making simulated rasters, as well as for predicting real-world rasters
#' cols = 100
#' rows = 100
#'
#' ### Create simulated environmental rasters  (sim data only) ####
#' s_var1r <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = 20)
#' s_var1 <- rescale_sim_raster(s_var1r, 50, 10)
#' s_var2 <- rescale_sim_raster(s_var1r, 3, 2)
#' b_var1r <- NLMR::nlm_gaussianfield(cols,rows,autocorr_range = 20)
#' b_var1 <- rescale_sim_raster(b_var1r, 100, 20)
#' b_var2 <- rescale_sim_raster(b_var1r, 6, 3)
#' ### Create a site-present trend surface  (sim data only)
#' trend_coords <- sim_trend(cols, rows, n = 3)
#' coords <- trend_coords$coords
#' trend <- trend_coords$trend
#' inv_trend <- abs(1-trend)
#' var1 <- (s_var1 * trend) + (b_var1 * inv_trend)
#' var2 <- (s_var2 * trend) + (b_var2 * inv_trend)
#' #### end simulated data creation ####
#'
#' ### Create raster stack of predictor variables
#' pred_var_stack <- raster::stack(var1, var2)
#' names(pred_var_stack) <- c("var1","var2")
#' ### scale rasters to training data
#' pred_var_stack_scaled <- scale_prediction_rasters(pred_var_stack, params, verbose = 0)
#' ### Predict raster (single chunk, not in parallel)
#' pred_rast <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params, split = FALSE, ppside = NULL,
#'                                 progress = FALSE, parallel = FALSE)
#' }
#'
scale_prediction_rasters <- function(pred_var_stack, params, verbose = 1){
  pred_var_stack_scaled <- pred_var_stack
  for(i in seq_len(dim(pred_var_stack)[3])){
    var_name <- names(pred_var_stack[[i]])
    if(verbose == 1){
      cat("Normalizing:", var_name, "\n")
    }
    pred_var_stack_scaled[[var_name]] <- raster::scale(pred_var_stack[[var_name]],
                                                       center = params$means[var_name],
                                                       scale  = params$sds[var_name])
  }
  return(pred_var_stack_scaled)
}
