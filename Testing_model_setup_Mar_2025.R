
library(dplyr)
library(purrr)
library(here)
library(terra)
library(sf)
library(ecospat)
library(usdm)
library(ggplot2)
library(tidyterra)
library(dismo)
library(predicts)
library(blockCV)
library(scales)
library(mgcv)
library(randomForest)


# Source helper functions -------------------------------------------------

source(here("Scripts/Additional_functions.R"))

# Set group ---------------------------------------------------------------

group = "kelp_gulls"



# Set scenario ------------------------------------------------------------

#scenario = "dispersal_limited"
scenario = "enviro_limited"



# Set outpath -------------------------------------------------------------

outpath <- here("Outputs", group)


# Load biological records -------------------------------------------------

bio.records <- read.csv(here("Data/Biological_records/ICEFREE_GBIF_PTM_XY.csv"), header = T) %>% 
  filter(grepl(group, PTM_NAME)) 
  

# Load domain -------------------------------------------------------------


# Ice-Free Union layer sourced from: <http://dx.doi.org/doi:10.26179/7mnh-j215>

ice_free <- rast(here("Data/Environmental_predictors/ice_free_upsamp_1km.tif"))



# Assign bio presences to their ice-free grid cell xy---------------------


# Get index of ice-free cells that bio points intersect
iceID <- terra::extract(ice_free, 
                        bio.records[,c("x", "y")], 
                        xy = T)

# Bio records become the centre of the ice-free grid cell they fall within
bio.records <- iceID %>%
  dplyr::select(x, y) %>% 
  mutate(Presence = 1)


# Remove duplicates (records that happen to be in the same cell)
bio.records <- bio.records %>% 
  dplyr::distinct(x, y, Presence, .keep_all = T) 


# Load environmental covariates for target group --------------------------

group_covs <- source(here("Scripts/Defining_predictors_per_group.R"))$value

# Enviro. variable for functional group
cov_names <- group_covs[[group]]

# Load the variable names
env_var_loc <- source(here("Scripts/Defining_covariate_locations.R"))$value

# Call variables, stack and set names
covs <- env_var_loc[cov_names] %>% 
  map(~rast(.x)) %>% 
  rast() %>%  
  setNames(cov_names)


# Nearest-neighbour domain selection --------------------------------------

if(group %in% c("greater_sheathbill", "antarctic_shag", "antarctic_tern") & scenario == "dispersal_limited") {
  
  bio.records.sf <- st_as_sf(bio.records, coords = c("x", "y"), crs = "EPSG:3031")
  
  near <- st_nearest_feature(bio.records.sf)
  near <- bio.records.sf[near,]
  
  dist.near <- as.vector(st_distance(bio.records.sf, near, by_element = T))
  #dist.near <- quantile(dist.near, .95)
  dist.near <- max(dist.near)
  
  # Buffer records by the 95th percentile of record locations
  domain.mask <- st_buffer(bio.records.sf, dist.near) %>% 
    st_union() %>% 
    st_cast("POLYGON") 
  
  ice_free <- terra::mask(ice_free, vect(domain.mask)) %>% 
    crop(ext(vect(domain.mask)))
  
  covs <- terra::mask(covs, vect(domain.mask)) %>% 
    crop(ext(vect(domain.mask)))
  
  st_write(domain.mask, here("Data/Environmental_predictors/domain_mask_colobanthus.shp"))
  
  
}


# Dispersal limit domain selection ----------------------------------------

if(group %in% c("algae", "bankform_moss", "crust_lichens", "dry_loving_moss", "dry_nems", "Entomobryomorpha_springtails", "frufol_lichens", "leaf_liver", "Mesostigmata", "nrt", "Poduromorpha_springtails", "Prostigmata", "Sarcoptiformes", "wet_loving_moss") & scenario == "dispersal_limited") {
  
  
  # Get dispersal limit distance for that functional group
  
  source(here("Scripts/Defining_dispersal_limits_per_group.R"))
  
  
}


# Set domain to just peninsula or keep as continent-wide ---------------------

# 
# if(group %in% c("adelie_penguins", "chinstrap_penguins", "gentoo_penguins", "greater_sheathbill", "kelp_gulls", "skuas")) {
#   
#   # Load the Peninsula boundary
#   
#   peninsula <- st_read(here("Data/Environmental_predictors/Peninsula_Continent_Boundary.shp")) %>% 
#     vect(.)
#   
#   # Crop ice-free land and covariates to Peninsula
#   ice_free <- terra::mask(ice_free, peninsula)
#   covs <- terra::mask(covs, peninsula)
#     
#   ext <- c(-3334647.00583271,-1055542.13095899, 213869.492861746, 2592485.26286439)
#   
#   ice_free <- crop(ice_free, ext)
#   covs <- crop(covs, ext)
# 
# }

############################################
# Check variable correlation ----------------------------------------------
############################################

# 
# ## TESTS
# 
# cov_names <- c("dist_veg","dist_coast")
# 
# ecospat.cor.plot(covs)
# 
# usdm::vif(covs)
# 
# usdm::vifstep(covs)
# 
# corplots <- ENMTools::raster.cor.plot(envs_terra)
# corplots$cor.mds.plot
# 
# corplots$cor.heatmap
# 
# # PLOT CORRELATION OF VARIABLES OVER SPACE
# # Function to calculate correlation in a moving window
# focal_correlation <- function(x, method = "spearman") {
#   if (any(is.na(x))) return(NA)
#   cor(x[1:(length(x)/2)], x[(length(x)/2 + 1):length(x)], method = method)
# }
# 
# # Apply the moving window correlation
# cor_raster <- focal(covs, w = matrix(1, 3, 3), fun = function(x) focal_correlation(x, method="spearman"))
# 
# # Plot
# plot(cor_raster, main="Local Spearman Correlation")
# 


############################################
# Background sampling -----------------------------------------------------
############################################

# Count number of non-NA cells in the ice-free raster
global(ice_free, fun = "notNA")


# Set the location and number of background points

# Random sampling at first

background <- predicts::backgroundSample(mask = ice_free, 
                                         n = 50000,
                                         tryf = 100)

background <- vect(background, crs = "EPSG:3031")

# Convert background points (SpatVector) to data frame
background_df <- as.data.frame(geom(background))



# Plot presence and background points in ice-free layer -------------------

bio.recordsVect <- vect(bio.records, geom=c("x", "y"), crs = "EPSG:3031")

# plet(ice_free, "rock_union1", tiles = "", col = "gray40") %>%
#   points(background, col = "white", cex = 0.5, popup = T) %>%
#   points(bio.recordsVect, col = "red", cex = 0.5, popup = T)

 
# 
# # If I wanted to plot multiple variables
# test <- c(ice_free, covs)
# 
# plet(test, c(1, 3), tiles = "", collapse = F)


#writeVector(background, here(paste0("Data/Biological_records/", group, "/background_points.shp")), overwrite = T)


# Add presences and background --------------------------------------------


background_df <- background_df[,c("x", "y")] %>% 
  mutate(Presence = 0)
  
pr_bg <- rbind(bio.records, background_df)

# reset rownames for spatial tuning
rownames(pr_bg) <- NULL


# Extract enviro. covs for training ---------------------------------------

train_PB_covs <- terra::extract(covs, pr_bg[,c("x", "y")], xy = T)
train_PB_covs <- cbind(train_PB_covs, pr_bg["Presence"])

# Remove rows where there's values missing from at least one covariate

print(paste0("RECORDS FROM ", nrow(train_PB_covs) - sum(complete.cases(train_PB_covs)), " ROWS IN TRAINING DATA REMOVED DUE TO MISSING COVARIATE VALUES"))

train_PB_covs <- train_PB_covs[complete.cases(train_PB_covs), ] 
# Reset the row IDs to adjust for removed rows
rownames(train_PB_covs) <- NULL
train_PB_covs <- dplyr::select(train_PB_covs, -ID)

# Extract enviro. covs for *current* prediction ------------------------------

pred_cur_covs <- as.data.frame(covs, xy = T)

print(paste0("RECORDS FROM ", nrow(pred_cur_covs) - sum(complete.cases(pred_cur_covs)), " ROWS IN PREDICTION DATA REMOVED DUE TO MISSING COVARIATE VALUES"))

# Remove rows with NA in any covariates
pred_cur_covs <- pred_cur_covs[complete.cases(pred_cur_covs), ]
# Reset the row IDs to adjust for removed rows
rownames(pred_cur_covs) <- NULL


# ggplot() +
#   geom_spatraster(data = ice_free) + 
#   scale_fill_viridis_c(na.value = NA) +  
#   geom_point(data = bio.records, aes(x = x, y = y), color = "red", size = 0.5) +  
#   geom_point(data = bio.recordssave, aes(x = x, y = y), color = "blue", size = 0.5) +  
#   theme_minimal()


# Extract enviro. covs for *future* prediction ------------------------------

# pred_fut_covs <- as.data.frame(covs, xy = T)
# 
# print(paste0("RECORDS FROM ", nrow(pred_fut_covs) - sum(complete.cases(pred_fut_covs)), " ROWS IN PREDICTION DATA REMOVED DUE TO MISSING COVARIATE VALUES"))
# 
# # Remove rows with NA in any covariates
# pred_fut_covs <- pred_fut_covs[complete.cases(pred_fut_covs), ]
# # Reset the row IDs to adjust for removed rows
# rownames(pred_fut_covs) <- NULL
# pred_fut_covs <- dplyr::select(pred_fut_covs, -ID)



############################################
# Spatial cross-validation fold creation ----------------------------------
############################################

train_PB_sf <- st_as_sf(train_PB_covs[, c("x", "y", "Presence")], coords = c("x", "y"), crs = "EPSG:3031")

# Generate spatial blocks
spblock <- cv_spatial(x = train_PB_sf, 
                      column = "Presence",
                      r = NULL,
                      size = 100000, # Size of the blocks in metres
                      k = 5,
                      hexagon = TRUE,
                      selection = "random",
                      iteration = 100, # to find evenly-dispersed folds
                      biomod2 = FALSE)

cv_plot(cv = spblock,
        x = train_PB_sf,
        points_alpha = 0.5,
        nrow = 2)

## TO DO - HAVE SOME SORT OF PROTOCOL FOR COUNTING NUMBER OF RECORDS PER BLOCK

# # skip species with zero or low records
# if(min(spblock$records) < 1) next
# if(any(spblock$records$train_1 < 5)) next
# if(sum(spblock$records$test_1) < 10) next

# subset folds
spfolds <- spblock$folds_list


###########################################
# Run models for every fold of the spatial cross-validation ------------------
##########################################

# Set a single fold for now:
f = 1

for(f in seq_along(spfolds)) {
  

  # Subset the training and testing data (spatial cross validation) (for the fth fold)
  
  train_PB_covs_scv <- train_PB_covs[spfolds[[f]][[1]], ]
  test_PB_covs_scv <- train_PB_covs[spfolds[[f]][[2]], ]
  
  
  ############################################
  # Normalise covariates for some models  ---------------------------------
  ############################################
  
  # Models that don't require it: 
  # Lasso GLM
  # MaxEnt
  
  # Models that do require it:
  # BRT
  # RF - downsampled
  # GAM
  
  train_PB_covs_scv_norm <- train_PB_covs_scv
  test_PB_covs_scv_norm <- test_PB_covs_scv
  pred_cur_covs_norm <- pred_cur_covs
  #pred_fut_covs_norm <- pred_fut_covs
  
  
  for (v in cov_names) {
    
    meanv <- mean(train_PB_covs_scv[, v])
    sdv <- sd(train_PB_covs_scv[, v])
    train_PB_covs_scv_norm[, v] <- (train_PB_covs_scv[, v] - meanv) / sdv
    test_PB_covs_scv_norm[, v] <- (test_PB_covs_scv[, v] - meanv) / sdv
    pred_cur_covs_norm[, v] <- (pred_cur_covs[, v] - meanv) / sdv
    #pred_fut_covs_norm[, v] <- (pred_fut_covs[, v] - meanv) / sdv
    
  }
  
  
  ############################################
  # Calculating the case weights (down-weighting)  --------------------------
  ############################################
  
  prNum <- as.numeric(table(train_PB_covs_scv$Presence)["1"]) # number of presences
  bgNum <- as.numeric(table(train_PB_covs_scv$Presence)["0"]) # number of backgrounds
  wt <- ifelse(train_PB_covs_scv$Presence == 1, 1, prNum / bgNum) # down-weighting
  
  
  ############################################
  # Maxent default ----------------------------------------------------------
  ############################################
  
  #*Downloaded 17th March 2025*
  
  # Set model type
  mod.type = "Maxent"
  
    # Set type of covariates (normalised or not)
  
  if(mod.type %in% c("GAM", "BRT", "RF")) {
    train_PB_covs_scv <- train_PB_covs_scv_norm
    test_PB_covs_scv <- test_PB_covs_scv_norm
    pred_cur_covs <- pred_cur_covs_norm
    pred_fut_covs <- pred_fut_covs_norm
  }
  
  # Make folder:
  
  dir.create(file.path(outpath, "Maxent_outputs"), showWarnings = F)
  
  # x is a dataframe with the covariates, each row is a PB point
  # p is the occurrences, 0s or 1s for each PB point
  
  maxent.mod <- NULL
  maxent.mod <- dismo::maxent(x = train_PB_covs_scv[, cov_names],
                              p = train_PB_covs_scv[["Presence"]],
                              removeDuplicates = FALSE,
                              path = file.path(outpath, "Maxent_outputs"),
                              args = c("nothreshold"))
  
  # Testing set prediction
  pred_test.mxt <- dismo::predict(maxent.mod, test_PB_covs_scv, args = "doclamp=false")
  pred_test.mxt <- cbind(test_PB_covs_scv, pred_test.mxt)
  colnames(pred_test.mxt)[grepl("pred", colnames(pred_test.mxt))] <- "pred"
  
  # Evaluate prediction on test set
  eval <- evaluate_prediction(pred_test.mxt)
  
  eval_df.mxt <- data.frame(group = group, 
                            fold = f, 
                            model = "MaxEnt", 
                            eval)
  
  # Save fold outputs and evaluation from prediction to test dataset
  write.csv(pred_test.mxt, file = file.path(outpath, "Maxent_outputs", paste0("Prediction_TEST_fold_", f, ".csv")))
  write.csv(eval_df.mxt, file = file.path(outpath, "Maxent_outputs", paste0("Evaluation_metric_fold_", f, ".csv")))
  
  
  # Current prediction
  pred_cur.mxt <- dismo::predict(maxent.mod, pred_cur_covs, args = "doclamp=false")
  pred_cur.mxt <- cbind(pred_cur_covs, pred_cur.mxt)
  colnames(pred_cur.mxt)[grepl("pred", colnames(pred_cur.mxt))] <- "pred"
  # pred_cur.rast <- rast(pred_cur[, c("x", "y", "pred")], type = "xyz", crs = "EPSG:3031")
  
  # Future prediction
  # pred_fut.mxt <- dismo::predict(maxent.mod, pred_fut_covs, args = "doclamp=false")
  
  # plet(pred.rast, "pred", tiles = "")
  # 
  # plet(pred.rast, "pred", tiles = "") %>% 
  #   points(background, col = "white", cex = 0.5, popup = T) %>% 
  #   points(bio.recordsVect, col = "red", cex = 0.5, popup = T)
  # 
  
  
  ############################################
  # Maxent tuned with spatial block cross validation ---------------------------
  ############################################
  
  # Set type of model
  mod.type = "Maxent-tuned"
  
  # Set type of covariates (normalised or not)
  
  if(mod.type %in% c("GAM", "BRT", "RF")) {
    train_PB_covs_scv <- train_PB_covs_scv_norm
    test_PB_covs_scv <- test_PB_covs_scv_norm
    pred_cur_covs <- pred_cur_covs_norm
    #pred_fut_covs <- pred_fut_covs_norm
  }
  
  # Make folder:
  
  dir.create(file.path(outpath, "Maxent_tuned_outputs"), showWarnings = F)
  
  
  # Setup for saving MaxEnt tuned parameters
  tune_df <- data.frame(group = NA, fold = NA, model = NA, reg = NA, args = NA)
  
  n <- 0
  z <- 0
  
  
  # This just pulls out the *testing* data for all the folds except fold f
  # This is the same as our training data for fold f, except that it's pre-partitioned into our 4 folds
  
  myspfolds <- foldsID(fld = spfolds)
  myspfolds <- myspfolds[-f]
  
  library(precrec)
  
  param_optim <- NULL
  param_optim <- maxent_param(data = train_PB_covs, # training data from all folds (for indexing)
                              folds = myspfolds, # testing data from the 4 other folds
                              filepath = file.path(outpath, "Maxent_tuned_outputs"))
  
  
  # save the tuned parameters
  tune_df[z <- z + 1, ] <- c(group, f, "MaxEnt-tuned", param_optim[1], param_to_txt(param_optim))
  
  # x is a dataframe with the covariates, each row is a PB point
  # p is the occurrences, 0s or 1s for each PB point
  
  maxent.tune.mod <- NULL
  maxent.tune.mod <- dismo::maxent(x = train_PB_covs_scv[, cov_names],
                                   p = train_PB_covs_scv[["Presence"]],
                                   removeDuplicates = FALSE,
                                   path = file.path(outpath, "Maxent_tuned_outputs"),
                                   args = param_optim)
  
  # Testing set prediction
  pred_test.mxttune <- dismo::predict(maxent.tune.mod, test_PB_covs_scv, args = "doclamp=false")
  pred_test.mxttune <- cbind(test_PB_covs_scv, pred_test.mxttune)
  colnames(pred_test.mxttune)[grepl("pred", colnames(pred_test.mxttune))] <- "pred"
  
  # Evaluate prediction on test set
  eval <- evaluate_prediction(pred_test.mxttune)
  
  eval_df.mxttune <- data.frame(group = group, 
                            fold = f, 
                            model = "MaxEnt-tuned", 
                            eval)
  
  # Save fold outputs and evaluation from prediction to test dataset
  write.csv(pred_test.mxttune, file = file.path(outpath, "Maxent_tuned_outputs", paste0("Prediction_TEST_fold_", f, ".csv")))
  write.csv(pred_test.mxttune, file = file.path(outpath, "Maxent_tuned_outputs", paste0("Evaluation_metric_fold_", f, ".csv")))
  
  
  # Current prediction
  pred_cur.mxttune <- dismo::predict(maxent.mod, pred_cur_covs, args = "doclamp=false")
  pred_cur.mxttune <- cbind(pred_cur_covs, pred_cur.mxttune)
  colnames(pred_cur.mxttune)[grepl("pred", colnames(pred_cur.mxttune))] <- "pred"
  # pred_cur.rast <- rast(pred_cur[, c("x", "y", "pred")], type = "xyz", crs = "EPSG:3031")
  
  # Future prediction
  # pred_fut.mxttune <- dismo::predict(maxent.mod, pred_fut_covs, args = "doclamp=false")
  
  # plet(pred.rast, "pred", tiles = "")
  # 
  # plet(pred.rast, "pred", tiles = "") %>% 
  #   points(background, col = "white", cex = 0.5, popup = T) %>% 
  #   points(bio.recordsVect, col = "red", cex = 0.5, popup = T)
  
  
  
  ############################################
  # LASSO regression --------------------------------------------------------
  ############################################
  
  # Set type of model
  mod.type = "LASSO"
  
  # Set type of covariates (normalised or not)
  
  if(mod.type %in% c("GAM", "BRT", "RF")) {
    train_PB_covs_scv <- train_PB_covs_scv_norm
    test_PB_covs_scv <- test_PB_covs_scv_norm
    pred_cur_covs <- pred_cur_covs_norm
    #pred_fut_covs <- pred_fut_covs_norm
  }
  
  # Make folder:
  
  dir.create(file.path(outpath, "LASSO_outputs"), showWarnings = F)
  
  # Convert input data to a matrix
  # Make quadratic terms for glmnet package
  
  # To make the orthogonal quadratic features for glmnet package (see more detail in the supplementary material
  # of Guillera-Arroita et al. 2014), the make_quadratic function in myspatial package is used. The object
  # generated by this function can be used in the generic predict() function to apply the transformation on the
  # training and testing datasets and later used in predicting to rasters. The package myspatial is archived in
  # GitHub and can be installed using the following code.
  
  # installing the package from github
  # remotes::install_github("rvalavi/myspatial")
  
  library(glmnet)
  library(myspatial)
  
  quad_obj <- make_quadratic(train_PB_covs_scv, cols = cov_names)
  
  # now we can predict this quadratic object on the training and testing data
  # this make two columns for each covariates used in the transformation
  train_quad <- predict(quad_obj, newdata = train_PB_covs_scv)
  test_quad <- predict(quad_obj, newdata = test_PB_covs_scv)
  pred_cur_quad <- predict(quad_obj, newdata = pred_cur_covs)
  # pred_fut_quad <- predict(quad_obj, newdata = pred_fut_covs)
  
  # convert the data.frames to sparse matrices
  # select all quadratic (and non-quadratic) columns, except the y (occ)
  
  new_vars <- names(train_quad)[names(train_quad) != c("Presence")]
  train_sparse <- sparse.model.matrix(~. -1, train_quad[, new_vars])
  test_sparse <- sparse.model.matrix(~. -1, test_quad[, new_vars])
  pred_cur_sparse <- sparse.model.matrix(~. -1, pred_cur_quad[, new_vars])
  # pred_fut_sparse <- sparse.model.matrix(~. -1, pred_fut_quad[, new_vars])
  
   # Fitting the lasso GLM
  
  lasso <- glmnet(x = training_sparse,
                  y = training_quad$Presence,
                  family = "binomial",
                  alpha = 1, # here 1 means fitting lasso
                  weights = wt)
  
  plot(lasso, xvar = "lambda", label = TRUE)
  
  lasso_cv <- cv.glmnet(x = training_sparse,
                        y = training_quad$Presence,
                        family = "binomial",
                        alpha = 1, # here 1 means fitting lasso
                        weights = wt,
                        nfolds = 10) # number of folds for cross-validation
  
  plot(lasso_cv)
  
  # Choose to use the lambda (regularisation parameter) that is 1SD from the minimum deviance (to avoid over-fitting if use the minimum deviance lambda)
  
  pred_test.lasso <- predict(lasso_cv, test_sparse, s = "lambda.1se", type = "response")
  pred_test.lasso <- cbind(test_PB_covs_scv, pred_test.lasso)
  colnames(pred_test.lasso)[grepl("lambda.1se", colnames(pred_test.lasso))] <- "pred"
  
  # Evaluate prediction on test set
  eval <- evaluate_prediction(pred_test.lasso)
  
  eval_df.lasso <- data.frame(group = group,
                              fold = f,
                              model = "LASSO",
                              eval)
  
  # Save fold outputs and evaluation from prediction to test dataset
  write.csv(pred_test.lasso, file = file.path(outpath, "LASSO_outputs", paste0("Prediction_TEST_fold_", f, ".csv")))
  write.csv(eval_df.lasso, file = file.path(outpath, "LASSO_outputs", paste0("Evaluation_metric_fold_", f, ".csv")))
  
  # Current prediction
  pred_cur.lasso <- predict(lasso_cv, pred_cur_sparse, s = "lambda.1se", type = "response")
  pred_cur.lasso <- cbind(pred_cur_covs, pred_cur.lasso)
  colnames(pred_cur.lasso)[grepl("lambda.1se", colnames(pred_cur.lasso))] <- "pred"
  
  # Future prediction
  # pred_fut.lasso <- predict(lasso_cv, pred_fut_sparse, s = "lambda.1se", type = "response")
  # pred_fut.lasso <- cbind(pred_fut_covs, pred_fut.lasso)
  # colnames(pred_fut.lasso)[grepl("lambda.1se", colnames(pred_fut.lasso))] <- "pred"
  
  
  
  ############################################
  # GAM ---------------------------------------------------------------------
  ############################################
  
  # Set type of model
  mod.type = "GAM"
  
  # Set type of covariates (normalised or not)
  
  if(mod.type %in% c("GAM", "BRT", "RF")) {
    train_PB_covs_scv <- train_PB_covs_scv_norm
    test_PB_covs_scv <- test_PB_covs_scv_norm
    pred_cur_covs <- pred_cur_covs_norm
    #pred_fut_covs <- pred_fut_covs_norm
  }
  
  # Make folder:
  
  dir.create(file.path(outpath, "GAM_outputs"), showWarnings = F)
  
  # Building GAM formula
  myform <- paste("occ ~", paste(paste0("s(", cov_names, ")"), collapse = "+"))
  
  # Rename "Presence" to "occ"
  train_PB_covs_scv.gam <- train_PB_covs_scv %>% 
    rename(occ = Presence)
  
  gam <- mgcv::gam(formula = as.formula(myform), 
                       data = train_PB_covs_scv.gam,
                       family = binomial(link = "logit"),
                       weights = wt,
                       method = "REML")

  pred_test.gam <- predict(gam, test_PB_covs_scv, type = "response")
  pred_test.gam <- cbind(test_PB_covs_scv, pred_test.gam)
  colnames(pred_test.gam)[grepl("pred", colnames(pred_test.gam))] <- "pred"
  
  # Evaluate prediction on test set
  eval <- evaluate_prediction(pred_test.gam)
  
  eval_df.gam <- data.frame(group = group,
                              fold = f,
                              model = "GAM",
                              eval)
  
  # Save fold outputs and evaluation from prediction to test dataset
  write.csv(pred_test.gam, file = file.path(outpath, "GAM_outputs", paste0("Prediction_TEST_fold_", f, ".csv")))
  write.csv(eval_df.gam, file = file.path(outpath, "GAM_outputs", paste0("Evaluation_metric_fold_", f, ".csv")))
 
  # Current prediction
  pred_cur.gam <- predict(gam, pred_cur_covs, type = "response")
  pred_cur.gam <- cbind(pred_cur_covs, pred_cur.gam)
  colnames(pred_cur.gam)[grepl("pred", colnames(pred_cur.gam))] <- "pred"
  
  # Future prediction
  # pred_fut.gam <- predict(gam, pred_fut_covs, type = "response")
  # pred_fut.gam <- cbind(pred_fut_covs, pred_fut.gam)
  # colnames(pred_fut.gam)[grepl("pred", colnames(pred_fut.gam))] <- "pred"
  

  ############################################
  # Boosted regression tree -------------------------------------------------
  ############################################
 
  # Set type of model
  mod.type = "BRT"
  
  # Set type of covariates (normalised or not)
  
  if(mod.type %in% c("GAM", "BRT", "RF")) {
    train_PB_covs_scv <- train_PB_covs_scv_norm
    test_PB_covs_scv <- test_PB_covs_scv_norm
    pred_cur_covs <- pred_cur_covs_norm
    #pred_fut_covs <- pred_fut_covs_norm
  }
  
  
  # Make folder:
  
  dir.create(file.path(outpath, "BRT_outputs"), showWarnings = F)
  
  # Initialise the tracker
  brt <- NULL
  b <- 0 
  
  # Attempt the model up to 30 times, with adjustments to parameterisation. 
  while(is.null(brt)){
    b <- b + 1
    if(b < 11){
      ntrees <- 50
      lrate <- 0.001
    } else if(b < 21){
      lrate <- 0.0001
    } else if(b < 31){
      ntrees <- 25
      lrate <- 0.0001
    } else{
      break # After 30 tries, give up and move to next model
    }
    
    brt <- dismo::gbm.step(data = train_PB_covs_scv,
                           gbm.x = which(colnames(train_PB_covs_scv) %in% cov_names), # Column indices for covs
                           gbm.y = which(colnames(train_PB_covs_scv) == "Presence"), # Column index for Presence
                           family = "bernoulli",
                           tree.complexity = ifelse(prNum < 50, 1, 5),
                           learning.rate = lrate,
                           bag.fraction = 0.75,
                           max.trees = 10000,
                           n.trees = ntrees,
                           n.folds = 5, # 5-fold cross-validation
                           site.weights = wt,
                           silent = TRUE) # Avoid printing cv results
    
  }
  
  if(is.null(brt)) {
    
    print(paste0("BRT model failed to converge"))
    
    
  } else { # If model converged:
    
    # The optimal number of trees (estimated from k-fold CV)
    brt$gbm.call$best.trees
    
    # Predict with best number of trees
    pred_test.brt <- dismo::predict(brt, test_PB_covs_scv, n.trees = brt$gbm.call$best.trees, type = "response")
    pred_test.brt <- cbind(test_PB_covs_scv, pred_test.brt)
    colnames(pred_test.brt)[grepl("pred", colnames(pred_test.brt))] <- "pred"
    
    # Evaluate prediction on test set
    eval <- evaluate_prediction(pred_test.brt)
    
    eval_df.brt <- data.frame(group = group,
                              fold = f,
                              model = "BRT",
                              eval)
    
    # Save fold outputs and evaluation from prediction to test dataset
    write.csv(pred_test.brt, file = file.path(outpath, "BRT_outputs", paste0("Prediction_TEST_fold_", f, ".csv")))
    write.csv(eval_df.brt, file = file.path(outpath, "BRT_outputs", paste0("Evaluation_metric_fold_", f, ".csv")))
    
    # Current prediction
    pred_cur.brt <- dismo::predict(brt, pred_cur_covs, n.trees = brt$gbm.call$best.trees, type = "response")
    pred_cur.brt <- cbind(pred_cur_covs, pred_cur.brt)
    colnames(pred_cur.brt)[grepl("pred", colnames(pred_cur.brt))] <- "pred"
    # pred_cur.rast <- rast(pred_cur[, c("x", "y", "pred")], type = "xyz", crs = "EPSG:3031")
    
    # Future prediction
    # pred_fut.brt <- dismo::predict(brt, pred_fut_covs, n.trees = brt$gbm.call$best.trees, type = "response")
    # pred_fut.brt <- cbind(pred_fut_covs, pred_fut.brt)
    # colnames(pred_fut.brt)[grepl("pred", colnames(pred_fut.brt))] <- "pred"
    
    
  }
  
  
  ############################################
  # Random forest -----------------------------------------------------------
  ############################################
  
  # Set type of model
  mod.type = "RF"
  
  # Set type of covariates (normalised or not)
  
  if(mod.type %in% c("GAM", "BRT", "RF")) {
    train_PB_covs_scv <- train_PB_covs_scv_norm
    test_PB_covs_scv <- test_PB_covs_scv_norm
    pred_cur_covs <- pred_cur_covs_norm
    #pred_fut_covs <- pred_fut_covs_norm
  }
  
  # Make folder:
  
  dir.create(file.path(outpath, "RF_outputs"), showWarnings = F)
  
  # Convert the response to factor for producing class relative likelihood
  train_PB_covs_scv$Presence <- as.factor(train_PB_covs_scv$Presence)
  
  # For down-sampling, the number of background (0s) in each bootstrap sample should the same as presences
  # (1s). For this, we use sampsize argument to do this.
  smpsize <- c("0" = prNum, "1" = prNum)
  
  rf <- randomForest::randomForest(formula = Presence ~.,
                                   data = train_PB_covs_scv,
                                   ntree = 1000,
                                   sampsize = smpsize,
                                   replace = T)
                     
  plot(rf, main = "RF down-sampled")
  
  pred_test.rf <- predict(rf, test_PB_covs_scv, type = "prob")[, "1"]
  pred_test.rf <- cbind(test_PB_covs_scv, pred_test.rf)
  colnames(pred_test.rf)[grepl("pred", colnames(pred_test.rf))] <- "pred"
  
  # Evaluate prediction on test set
  eval <- evaluate_prediction(pred_test.rf)
  
  eval_df.rf <- data.frame(group = group,
                              fold = f,
                              model = "RF",
                              eval)
  
  # Save fold outputs and evaluation from prediction to test dataset
  write.csv(pred_test.rf, file = file.path(outpath, "RF_outputs", paste0("Prediction_TEST_fold_", f, ".csv")))
  write.csv(eval_df.rf, file = file.path(outpath, "RF_outputs", paste0("Evaluation_metric_fold_", f, ".csv")))
  
  # Current prediction
  pred_cur.rf <- predict(rf, pred_cur_covs, type = "prob")[, "1"]
  pred_cur.rf <- cbind(pred_cur_covs, pred_cur.rf)
  colnames(pred_cur.rf)[grepl("pred", colnames(pred_cur.rf))] <- "pred"
  
  # Future prediction
  # pred_fut.rf <- predict(rf, pred_fut_covs, type = "prob")[, "1"]
  # pred_fut.rf <- cbind(pred_fut_covs, pred_fut.rf)
  # colnames(pred_fut.rf)[grepl("pred", colnames(pred_fut.rf))] <- "pred"
  
  
  ############################################
  # Ensemble of all model predictions ---------------------------------------
  ############################################
  
  pred_cur_ens <- data.frame(maxent = scales::rescale(pred_cur.mxt[["pred"]], to = c(0,1)),
                             maxent.tune = scales::rescale(pred_cur.mxttune[["pred"]], to = c(0,1)),
                             lasso = scales::rescale(pred_cur.lasso[["pred"]], to = c(0,1)),
                             gam = scales::rescale(pred_cur.gam[["pred"]], to = c(0,1)),
                             brt = scales::rescale(pred_cur.brt[["pred"]], to = c(0,1)),
                             rf = scales::rescale(pred_cur.rf[["pred"]], to = c(0,1))) %>% 
    rowMeans()
  
  pred_fut_ens <- data.frame(maxent = scales::rescale(pred_fut.mxt[["pred"]], to = c(0,1)),
                             maxent.tune = scales::rescale(pred_fut.mxttune[["pred"]], to = c(0,1)),
                             lasso = scales::rescale(pred_fut.lasso[["pred"]], to = c(0,1)),
                             gam = scales::rescale(pred_fut.gam[["pred"]], to = c(0,1)),
                             brt = scales::rescale(pred_fut.brt[["pred"]], to = c(0,1)),
                             rf = scales::rescale(pred_fut.rf[["pred"]], to = c(0,1))) %>% 
    rowMeans()
  
  write.csv(pred_cur_ens, file = here("Outputs", group, paste0("Prediction_ensemble_current_fold_", f, ".csv")))
  write.csv(pred_fut_ens, file = here("Outputs", group, paste0("Prediction_ensemble_future_fold_", f, ".csv")))
  
}


############################################
# Fitting models to complete dataset ---------------------------------------
############################################




# ARCHIVE -----------------------------------------------------------------

# # Normalise covariates and **save mean and sd for normalising prediction covs**
# 
# mean <- global(covs, fun = "mean", na.rm = T)[,1]
# 
# sd <- global(covs, fun = "sd", na.rm = T)[,1]
# 
# covs_norm <- (covs - mean) / sd


