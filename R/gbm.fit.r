#' @export gbm.fit
gbm.fit <- function(x,y,
                    offset = NULL,
                    distribution = "Bernoulli",
                    w = NULL,
                    var.monotone = NULL,
                    n.trees = 100,
                    interaction.depth = 1,
                    n.minobsinnode = 10,
                    shrinkage = 0.001,
                    bag.fraction = 0.5,
                    nTrain = NULL,
                    train.fraction = NULL,
                    mFeatures = NULL,
                    keep.data = TRUE,
                    verbose = TRUE,
                    var.names = NULL,
                    response.name = "y",
                    group = NULL,
                    prior.node.coeff.var = 1000,
                    strata = NULL, obs.id = 1:nrow(x)) {
  # Highlight new API
  warning("gbm is depracated - using gbm2...")
  
  # Reconstruct data
  data <- data.frame(y, x)
  if(is.null(var.names)) {
    var.names <- getVarNames(x)
  }
  if(length(response.name) != ncol(y)) {
      response.name <- paste("X", seq_len(ncol(y)) ,sep="")
  }
  colnames(data) <- c(response.name, var.names)
  
  # Calculate nTrain if necessary 
  if(!is.null(nTrain) && !is.null(train.fraction)) {
    stop("Parameters 'nTrain' and 'train.fraction' cannot both be specified")
    
  } else if(!is.null(train.fraction)) {
    warning("Parameter 'train.fraction' of gbm.fit is deprecated, please specify 'nTrain' instead")
    nTrain <- floor(train.fraction*length(unique(obs.id)))
    
  } else if(is.null(nTrain)) {
    # both undefined, use all training data
    nTrain <- length(unique(obs.id))
  }
  
  # Set distribution object - put in all possible additional parameters (this will generate warnings)
  if (is.character(distribution)){ distribution <- list(name=distribution) }
  dist_obj <- gbm_dist(distribution$name, ties=tied.times.method, strata=strata,
                       group=group, metric=distribution$metric,
                       max.rank=distribution$max.rank, prior_node_coeff_var=prior.node.coeff.var,
                       alpha=distribution$alpha, df=distribution$df, power=distribution$power)
  
  # Set up training parameters
  if(is.null(mFeatures)) mFeatures <- ncol(x)
  params <- training_params(num_trees=n.trees, interaction_depth=interaction.depth, min_num_obs_in_node=n.minobsinnode, 
                            shrinkage=shrinkage, bag_fraction=bag.fraction, id=obs.id,
                            num_train=nTrain, num_features=mFeatures)
  
  # Call gbm2 - formula defined from the data
  gbm_fit_obj <- gbm2(formula(data), distribution=dist_obj, data, weights=w, offset=offset,
                      train_params= params, var_monotone=var.monotone, var_names=var.names,
                      cv_folds=1, cv_class_stratify=FALSE, fold_id=NULL,
                      keep_gbm_data=keep.data, is_verbose=verbose)
  
  return(gbm_fit_obj)
}