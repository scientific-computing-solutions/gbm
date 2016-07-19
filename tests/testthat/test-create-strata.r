####################
# Author: James Hickey
#
# Series of test to check the functionality that creates strata
#
####################

context("Testing Strata creation:")
test_that("Strata creation function requires GBMData and GBMDist objects", {
  # Require Surv to be available
  require(survival)
  
  # create some data
  set.seed(1)
  N <- 3000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
  tt.surv <- rexp(N,exp(f))
  tt.cens <- rexp(N,0.5)
  delta <- as.numeric(tt.surv <= tt.cens)
  tt <- apply(cbind(tt.surv,tt.cens),1,min)
  
  # throw in some missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  offset <- rep(0, N/2)
  
  Resp <- Surv(tt, delta)
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset)
  train_p <- training_params(id=rep(1, N), num_train = 1, num_features = 3)
  dist <- gbm_dist("CoxPH")
  
  # When not a GBMData object or GBMDist
  copy_data <- data
  copy_dist <- dist
  attr(data, "class") <- "NOTData"
  attr(dist, "class") <- "NOTDist"
  
  # Then an error is thrown
  expect_error(create_strata(data, train_p, copy_dist))
  expect_error(create_strata(copy_data, train_p, dist))
})

test_that("Distribution object remains unchanged if not CoxPH - strata remain undefined", {
  # Require Surv to be available
  require(survival)
  
  # create some data
  set.seed(1)
  N <- 3000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
  tt.surv <- rexp(N,exp(f))
  tt.cens <- rexp(N,0.5)
  delta <- as.numeric(tt.surv <= tt.cens)
  tt <- apply(cbind(tt.surv,tt.cens),1,min)
  
  # throw in some missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  offset <- rep(0, N/2)
  
  Resp <- Surv(tt, delta)
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset)
  train_p <- training_params(id=rep(1, N), num_train = 1, num_features = 3)
  
  # GIVEN Dist Not COXPH
  dist <- gbm_dist("AdaBoost")
  
  # When strata created
  copy_dist <- dist
  dist <- create_strata(data, train_p, dist)
  
  # Then dist object is unchanged
  expect_equal(copy_dist, dist)
})

test_that("Creating strata fills strata, time_order and sorted fields - CoxPH", {
  # Require Surv to be available
  require(survival)
  
  # create some data
  set.seed(1)
  N <- 3000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
  tt.surv <- rexp(N,exp(f))
  tt.cens <- rexp(N,0.5)
  delta <- as.numeric(tt.surv <= tt.cens)
  tt <- apply(cbind(tt.surv,tt.cens),1,min)
  
  # throw in some missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  offset <- rep(0, N/2)
  
  Resp <- Surv(tt, delta)
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset)
  train_p <- training_params(id=c(rep(1, N/2), rep(2, N/2)), num_train = 1, num_features = 3)
  
  # GIVEN Dist - COXPH
  dist <- gbm_dist("CoxPH")
  
  # When strata created
  dist <- create_strata(data, train_p, dist)
  
  # Then dist object is changed
  expect_equal(length(dist$strata), N)
  expect_equal(length(dist$time_order), N)
  expect_equal(nrow(dist$sorted), N)
})

test_that("Strata can only be created when response is Survival object - CoxPH", {
  # Require Surv to be available
  require(survival)
  
  # create some data
  set.seed(1)
  N <- 3000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
  tt.surv <- rexp(N,exp(f))
  tt.cens <- rexp(N,0.5)
  delta <- as.numeric(tt.surv <= tt.cens)
  tt <- apply(cbind(tt.surv,tt.cens),1,min)
  
  # throw in some missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  offset <- rep(0, N/2)
  
  Resp <- Surv(tt, delta)
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset)
  train_p <- training_params(id=c(rep(1, N/2), rep(2, N/2)), num_train = 1, num_features = 3)
  
  # GIVEN Dist - COXPH
  dist <- gbm_dist("CoxPH")
  
  # When response is not a Survival object
  attr(data$y, "type") <- "Not Right or Counting"
 
  # Then error thrown when creating strata
  expect_error(create_strata(data, train_p, dist))
})

test_that("If strata field in distribution object is NULL, all data are put in same strata - CoxPH", {
  # Require Surv to be available
  require(survival)
  
  # create some data
  set.seed(1)
  N <- 3000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
  tt.surv <- rexp(N,exp(f))
  tt.cens <- rexp(N,0.5)
  delta <- as.numeric(tt.surv <= tt.cens)
  tt <- apply(cbind(tt.surv,tt.cens),1,min)
  
  # throw in some missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  offset <- rep(0, N/2)
  
  Resp <- Surv(tt, delta)
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset)
  train_p <- training_params(id=c(rep(1, N/2), rep(2, N/2)), num_train = 1, num_features = 3)
  
  # GIVEN Dist - COXPH
  dist <- gbm_dist("CoxPH")
  
  # When creating strata with strata initialized to NULL
  dist$strata <- NULL
  dist <- create_strata(data, train_p, dist)
  
  # Then all examples put in same strata
  expect_equal(dist$strata, N/2)
})

test_that("The training responses are sorted according to strata and this order is stored in time_order", {
  # Require Surv to be available
  require(survival)
  
  # create some data
  set.seed(1)
  N <- 3000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
  tt.surv <- rexp(N,exp(f))
  tt.cens <- rexp(N,0.5)
  delta <- as.numeric(tt.surv <= tt.cens)
  tt <- apply(cbind(tt.surv,tt.cens),1,min)
  
  # throw in some missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  offset <- rep(0, N/2)
  
  Resp <- Surv(tt, delta)
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset)
  train_p <- training_params(id=c(rep(1, N/2), rep(2, N/2)), num_train = 1, num_features = 3)
  
  # GIVEN Dist - COXPH
  dist <- gbm_dist("CoxPH")
  
  # When strata are created
  dist <- create_strata(data, train_p, dist)
    
  # Then the training responses are sorted according to strata
  expect_equal(order(-data$y[1:N/2, 1]), dist$time_order)
})

test_that("Strata are sorted according to observation id if provided and strata is not NULL- CoxPH", {
  # Require Surv to be available
  require(survival)
  
  # create some data
  set.seed(1)
  N <- 3000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
  tt.surv <- rexp(N,exp(f))
  tt.cens <- rexp(N,0.5)
  delta <- as.numeric(tt.surv <= tt.cens)
  tt <- apply(cbind(tt.surv,tt.cens),1,min)
  
  # throw in some missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  offset <- rep(0, N/2)
  
  Resp <- Surv(tt, delta)
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset)
  train_p <- training_params(id=c(rep(1, N/2), rep(2, N/2)), num_train = 1, num_features = 3)
  
  # GIVEN Dist - COXPH
  dist <- gbm_dist("CoxPH")
  strata <- sample(1:5, N, replace=TRUE)
  dist$strata <- strata
  
  # When strata are created
  dist <- create_strata(data, train_p, dist)
  
  # Then ordered by id
  expect_equal(strata[order(train_p$id)], dist$strata)
})
