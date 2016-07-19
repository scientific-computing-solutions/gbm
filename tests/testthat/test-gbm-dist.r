####################
# Author: James Hickey
#
# Series of test to validate the GBMDist objects.
#
####################


########## Definition ###############

context("Testing GBMDist Object Definition:")

test_that("Default distribution is Gaussian", {
  gbm_dist_obj <- gbm_dist()
  expect_equal(gbm_dist_obj$name, "Gaussian")
})

test_that("Check AdaBoost Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="AdaBoost")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("AdaBoostGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Bernoulli Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Bernoulli")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("BernoulliGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check CoxPH Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="CoxPH")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("CoxPHGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Gamma Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Gamma")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("GammaGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Gaussian Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Gaussian")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("GaussianGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Huberized Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Huberized")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("HuberizedGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Laplace Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Laplace")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("LaplaceGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Pairwise Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Pairwise")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("PairwiseGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Poisson Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Poisson")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("PoissonGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Quantile Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Quantile")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("QuantileGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check TDist Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="TDist")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("TDistGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Tweedie Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Tweedie")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("TweedieGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check AdaBoost distribution has the right elements - in right order", {
  gbm_dist_obj <- gbm_dist()
  expect_equal(names(gbm_dist_obj), c("name", "reorder", "ordering"))
})

test_that("Check AdaBoost distribution has the right elements - in right order", {
  gbm_dist_obj <- gbm_dist()
  expect_equal(names(gbm_dist_obj), c("name", "reorder", "ordering"))
})

test_that("Check Bernoulli distribution has the right elements - in right order", {
  gbm_dist_obj <- gbm_dist()
  expect_equal(names(gbm_dist_obj), c("name", "reorder", "ordering"))
})

test_that("Check CoxPH distribution has the right elements - order", {
  gbm_dist_obj <- gbm_dist()
  expect_equal(names(gbm_dist_obj), c("name", "reorder", "ordering"))
})

test_that("Check Gamma distribution has the right elements - in right order", {
  gbm_dist_obj <- gbm_dist()
  expect_equal(names(gbm_dist_obj), c("name", "reorder", "ordering"))
})

test_that("Check Gaussian distribution has the right elements - in right order", {
  gbm_dist_obj <- gbm_dist(name="Gaussian")
  expect_equal(names(gbm_dist_obj), c("name", "reorder", "ordering"))
})

test_that("Check Huberized distribution has the right elements - in right order", {
  gbm_dist_obj <- gbm_dist(name="Huberized")
  expect_equal(names(gbm_dist_obj), c("name", "reorder", "ordering"))
})

test_that("Check Laplace distribution has the right elements - in right order", {
  gbm_dist_obj <- gbm_dist(name="Laplace")
  expect_equal(names(gbm_dist_obj), c("name", "reorder", "ordering"))
})

test_that("Check Pairwise distribution has the right elements - in right order", {
  gbm_dist_obj <- gbm_dist(name="Pairwise")
  expect_equal(names(gbm_dist_obj), c("name", "reorder", "ordering", "metric", "group", "max.rank"))
})

test_that("Check Poisson distribution has the right elements - in right order", {
  gbm_dist_obj <- gbm_dist(name="Poisson")
  expect_equal(names(gbm_dist_obj), c("name", "reorder", "ordering"))
})

test_that("Check Quantile distribution has the right elements - in right order", {
  gbm_dist_obj <- gbm_dist(name="Quantile")
  expect_equal(names(gbm_dist_obj), c("name", "reorder", "ordering", "alpha"))
})

test_that("Check TDist distribution has the right elements - in right order", {
  gbm_dist_obj <- gbm_dist(name="TDist")
  expect_equal(names(gbm_dist_obj), c("name", "reorder", "ordering", "df"))
})

test_that("Check Tweedie distribution has the right elements - in right order", {
  gbm_dist_obj <- gbm_dist(name="Tweedie")
  expect_equal(names(gbm_dist_obj), c("name", "reorder", "ordering", "power"))
})

###### Warnings ######
context("Check warnings are thrown when too many arguments on construction of object")
test_that("Check warning for too many arguments - AdaBoost", {
  expect_warning(gbm_dist(name="AdaBoost", extra=1.0))
})

test_that("Check warning for too many arguments - Bernoulli", {
  expect_warning(gbm_dist(name="Bernoulli", extra=1.0))
})

test_that("Check warning for too many arguments - CoxPH", {
  expect_warning(gbm_dist(name="CoxPH", extra=1.0))
})

test_that("Check warning for too many arguments - Gamma", {
  expect_warning(gbm_dist(name="Gamma", extra=1.0))
})

test_that("Check warning for too many arguments - Gaussian", {
  expect_warning(gbm_dist(name="Gaussian", extra=1.0))
})

test_that("Check warning for too many arguments - Huberized", {
  expect_warning(gbm_dist(name="Huberized", extra=1.0))
})

test_that("Check warning for too many arguments - Laplace", {
  expect_warning(gbm_dist(name="Laplace", extra=1.0))
})

test_that("Check warning for too many arguments - Pairwise", {
  expect_warning(gbm_dist(name="Pairwise", extra=1.0))
})

test_that("Check warning for too many arguments - Poisson", {
  expect_warning(gbm_dist(name="Poisson", extra=1.0))
})

test_that("Check warning for too many arguments - Quantile", {
  expect_warning(gbm_dist(name="Quantile", extra=1.0))
})

test_that("Check warning for too many arguments - TDist", {
  expect_warning(gbm_dist(name="TDist", extra=1.0))
})

test_that("Check warning for too many arguments - Tweedie", {
  expect_warning(gbm_dist(name="Tweedie", extra=1.0))
})

##### Error checking ##### 
context("Check expect errors on construction if incorrect parameters provided")
test_that("Error thrown if unsupported distribution selected", {
  expect_error(gbm_dist("No sense in believing this will construct"))
})

test_that("Error thrown if 'ties' parameter is not a string- CoxPH", {
  expect_error(gbm_dist(name="CoxPH", ties=1.0))
  expect_error(gbm_dist(name="CoxPH", ties=Inf))
  expect_error(gbm_dist(name="CoxPH", ties=NA))
  expect_error(gbm_dist(name="CoxPH", ties=NULL))
})

test_that("Error thrown if strata not a vector of integers - CoxPH", {
  expect_error(gbm_dist(name="CoxPH", strata=c(1.2, 1.4, 1.5)))
  expect_error(gbm_dist(name="CoxPH", strata=NULL))
  expect_error(gbm_dist(name="CoxPH", strata="String"))
  expect_error(gbm_dist(name="CoxPH", strata=-0.1))
  expect_error(gbm_dist(name="CoxPH", strata=Inf))
})

test_that("Error thrown if sorted not a vector of integers - CoxPH", {
  expect_error(gbm_dist(name="CoxPH", sorted=c(1.2, 1.4, 1.5)))
  expect_error(gbm_dist(name="CoxPH", sorted=NULL))
  expect_error(gbm_dist(name="CoxPH", sorted="String"))
  expect_error(gbm_dist(name="CoxPH", sorted=-0.1))
  expect_error(gbm_dist(name="CoxPH", sorted=Inf))
})

test_that("Error thrown if prior coefficient of variation if not a finite double - CoxPH", {
  expect_error(gbm_dist(name="CoxPH", prior_node_coeff=Inf))
  expect_error(gbm_dist(name="CoxPH", prior_node_coeff="Nope"))
  expect_error(gbm_dist(name="CoxPH", prior_node_coeff=c(1.2, 3.4)))
})

test_that("Error thrown if max.rank is not a finite double greater than 0.0 - Pairwise", {
  expect_error(gbm_dist(name="Pairwise", metric="ndcg", max.rank=-0.1))
  expect_error(gbm_dist(name="Pairwise", metric="ndcg", max.rank=-0.1))
  expect_error(gbm_dist(name="Pairwise", metric="ndcg", max.rank="Stuff"))
  expect_error(gbm_dist(name="Pairwise", metric="ndcg", max.rank=c(1.0, 2.0)))
  expect_error(gbm_dist(name="Pairwise", metric="ndcg", max.rank=Inf))
  expect_error(gbm_dist(name="Pairwise", metric="ndcg", max.rank=NA))
  expect_error(gbm_dist(name="Pairwise", metric="ndcg", max.rank=NULL))
  
  expect_error(gbm_dist(name="Pairwise", metric="mrr", max.rank=-0.1))
  expect_error(gbm_dist(name="Pairwise", metric="mrr", max.rank=-0.1))
  expect_error(gbm_dist(name="Pairwise", metric="mrr", max.rank="Stuff"))
  expect_error(gbm_dist(name="Pairwise", metric="mrr", max.rank=c(1.0, 2.0)))
  expect_error(gbm_dist(name="Pairwise", metric="mrr", max.rank=Inf))
  expect_error(gbm_dist(name="Pairwise", metric="mrr", max.rank=NA))
  expect_error(gbm_dist(name="Pairwise", metric="mrr", max.rank=NULL))
})

test_that("Error thrown if max.rank is non-zero for conc or map - Pairwise", {
  expect_error(gbm_dist(name="Pairwise", metric="conc", max.rank=1.0))
  expect_error(gbm_dist(name="Pairwise", metric="map", max.rank=1.0))
})

test_that("Error thrown if group is not a string - Pairwise", {
  expect_error(gbm_dist(name="Pairwise", group=1.0))
  expect_error(gbm_dist(name="Pairwise", group=Inf))
  expect_error(gbm_dist(name="Pairwise", group=NA))
  expect_error(gbm_dist(name="Pairwise", group=NULL))
  expect_error(gbm_dist(name="Pairwise", group=c("Group1", "Group2")))
})

test_that("Error thrown if metric is not: ndcg, map, mrr or conc - Pairwise", {
  expect_error(gbm_dist(name="Pairwise", metric="Made-up"))
})

test_that("Error thrown if alpha specified is not a finite double between 0.0 and 1.0 - Quantile", {
  expect_error(gbm_dist(name="Quantile", alpha=2.0))
  expect_error(gbm_dist(name="Quantile", alpha=-0.01))
  expect_error(gbm_dist(name="Quantile", alpha="Character"))
  expect_error(gbm_dist(name="Quantile", alpha=Inf))
  expect_error(gbm_dist(name="Quantile", alpha=c(0.5, 0.1)))
  expect_error(gbm_dist(name="Quantile", alpha=NA))
  expect_error(gbm_dist(name="Quantile", alpha=NULL))
})

test_that("Error thrown if degrees of freedom specified is not a finite double > 0.0 - TDist", {
  expect_error(gbm_dist(name="TDist", df=-0.01))
  expect_error(gbm_dist(name="TDist", df="Character"))
  expect_error(gbm_dist(name="TDist", df=Inf))
  expect_error(gbm_dist(name="TDist", df=c(0.5, 0.1)))
  expect_error(gbm_dist(name="TDist", df=NA))
  expect_error(gbm_dist(name="TDist", df=NULL))
})

test_that("Error thrown if power specified is not a finite double > 0.0 - Tweedie", {
  expect_error(gbm_dist(name="Tweedie", power=-0.01))
  expect_error(gbm_dist(name="Tweedie", power="Character"))
  expect_error(gbm_dist(name="Tweedie", power=Inf))
  expect_error(gbm_dist(name="Tweedie", power=c(0.5, 0.1)))
  expect_error(gbm_dist(name="Tweedie", power=NA))
  expect_error(gbm_dist(name="Tweedie", power=NULL))
})

##### Default Parameters #####
context("Check default values of fields")
test_that("AdaBoost has no default ordering and reorder is FALSE", {
  expect_true(is.null(gbm_dist(name="AdaBoost")$ordering))
  expect_false(gbm_dist(name="AdaBoost")$reorder)
})

test_that("Bernoulli has no default ordering and reorder is FALSE", {
  expect_true(is.null(gbm_dist(name="Bernoulli")$ordering))
  expect_false(gbm_dist(name="Bernoulli")$reorder)
})

test_that("CoxPH has no default ordering and reorder is TRUE", {
  expect_true(is.null(gbm_dist(name="CoxPH")$ordering))
  expect_true(gbm_dist(name="CoxPH")$reorder)
})

test_that("Gamma has no default ordering and reorder is FALSE", {
  expect_true(is.null(gbm_dist(name="Gamma")$ordering))
  expect_false(gbm_dist(name="Gamma")$reorder)
})

test_that("Gaussian has no default ordering and reorder is FALSE", {
  expect_true(is.null(gbm_dist(name="Gaussian")$ordering))
  expect_false(gbm_dist(name="Gaussian")$reorder)
})

test_that("Laplace has no default ordering and reorder is FALSE", {
  expect_true(is.null(gbm_dist(name="Laplace")$ordering))
  expect_false(gbm_dist(name="Laplace")$reorder)
})

test_that("Pairwise has no default ordering and reorder is TRUE", {
  expect_true(is.null(gbm_dist(name="Pairwise")$ordering))
  expect_true(gbm_dist(name="Pairwise")$reorder)
})

test_that("Poisson has no default ordering and reorder is FALSE", {
  expect_true(is.null(gbm_dist(name="Poisson")$ordering))
  expect_false(gbm_dist(name="Poisson")$reorder)
})

test_that("Quantile has no default ordering and reorder is FALSE", {
  expect_true(is.null(gbm_dist(name="Quantile")$ordering))
  expect_false(gbm_dist(name="Quantile")$reorder)
})

test_that("TDist has no default ordering and reorder is FALSE", {
  expect_true(is.null(gbm_dist(name="TDist")$ordering))
  expect_false(gbm_dist(name="TDist")$reorder)
})

test_that("Tweedie has no default ordering and reorder is FALSE", {
  expect_true(is.null(gbm_dist(name="Tweedie")$ordering))
  expect_false(gbm_dist(name="Tweedie")$reorder)
})

test_that("CoxPH - defaults to 'efron', a prior coeff var of 1000, with NAs for sorted and strata", {
  expect_true(is.na(gbm_dist(name="CoxPH")$strata))
  expect_true(is.na(gbm_dist(name="CoxPH")$sorted))
  expect_equal(gbm_dist(name="CoxPH")$prior_node_coeff, 1000)
  expect_equal(gbm_dist(name="CoxPH")$ties, "efron")
})

test_that("Pairwise params default to - 'ndcg', max.rank=0 and group='query'", {
  expect_equal(gbm_dist(name="Pairwise")$metric, "ndcg")
  expect_equal(gbm_dist(name="Pairwise")$max.rank, 0)
  expect_equal(gbm_dist(name="Pairwise")$group, "query")
})

test_that("Quantile alpha defaults to 0.25", {
  expect_equal(gbm_dist(name="Quantile"), 0.25)
})

test_that("TDist df defaults to 4", {
  expect_equal(gbm_dist(name="TDist")$df, 4)
})

test_that("Tweedie defaults to dist with power = 1.5", {
  expect_equal(gbm_dist(name="Tweedie")$power, 1.5)
})