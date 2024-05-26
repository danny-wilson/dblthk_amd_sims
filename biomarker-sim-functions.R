# Load required libraries
library(glmnet)
library(combinat)
library(mvtnorm)
library(harmonicmeanp)
library(hommel)
# Load required source code
source("summary_mvMR_BF.R")

# Check my edited version has been loaded
stopifnot(exists("beta2_summary"))

# Class for data for analysis
setClass("bmsim_data",
	representation(
		id = "character",
		y = "matrix",
		x = "matrix",
		n = "integer",
		m = "integer",
		y.name = "character",
		x.names = "character"
	)
)

setClass("bmsim_pvalueTests",
	representation(
		method = "character",
		marginal.neglog10padj = "numeric",
		pairwise.neglog10padj = "numeric",
		truenull.neglog10padj = "numeric",
		headline.neglog10padj = "numeric"
	)
)

# Class for standardized analysis results
setClass("bmsim_analysisResults",
	representation(
		analysis = "character",
		data = "character",
		m = "integer",
		names = "character",
		estimate = "numeric",
		stderror = "numeric",
		# Frequentist Bonferroni-adjusted
		signif.neglog10padj = "numeric",
		pairwise.signif.neglog10padj = "numeric",
		truenull.signif.neglog10padj = "numeric",
		headline.signif.neglog10padj = "numeric",
		# Frequentist BH
		# Frequentist HMP
		# Frequentist Multilevel Simes
		# Frequentist Hommel
		# Frequentist Cauchy
		# Frequentist E-values
		pvalueTests = "list",
		# Bayesian
		signif.log10po = "numeric",
		pairwise.signif.log10po = "numeric",
		truenull.signif.log10po = "numeric",
		headline.signif.log10po = "numeric",
		time.secs = "numeric"
	)
)

# Class for standardized performance metrics
setClass("bmsim_test_performance",
	representation(
		# Frequentist test performance
		bias = "numeric",
		typeI = "numeric",
		typeII = "numeric",
		call = "numeric",
		familywiseI = "numeric",
		familywiseII = "numeric",
		pfdr = "numeric",
		fdr = "numeric",
		fndr = "numeric",
		fpr = "numeric",
		fnr = "numeric",
		tp = "numeric",
		fp = "numeric",
		tn = "numeric",
		fn = "numeric"
	),
	prototype(
		# Frequentist test performance
		bias = double(0),
		typeI = double(0),
		typeII = double(0),
		call = double(0),
		familywiseI = double(0),
		familywiseII = double(0),
		pfdr = double(0),
		fdr = double(0),
		fndr = double(0),
		fpr = double(0),
		fnr = double(0),
		tp = double(0),
		fp = double(0),
		tn = double(0),
		fn = double(0)
	)
)

# Class for standardized performance metrics
setClass("bmsim_performance",
	representation(
		analysis = "character",
		data = "character",
		m = "integer",
		names = "character",
		# Estimator performance
		estimate.bias = "numeric",
		estimate.L1 = "numeric",
		estimate.L2 = "numeric",
		# Std Error performance
		stderr.coverage = "numeric",
		# Frequentist test performance
		freqt.marginal = "bmsim_test_performance",
		freqt.pairwise = "bmsim_test_performance",
		freqt.truenull = "bmsim_test_performance",
		freqt.headline = "bmsim_test_performance",
		# p-value test performance,
		pvalueTests.marginal = "list",
		pvalueTests.pairwise = "list",
		pvalueTests.truenull = "list",
		pvalueTests.headline = "list",
		# Bayesian test performance
		bayes.marginal = "bmsim_test_performance",
		bayes.pairwise = "bmsim_test_performance",
		bayes.truenull = "bmsim_test_performance",
		bayes.headline = "bmsim_test_performance",
		# Prediction performance
		prediction.L1 = "numeric",
		prediction.L2 = "numeric",
		# Computation time (seconds)
		time.secs = "numeric"
	),
	prototype(
		analysis = character(0),
		data = character(0),
		m = integer(0),
		names = character(0),
		# Estimator performance
		estimate.bias = double(0),
		estimate.L1 = double(0),
		estimate.L2 = double(0),
		# Std Error performance
		stderr.coverage = double(0),
		# Frequentist test performance
		freqt.marginal = new("bmsim_test_performance"),
		freqt.pairwise = new("bmsim_test_performance"),
		freqt.truenull = new("bmsim_test_performance"),
		freqt.headline = new("bmsim_test_performance"),
		# p-value test performance,
		pvalueTests.marginal = list(),
		pvalueTests.pairwise = list(),
		pvalueTests.truenull = list(),
		pvalueTests.headline = list(),
		# Bayesian test performance
		bayes.marginal = new("bmsim_test_performance"),
		bayes.pairwise = new("bmsim_test_performance"),
		bayes.truenull = new("bmsim_test_performance"),
		bayes.headline = new("bmsim_test_performance"),
		# Prediction performance
		prediction.L1 = double(0),
		prediction.L2 = double(0),
		# Computation time (seconds)
		time.secs = double(0)
	)
)

# Functions to validate classes
setGeneric("validate", function(obj) standardGeneric("validate"))
setMethod("validate", "bmsim_data", function(obj) {
	stopifnot(is(obj, "bmsim_data"))
	stopifnot(nrow(obj@y)==obj@n)
	stopifnot(ncol(obj@y)==1)
	stopifnot(nrow(obj@x)==obj@n)
	stopifnot(ncol(obj@x)==obj@m)
	stopifnot(length(obj@y.name)==1)
	stopifnot(length(obj@x.names)==obj@m)
	stopifnot(colnames(obj@y)==obj@y.name)
	stopifnot(colnames(obj@x)==obj@x.names)
	stopifnot(is.na(match(obj@y.name, obj@x.names)))
	return(TRUE)
})
setMethod("validate", "bmsim_analysisResults", function(obj) {
	stopifnot(is(obj, "bmsim_analysisResults"))
	stopifnot(length(obj@analysis)==1)
	stopifnot(length(obj@data)==1)
	stopifnot(length(obj@m)==1)
	stopifnot(length(obj@names)==obj@m)
	stopifnot(length(obj@estimate)==obj@m)
	stopifnot(length(obj@stderror)==obj@m)
	stopifnot(length(obj@signif.neglog10padj)==obj@m)
	stopifnot(length(obj@signif.log10po)==obj@m)
	stopifnot(length(obj@time.secs)==1)
	return(TRUE)
})

# Prediction function
setMethod("predict", "bmsim_analysisResults", function(object, newdata) {
	stopifnot(is(newdata, "bmsim_data"))
	validate(newdata)
	stopifnot(object@m==newdata@m)
	as.vector(newdata@x %*% matrix(object@estimate, ncol=1))
})

#   Define functions to compute performance metrics:
#     Point estimates:
#		Bias
#       Mean absolute error (L1)
#       Root mean square error loss (L2)
#     95% interval estimates:
#       Coverage
#     Hypothesis tests:
#       Type I error loss
#       Familywise type I error loss
#       Type II error loss
#       Familywise type II error loss
#       False discovery rate
#       False non-discovery rate
#		False positive rate
#		False negative rate
#     Prediction:
#       Mean absolute error (L1)
#       RMSE (L2)
#     Compute:
#       Compute time

performance.bias = function(estimate, parameter) estimate-parameter
performance.L1 = function(estimate, parameter) mean(abs(estimate-parameter))
performance.L2 = function(estimate, parameter) sqrt(mean((estimate-parameter)^2))
performance.stderr.coverage = function(estimate, stderr, parameter, alpha=0.05) {
	1*(qnorm(alpha/2, estimate, stderr) <= parameter & parameter <= qnorm(1-alpha/2, estimate, stderr))
}
performance.typeI = function(test.signif, param.signif) (1-param.signif)*test.signif
performance.typeII = function(test.signif, param.signif) param.signif*(1-test.signif)
performance.familywise = function(elementwise) 1*(sum(elementwise)>0)
performance.error.rate = function(num, den) sum(num)/pmax(1, sum(den))
na.if = function(cond) ifelse(cond, as.double(NA), as.double(1))

calc.test.performance = function(test.signif, param.signif, param.names=NULL) {
	stopifnot(length(test.signif)==length(param.signif))
	if(!is.null(param.names)) stopifnot(length(test.signif)==length(param.names))
	typeI = performance.typeI(test.signif, param.signif)
	typeII = performance.typeII(test.signif, param.signif)
	ret = new("bmsim_test_performance",
		bias = performance.bias(test.signif, param.signif),
		typeI = typeI,
		typeII = typeII,
		call = test.signif,
		familywiseI = performance.familywise(typeI) * na.if(all(param.signif==1)),
		familywiseII = performance.familywise(typeII) * na.if(all(param.signif==0)),
		pfdr = performance.error.rate(typeI, test.signif) * na.if(all(test.signif==0)),
		fdr = performance.error.rate(typeI, test.signif),
		fndr = performance.error.rate(typeII, 1-test.signif),
		fpr = performance.error.rate(typeI, 1-param.signif) * na.if(all(param.signif==1)),
		fnr = performance.error.rate(typeII, param.signif) * na.if(all(param.signif==0)),
		tp = sum( (1-typeI)*test.signif),
		fp = sum(    typeI *test.signif),
		tn = sum((1-typeII)*(1-test.signif)),
		fn = sum(   typeII *(1-test.signif))
	)
	if(!is.null(param.names)) {
		names(ret@bias) <- param.names
		names(ret@typeI) <- param.names
		names(ret@typeII) <- param.names
		names(ret@call) <- param.names
	}
	return(ret)
}

setGeneric("performance", function(obj, ...) standardGeneric("performance"))
setMethod("performance", "bmsim_analysisResults", function(obj, parameter, newdata=NULL, thresh.neglog10padj=NA, thresh.log10po=NA) {
	validate(obj)
	stopifnot(length(parameter)==obj@m)
	if(!is.null(newdata)) {
		stopifnot(newdata@m==obj@m)
		validate(newdata)
	}
	
	ret = new("bmsim_performance",
		analysis = obj@analysis,
		data = obj@data,
		m = obj@m,
		names = obj@names
	)
	
	# Estimate
	ret@estimate.bias = performance.bias(obj@estimate, parameter)
	ret@estimate.L1 = performance.L1(obj@estimate, parameter)
	ret@estimate.L2 = performance.L2(obj@estimate, parameter)
	names(ret@estimate.bias) <- ret@names

	# Std Error
	ret@stderr.coverage = performance.stderr.coverage(obj@estimate, obj@stderror, parameter)
	names(ret@stderr.coverage) <- ret@names

	# Hypothesis tests: marginal
	param.signif = 1*(parameter!=0)
	
	freqt.signif = 1*(obj@signif.neglog10padj >= thresh.neglog10padj)
	ret@freqt.marginal = calc.test.performance(freqt.signif, param.signif, ret@names)

	bayes.signif = 1*(obj@signif.log10po >= thresh.log10po)
	ret@bayes.marginal = calc.test.performance(bayes.signif, param.signif, ret@names)

	ret@pvalueTests.marginal = lapply(obj@pvalueTests, function(test) {
		ret = new("bmsim_test_performance")
		if(length(test@marginal.neglog10padj)>0) ret = calc.test.performance(1*(test@marginal.neglog10padj >= thresh.neglog10padj), param.signif, names(test@marginal.neglog10padj))
		return(ret)
	})

	# Hypothesis tests: pairwise
	if(length(obj@pairwise.signif.neglog10padj)>0) {
		param.signif = outer(param.signif, param.signif, FUN="|")
		param.signif = 1*(param.signif[lower.tri(param.signif)])
		
		pairwise.freqt.signif = 1*(obj@pairwise.signif.neglog10padj >= thresh.neglog10padj)
		ret@freqt.pairwise = calc.test.performance(pairwise.freqt.signif, param.signif, names(obj@pairwise.signif.neglog10padj))

		pairwise.bayes.signif = 1*(obj@pairwise.signif.log10po >= thresh.log10po)
		ret@bayes.pairwise = calc.test.performance(pairwise.bayes.signif, param.signif, names(obj@pairwise.signif.log10po))
		
		ret@pvalueTests.pairwise = lapply(obj@pvalueTests, function(test) {
			ret = new("bmsim_test_performance")
			if(length(test@pairwise.neglog10padj)>0) ret = calc.test.performance(1*(test@pairwise.neglog10padj >= thresh.neglog10padj), param.signif, names(test@pairwise.neglog10padj))
			return(ret)
		})
	}

	# Hypothesis tests: truenull
	if(length(obj@truenull.signif.neglog10padj)>0) {
		# Scalar; members of the true null set are never truly significant
		param.signif = 0

		truenull.freqt.signif = 1*(obj@truenull.signif.neglog10padj >= thresh.neglog10padj)
		ret@freqt.truenull = calc.test.performance(truenull.freqt.signif, param.signif, names(obj@truenull.signif.neglog10padj))

		truenull.bayes.signif = 1*(obj@truenull.signif.log10po >= thresh.log10po)
		ret@bayes.truenull = calc.test.performance(truenull.bayes.signif, param.signif, names(obj@truenull.signif.log10po))
		
		ret@pvalueTests.truenull = lapply(obj@pvalueTests, function(test) {
			ret = new("bmsim_test_performance")
			if(length(test@truenull.neglog10padj)>0) ret = calc.test.performance(1*(test@truenull.neglog10padj >= thresh.neglog10padj), param.signif, names(test@truenull.neglog10padj))
			return(ret)
		})
	}

	# Hypothesis tests: headline
	if(length(obj@headline.signif.neglog10padj)>0) {
		param.signif = 1*any(parameter!=0)
		
		headline.freqt.signif = 1*(obj@headline.signif.neglog10padj >= thresh.neglog10padj)
		ret@freqt.headline = calc.test.performance(headline.freqt.signif, param.signif, names(obj@headline.signif.neglog10padj))

		headline.bayes.signif = 1*(obj@headline.signif.log10po >= thresh.log10po)
		ret@bayes.headline = calc.test.performance(headline.bayes.signif, param.signif, names(obj@headline.signif.log10po))
		
		ret@pvalueTests.headline = lapply(obj@pvalueTests, function(test) {
			ret = new("bmsim_test_performance")
			if(length(test@headline.neglog10padj)>0) ret = calc.test.performance(1*(test@headline.neglog10padj >= thresh.neglog10padj), param.signif, names(test@headline.neglog10padj))
			return(ret)
		})
	}

	# Prediction
	if(!is.null(newdata)) {
		prediction = predict(obj, newdata)
		ret@prediction.L1 = performance.L1(prediction, newdata@y)
		ret@prediction.L2 = performance.L2(prediction, newdata@y)
	} else {
		ret@prediction.L1 = as.double(NA)
		ret@prediction.L2 = as.double(NA)
	}
	
	# Computation time (seconds)
	ret@time.secs = obj@time.secs
	
	return(ret)
})

# Methods to calculate pvalueTests
calc.Bonferroni = function(p.unadj, nu, params=NULL) {
	stopifnot(nu>=length(p.unadj))
	if(!is.null(params)) {
		stopifnot(length(params)==length(p.unadj))
		stopifnot(all(!is.na(params)))
	}
	p.adj = p.adjust(p.unadj, method="bonferroni", n=nu)
	names(p.adj) <- names(p.unadj)
	p.pair = outer(p.adj, p.adj, pmin)
	names.pair = outer(names(p.unadj), names(p.unadj), Vectorize(function(x, y) paste(x, y, sep=" | ")))
	p.pair = p.pair[lower.tri(p.pair)]
	names(p.pair) <- names.pair[lower.tri(names.pair)]
	p.truenull = as.double(NA)
	if(!is.null(params)) {
		p.truenull = min(p.adj[params==0])
		names(p.truenull) <- paste(names(p.unadj)[params==0], collapse=" | ")
	}
	p.headline = min(p.adj)
	names(p.headline) <- paste(names(p.unadj), collapse=" | ")
	new("bmsim_pvalueTests",
		method = "Bonferroni adjustment",
		marginal.neglog10padj = -log10(p.adj),
		pairwise.neglog10padj = -log10(p.pair),
		truenull.neglog10padj = -log10(p.truenull),
		headline.neglog10padj = -log10(p.headline)
	)
}

# False discovery rate control; Benjamini & Hochberg (1995)
calc.BH = function(p.unadj, nu, params=NULL) {
	stopifnot(nu>=length(p.unadj))
	if(!is.null(params)) {
		stopifnot(length(params)==length(p.unadj))
		stopifnot(all(!is.na(params)))
	}
	p.pair = outer(p.unadj, p.unadj, Vectorize(function(x, y) min(p.adjust(c(x, y), method="BH", n=nu))))
	names.pair = outer(names(p.unadj), names(p.unadj), Vectorize(function(x, y) paste(x, y, sep=" | ")))
	p.pair = p.pair[lower.tri(p.pair)]
	names(p.pair) <- names.pair[lower.tri(names.pair)]
	p.truenull = as.double(NA)
	if(!is.null(params)) {
		p.truenull = min(p.adjust(p.unadj[params==0], method="BH", n=nu))
		names(p.truenull) <- paste(names(p.unadj)[params==0], collapse=" | ")
	}
	p.headline = min(p.adjust(p.unadj, method="BH", n=nu))
	names(p.headline) <- paste(names(p.unadj), collapse=" | ")
	new("bmsim_pvalueTests",
		method = "FDR control; Benjamini and Hochberg 1995",
		marginal.neglog10padj = -log10(p.adjust(p.unadj, method="BH", n=nu)),
		pairwise.neglog10padj = -log10(p.pair),
		truenull.neglog10padj = -log10(p.truenull),
		headline.neglog10padj = -log10(p.headline)
	)
}

# Harmonic mean p-value procedure; Wilson (2019); anticonservative for non-small p
calc.HMP = function(p.unadj, nu, min.p.unadj=1e-308, params=NULL) {
	stopifnot(nu>=length(p.unadj))
	if(!is.null(params)) {
		stopifnot(length(params)==length(p.unadj))
		stopifnot(all(!is.na(params)))
	}
	p.unadj = pmax(min.p.unadj, p.unadj)
	p.adj = sapply(p.unadj, Vectorize(function(p) p.hmp(c(p, rep(1, nu-1)), L=nu)))
	names(p.adj) <- names(p.unadj)
	p.pair = outer(p.unadj, p.unadj, Vectorize(function(x, y) p.hmp(c(x, y, rep(1, nu-2)), L=nu)))
	names.pair = outer(names(p.unadj), names(p.unadj), Vectorize(function(x, y) paste(x, y, sep=" | ")))
	p.pair = p.pair[lower.tri(p.pair)]
	names(p.pair) <- names.pair[lower.tri(names.pair)]
	p.truenull = as.double(NA)
	if(!is.null(params)) {
		p.truenull = p.hmp(c(p.unadj[params==0], rep(1, sum(params!=0))), L=nu)
		names(p.truenull) <- paste(names(p.unadj)[params==0], collapse=" | ")
	}
	p.headline = p.hmp(c(p.unadj, rep(1, nu-length(p.unadj))), L=nu)
	names(p.headline) <- paste(names(p.unadj), collapse=" | ")
	ret = new("bmsim_pvalueTests",
		method = "Harmonic mean p-value procedure; Wilson 2019",
		marginal.neglog10padj = -log10(p.adj),
		pairwise.neglog10padj = -log10(p.pair),
		truenull.neglog10padj = -log10(p.truenull),
		headline.neglog10padj = -log10(p.headline)
	)
	return(ret)
}

# Simes' test (analogous to harmonic mean p-value test)
p.Simes = function(p, w = NULL, L = NULL, w.sum.tolerance = 1e-6, multilevel = TRUE) {
  if(is.null(L) & multilevel) {
	warning("L not specified: for multilevel testing set L to the total number of individual p-values")
	L = length(p)
  }
  if(length(p) == 0) return(NA)
  if(length(p) > L) stop("The number of p-values cannot exceed L")
  if(is.null(w)) {
	w = rep(1/L,length(p))
  } else {
	if(any(w<0)) stop("No weights can be negative")
	if(length(w)!=length(p)) stop("When specified, length of w must equal length of p")
  }
  w.sum = sum(w)
  if(w.sum>1+w.sum.tolerance) {
	stop("Weights cannot exceed 1")
  }
  pstar = p/w
  return(c(p.Simes = min(sort(pstar)/(1:length(p)))))
}
# Multilevel Simes test; Wilson (2020)
calc.Simes = function(p.unadj, nu, params=NULL) {
	stopifnot(nu>=length(p.unadj))
	if(!is.null(params)) {
		stopifnot(length(params)==length(p.unadj))
		stopifnot(all(!is.na(params)))
	}
	p.pair = outer(p.unadj, p.unadj, Vectorize(function(x, y) p.Simes(c(x, y), L=nu)))
	names.pair = outer(names(p.unadj), names(p.unadj), Vectorize(function(x, y) paste(x, y, sep=" | ")))
	p.pair = p.pair[lower.tri(p.pair)]
	names(p.pair) <- names.pair[lower.tri(names.pair)]
	p.truenull = as.double(NA)
	if(!is.null(params)) {
		p.truenull = p.Simes(c(p.unadj[params==0], rep(1, sum(params!=0))), L=nu)
		names(p.truenull) <- paste(names(p.unadj)[params==0], collapse=" | ")
	}
	p.headline = p.Simes(p.unadj, L=nu)
	names(p.headline) <- paste(names(p.unadj), collapse=" | ")
	ret = new("bmsim_pvalueTests",
		method = "Multilevel Simes test; Wilson 2020",
		marginal.neglog10padj = -log10(p.unadj*nu),
		pairwise.neglog10padj = -log10(p.pair),
		truenull.neglog10padj = -log10(p.truenull),
		headline.neglog10padj = -log10(p.headline)
	)
	return(ret)
}

# Multilevel Simes test; Hommel (1988), efficient implementation by Goeman, Meijer, Krebs, Solari (2019)
calc.Hommel = function(p.unadj, nu, params=NULL) {
	m = length(p.unadj)
	stopifnot(nu>=m)
	if(!is.null(params)) {
		stopifnot(length(params)==length(p.unadj))
		stopifnot(all(!is.na(params)))
	}
	# Pad the p-values with the nu-m unobserved values, taking the worst case
	p.fill = c(p.unadj, rep(1.0, nu-m))
	p.adj = p.adjust(hommel(p.fill))[1:m]
	p.pair = outer(p.adj, p.adj, pmin)
	names.pair = outer(names(p.unadj), names(p.unadj), Vectorize(function(x, y) paste(x, y, sep=" | ")))
	p.pair = p.pair[lower.tri(p.pair)]
	names(p.pair) <- names.pair[lower.tri(names.pair)]
	p.truenull = as.double(NA)
	if(!is.null(params)) {
		p.truenull = min(p.adj[params==0])
		names(p.truenull) <- paste(names(p.unadj)[params==0], collapse=" | ")
	}
	p.headline = min(p.adj)
	names(headline) <- paste(names(p.unadj), collapse=" | ")
	ret = new("bmsim_pvalueTests",
		method = "Multilevel Simes test; Hommel 1988, Goeman et al (2019)",
		marginal.neglog10padj = -log10(p.adj),
		pairwise.neglog10padj = -log10(p.pair),
		truenull.neglog10padj = -log10(p.truenull),
		headline.neglog10padj = -log10(p.headline)
	)
	return(ret)
}

# Cauchy combination test
p.Cauchy = function(p, w=rep(1/length(p), length(p))) {
  t0 = sum(w*tan((0.5-p)*pi))
  0.5 - atan(t0)/pi
}
# Cauchy combination test; Liu and Xie (2020)
# Implemented as an expectation over excluded values: assuming they are less significant than
# the included values and follow a truncated uniform(0,1) distribution.
# Otherwise by forcing them to be 1, the combined p-value would be forced to be 0 or 1.
calc.Cauchy = function(p.unadj, nu, nsim=10000) {
	stopifnot(nu>=length(p.unadj))
	p.max = pmin(1, max(p.unadj))
	headline = mean(replicate(nsim, {p.Cauchy(c(p.unadj, runif(nu-length(p.unadj), p.max, 1)))}))
	ret = new("bmsim_pvalueTests",
		method = "Cauchy combination test; Liu and Xie 2020",
		marginal.neglog10padj = double(0),
		pairwise.neglog10padj = double(0),
		truenull.neglog10padj = double(0),
		headline.neglog10padj = -log10(headline)
	)
	names(ret@headline.neglog10padj) <- paste(names(p.unadj), collapse=" | ")
	return(ret)
}

# sum of logarithms function
logsum = function(x,na.rm=FALSE) {
	mx = max(x,na.rm=na.rm)
	log(sum(exp(x-mx),na.rm=na.rm))+mx
}
# p-value to e-value
logp2loge.value = function(logp, kappa=0.1) log(kappa)+(kappa-1)*logp
# e-value to p-value
e2p.value = function(eval) pmin(1, 1/eval)
# E-value test
p.evalue = function(p, L=length(p), kappa=0.1) {
	logBF = logp2loge.value(log(p))
	evalue = exp(logsum(logBF))/L
	# Under the null, the following is uniform or more conservative by Markov's inequality
	# assuming the expectation of the BF is 1 or less under the null
	p = e2p.value(evalue)
	return(p)
}
# Combined test for evalues
combined.evalues = function(log.evalues, K=length(log.evalues), simple=FALSE) {
	stopifnot(K>=length(log.evalues))
	if(simple) {
		evalue = exp(logsum(log.evalues))/K
		# Under the null, the following is uniform or more conservative by Markov's inequality
		# assuming the expectation of the BF is 1 or less under the null
	} else {
		# Implement Algorithm 1 of Vovk and Wang (2021)
		pi = order(log.evalues, decreasing=FALSE)
		log.e.ordered = log.evalues[pi]
		S = sapply(1:K, function(i) exp(logsum(log.e.ordered[1:i])))
		e.star = rep(as.double(NA), K)
		for(k in 1:K) {
			e.star[pi][k] = exp(log.e.ordered)[k]
			if(k>1) for(i in 1:(k-1)) {
				e = (exp(log.e.ordered)[k]+S[i])/(i+1)
				if(e < e.star[pi][k]) e.star[pi][k] = e
			}
		}
		evalue = e.star
	}
	# Apply the unique e to p calibrator
	p = e2p.value(evalue)
	return(p)
}
# An evalue procedure; based on Vovk and Wang (2021)
# Multilevel property based on specifying nu to be sum of all tests to be combined
calc.evalue = function(p.unadj, nu, kappa=0.1, params=NULL) {
	stopifnot(nu>=length(p.unadj))
	if(!is.null(params)) {
		stopifnot(length(params)==length(p.unadj))
		stopifnot(all(!is.na(params)))
	}
	p.adj = sapply(p.unadj, Vectorize(function(p) p.evalue(p, L=nu, kappa=kappa)))
	names(p.adj) <- names(p.unadj)
	p.pair = outer(p.unadj, p.unadj, Vectorize(function(x, y) p.evalue(c(x, y), L=nu, kappa=kappa)))
	names.pair = outer(names(p.unadj), names(p.unadj), Vectorize(function(x, y) paste(x, y, sep=" | ")))
	p.pair = p.pair[lower.tri(p.pair)]
	names(p.pair) <- names.pair[lower.tri(names.pair)]
	p.truenull = as.double(NA)
	if(!is.null(params)) {
		p.truenull = p.evalue(p.unadj[params==0], L=nu, kappa=kappa)
		names(p.truenull) <- paste(names(p.unadj)[params==0], collapse=" | ")
	}
	p.headline = p.evalue(c(p.unadj, rep(1, nu-length(p.unadj))), L=nu, kappa=kappa)
	names(p.headline) <- paste(names(p.unadj), collapse=" | ")
	ret = new("bmsim_pvalueTests",
		method = "An E-value combination test for p-values; based on Vovk and Wang 2021",
		marginal.neglog10padj = -log10(p.adj),
		pairwise.neglog10padj = -log10(p.pair),
		truenull.neglog10padj = -log10(p.truenull),
		headline.neglog10padj = -log10(p.headline)
	)
	return(ret)
}

# Multilevel procedure to convert BFs to p-values
# Warning: Multiple testing correction not applied
calc.evalue.BF2p = function(marginal.log10bf, pairwise.log10bf, truenull.log10bf, headline.log10bf, nu=length(marginal.log10bf)) {
	ret = new("bmsim_pvalueTests",
		method = "An E-value conversion of BFs to p-values; based on Vovk and Wang 2021",
		marginal.neglog10padj = -marginal.log10bf,
		pairwise.neglog10padj = -pairwise.log10bf,
		truenull.neglog10padj = -truenull.log10bf,
		headline.neglog10padj = -headline.log10bf
	)
	return(ret)
}

# Load test data 'amd_example' based on demo_AMD
load.data = function() {
	load("amd_example")
	
	betaX = amd_example$betaX
	amd_beta = amd_example$amd_beta
	amd_se = amd_example$amd_se
	rs = amd_example$annotate[,1]
	genes = amd_example$annotate[,7]
	rf = colnames(betaX)
	
	# There was one influential variant in the LIPC gene region and two outliers in the APOE and FUT2 gene region. We are removing these 3 data points for the following analysis.
	
	LIPC = which(genes == "LIPC")
	FUT2 = which(genes == "FUT2")
	APOE = which(genes == "APOE")
	exclude_vec = c(LIPC,FUT2,APOE)
	betaX=betaX[-exclude_vec,]
	amd_beta = amd_beta[-exclude_vec]
	amd_se = amd_se[-exclude_vec]
	rs = rs[-exclude_vec]
	genes = genes[-exclude_vec]
	
	# Next, we perform an inverse variance weighting (IVW) based on the standard error of the amd beta effect estimates prior to subsequent analysis.
	
	betaX_ivw = betaX / amd_se
	amd_beta_ivw = matrix(amd_beta / amd_se, ncol=1)
	biomarkers = gsub("^beta_", "", colnames(betaX_ivw))
	colnames(amd_beta_ivw) <- "amd"
	colnames(betaX_ivw) <- biomarkers
	
	data = new("bmsim_data",
		id = "demo_amd_example_145_snps",
		y = amd_beta_ivw,
		x = as.matrix(betaX_ivw),
		n = nrow(betaX_ivw),
		m = ncol(betaX_ivw),
		y.name = colnames(amd_beta_ivw),
		x.names = colnames(betaX_ivw)
	)
	validate(data)
	return(data)
}

# Conduct univariable MR
univariable.MR = function(data, nu=data@m) {
	validate(data)
	ret = new("bmsim_analysisResults",
		analysis = "univariable Mendelian randomization",
		data = data@id,
		m = data@m,
		names = data@x.names,
		estimate = rep(as.double(NA),data@m),
		stderror = rep(as.double(NA),data@m),
		signif.neglog10padj = rep(as.double(NA),data@m),
		signif.log10po = rep(as.double(NA),data@m)
	)
	start_time = Sys.time()
	for(j in 1:data@m) {
		dataf = data.frame(y=data@y[,1], x=data@x[,j])
		coefs = summary(lm(y ~ 0 + x, data=dataf))$coefficients
		ret@estimate[j] = coefs[1,"Estimate"]
		ret@stderror[j] = coefs[1,"Std. Error"]
		ret@signif.neglog10padj[j] = -log10(nu * coefs[1,"Pr(>|t|)"])
	}
	ret@time.secs = as.double(difftime(Sys.time(), start_time, units="secs"))
	validate(ret)
	return(ret)
}

# Filter test data based on univariable associations
filter.data = function(data, m=as.integer(15), results.univariable.MR=NULL) {
	validate(data)
	stopifnot(is(m,"integer"))
	if(is.null(results.univariable.MR)) {
		results.univariable.MR = univariable.MR(full.data)
	} else {
		validate(results.univariable.MR)
		stopifnot(results.univariable.MR@data == data@id)
		stopifnot(!any(is.na(results.univariable.MR@signif.neglog10padj)))
	}
	rk = rank(-results.univariable.MR@signif.neglog10padj, ties.method="first")
	keep = rk<=m
	ret = new("bmsim_data",
		id = paste0(data@id, "_", m, "_biomarkers"),
		y = data@y,
		x = data@x[,keep],
		n = data@n,
		m = m,
		y.name = data@y.name,
		x.names = data@x.names[keep]
	)
	validate(ret)
	return(ret)
}

# Exhaustive MR-BMA analysis
mr.bma.x = function(data, sigma=0.5, prior_prob=0.1, nsim=10, nu=data@m) {
	validate(data)
	start_time = Sys.time()
	
	rs = as.character(1:data@n)
	mr.input = new("mvMRInput", betaX = data@x, betaY = data@y, snps=rs, exposure=data@x.names, outcome = "y")
	mr.output = summary_mvMR_BF(mr.input, sigma=sigma, prior_prob=prior_prob, calc.se=TRUE)
	check_time = Sys.time()
	
	pval = rep(as.double(NA), data@m)
	if(nsim>0) {
		mr.input.sim = mr.input
		pp.sim = matrix(as.double(NA), nsim, data@m)
		for(sim in 1:nsim) {
			mr.input.sim@betaY = matrix(sample(mr.input@betaY), ncol=1)
			mr.output.sim = summary_mvMR_BF(mr.input.sim, sigma=sigma, prior_prob=prior_prob, calc.se=TRUE)
			pp.sim[sim,] = mr.output.sim@pp_marginal
		}
		pval = rowMeans(t(pp.sim)>=mr.output@pp_marginal)
	}
	
	# Bayesian model averaged result
	pp2log10po = function(pp, eps=1e-308) log10((eps + pp)/((eps + 1-pp)))
	ret = new("bmsim_analysisResults",
		analysis = "Bayesian model averaged multivariable Mendelian randomization with MR-BMA",
		data = data@id,
		m = data@m,
		names = data@x.names,
		estimate = mr.output@BMAve_Estimate,
		stderror = mr.output@BMAve_StdError,
		signif.neglog10padj = -log10(nu * pval),
		signif.log10po = pp2log10po(mr.output@pp_marginal),
		pvalueTests = list(
			"Evalue.BF2p" = calc.evalue.BF2p(pp2log10po(mr.output@pp_marginal) - log10(prior_prob/(1-prior_prob)), NA, NA, NA, nu)
		),
		time.secs = as.double(difftime(Sys.time(), start_time, units="secs"))
	)
	names(ret@estimate) <- ret@names
	names(ret@stderror) <- ret@names
	names(ret@signif.neglog10padj) <- ret@names
	names(ret@signif.log10po) <- ret@names
	validate(ret)
	
	# Bayesian MAP result
	ret.map = new("bmsim_analysisResults",
		analysis = "Bayesian model selection multivariable Mendelian randomization with MR-BMA",
		data = data@id,
		m = data@m,
		names = data@x.names,
		estimate = mr.output@BestModel_Estimate,
		stderror = mr.output@BestModel_StdError,
		signif.neglog10padj = rep(as.double(NA), len=length(mr.output@BestModel_Estimate)),
		signif.log10po = rep(as.double(NA), len=length(mr.output@BestModel_Estimate)),
		time.secs = as.double(difftime(check_time, start_time, units="secs"))
	)
	names(ret.map@estimate) <- ret.map@names
	names(ret.map@stderror) <- ret.map@names
	names(ret.map@signif.neglog10padj) <- ret.map@names
	names(ret.map@signif.log10po) <- ret.map@names
	validate(ret.map)

	return(list('bma'=ret, 'map'=ret.map, 'MRBF2'=mr.output))
}

# Post model selection, fit a single model with all add1 and drop1 p-values
add1drop1 = function(data, binary.inclusion.vector, nu=data@m, print.ssq=FALSE, e.value.kappa=0.1, analysis.name="Leave one out/add one in significance testing", params=NULL) {
	if(!is.null(params)) {
		stopifnot(length(params)==data@m)
		stopifnot(all(!is.na(params)))
	}
	start_time = Sys.time()
	dataf = data.frame(y=data@y, data@x)
	model.formula = paste(data@y.name, "~ 0")
	if(sum(binary.inclusion.vector)>0) {
		model.formula = paste(model.formula, paste(data@x.names[binary.inclusion.vector], collapse=" + "), sep=" + ")
	}
	fit = lm(as.formula(model.formula), data=dataf)
	fit.coef = summary(fit)$coefficients
	if(print.ssq) cat("sigma^2 =", summary(fit)$sigma^2, "\n")

	estimate = rep(as.double(0), data@m)
	estimate[binary.inclusion.vector] = fit.coef[,"Estimate"]
	stderror = rep(as.double(0), data@m)
	stderror[binary.inclusion.vector] = fit.coef[,"Std. Error"]

	pval = rep(as.double(NA), len=data@m)
	if(sum(binary.inclusion.vector)>0) {
		drop.p = as.matrix(drop1(fit, test="Chisq"))[,"Pr(>Chi)"][-1]
		pval[match(names(drop.p), data@x.names)] = drop.p
	}
	if(sum(!binary.inclusion.vector)>0) {
		other.vars = paste("~.", paste(data@x.names[!binary.inclusion.vector], collapse=" + "), sep=" + ")
		add.p = as.matrix(add1(fit, scope = formula(other.vars), test="Chisq"))[,"Pr(>Chi)"][-1]
		pval[match(names(add.p), data@x.names)] = add.p
	}
	
	ret = new("bmsim_analysisResults",
		analysis = analysis.name,
		data = data@id,
		m = data@m,
		names = data@x.names,
		estimate = estimate,
		stderror = stderror,
		signif.neglog10padj = -log10(nu * pval),
		signif.log10po = rep(as.double(NA), len=data@m),
		time.secs = as.double(difftime(Sys.time(), start_time, units="secs"))
	)
	names(ret@estimate) <- ret@names
	names(ret@stderror) <- ret@names
	names(ret@signif.neglog10padj) <- ret@names
	names(ret@signif.log10po) <- ret@names
	names(pval) <- ret@names
	ret@pvalueTests = list(
		"Bonf" = calc.Bonferroni(pval, nu, params),
		"BH" = calc.BH(pval, nu, params),
		"HMP" = calc.HMP(pval, nu, params),
		"Simes" = calc.Simes(pval, nu, params),
		"Hommel" = calc.Hommel(pval, nu, params),
		"Cauchy" = calc.Cauchy(pval, nu, params),
		"Evalue" = calc.evalue(pval, nu, kappa=e.value.kappa, params)
	)
	validate(ret)
	return(ret)
}

# Exhaustive doublethink analysis
# h and mu can be vectors
# By default, unit information prior (h=1) and 10% expected inclusion probability
# By default, nu is the number of variables in the data, but can be set bigger
doublethink.x = function(data, h=1, mu=.1/(1-.1), nu=data@m, e.value.kappa=0.1, params=NULL) {
	validate(data)
	start_time = Sys.time()
	# Process h and mu
	nanal = max(length(h), length(mu))
	if(nanal>1) {
		if(length(h)==1) h = rep(h, nanal)
		if(length(mu)==1) mu = rep(mu, nanal)
	}
	stopifnot(length(h)==nanal)
	stopifnot(length(mu)==nanal)
	hyper.names = paste0("mu = ", mu, "; h = ", h)
	if(!is.null(params)) {
		stopifnot(length(params)==data@m)
		stopifnot(all(!is.na(params)))
	}
	# Power parameter
	xi = h/(data@n + h)
	# Validate nu
	stopifnot(nu>=data@m)
	# Total number of additive models
	L = as.integer(2^data@m)
	# Model inclusion matrix
	s = t(sapply(0:(L-1), function(s) {
		(incl = intToBits(s)[data@m:1]==1)
	})); colnames(s) <- data@x.names
	# Temporary storage
	loglik = rep(as.double(NA), L)
	degfree = rep(as.double(NA), L)
	estimate = matrix(as.double(0), L, data@m)
	stderror = matrix(as.double(0), L, data@m)
	# Fit the models (computationally intensive step)
	dataf = data.frame(y=data@y, data@x)
	for(i in 1:L) {
		incl = c(TRUE, s[i,])
		model.formula = paste(data@y.name, "~ 0 + .")
		fit = lm(as.formula(model.formula), data=dataf[, incl, drop=FALSE])
		fit.coef = summary(fit)$coefficients
		loglik[i] = logLik(fit)
		degfree[i] = attributes(logLik(fit))$df
		estimate[i,incl[-1]] = fit.coef[,"Estimate"]
		stderror[i,incl[-1]] = fit.coef[,"Std. Error"]
	}
	# Re-centre all log-likelihoods and degrees of freedom against the grand null
	loglik = loglik-loglik[1]
	degfree = degfree-degfree[1]
	# Posterior odds and posterior probability
	logc = log(mu[1] * sqrt(xi[1]))
	logPO = degfree * logc + (1-xi[1]) * loglik
	PO.prop = exp(logPO - max(logPO))
	PP = matrix(PO.prop/sum(PO.prop), ncol=1)
	if(nanal>1) {
		for(j in 2:nanal) {
			logc = log(mu[j] * sqrt(xi[j]))
			logPO = degfree * logc + (1-xi[j]) * loglik
			PO.prop = exp(logPO - max(logPO))
			PP = cbind(PP, PO.prop/sum(PO.prop))
		}
	}
	colnames(PP) <- hyper.names
	end_time = Sys.time()
	# Model inclusion vector for pairwise tests
	s.pairwise = matrix(as.logical(NA), nrow=nrow(s), ncol=0)
	for(v1 in 1:(data@m-1)) {
		for(v2 in (v1+1):data@m) {
			s.pairwise = cbind(s.pairwise, s[,v1] | s[,v2])
			colnames(s.pairwise)[ncol(s.pairwise)] <- paste(colnames(s)[c(v1, v2)], collapse=" | ")
		}
	}
	# Model inclusion vector for one-degree of freedom pairwise tests
	s.pairwise.xor = matrix(as.logical(NA), nrow=nrow(s), ncol=0)
	for(v1 in 1:(data@m-1)) {
		for(v2 in (v1+1):data@m) {
			s.pairwise.xor = cbind(s.pairwise.xor, xor(s[,v1], s[,v2]))
			colnames(s.pairwise.xor)[ncol(s.pairwise.xor)] <- paste(colnames(s)[c(v1, v2)], collapse=" ^ ")
		}
	}
	# Model inclusion vector for set of true nulls
	s.truenull = matrix(as.logical(NA), nrow=nrow(s), ncol=1)
	if(!is.null(params)) {
		wh = which(params==0)
		if(length(wh)>0) {
			s.truenull = matrix(rowSums(s[,wh,drop=FALSE])>0, nrow=nrow(s), ncol=1)
			colnames(s.truenull) <- paste(colnames(s)[wh], collapse=" | ")
		}
	}
	# Model inclusion vector for set of true nulls: one degree of freedom tests
	s.truenull.xor = matrix(as.logical(NA), nrow=nrow(s), ncol=1)
	if(!is.null(params)) {
		wh = which(params==0)
		if(length(wh)>0) {
			s.truenull.xor = matrix(rowSums(s[,wh,drop=FALSE])==1, nrow=nrow(s), ncol=1)
			colnames(s.truenull.xor) <- paste(colnames(s)[wh], collapse=" ^ ")
		}
	}
	# Compute the return objects: model-averaged results
	doublethink.bma = list()
	doublethink.bma.preciser = list()
	doublethink.bma.1df.preciser = list()
	for(j in 1:nanal) {
		PO.marginal = colSums(PP[,j]*s)/colSums(PP[,j]*(1-s))
		PO.pairwise = colSums(PP[,j]*s.pairwise)/colSums(PP[,j]*(1-s.pairwise))
		PO.truenull = colSums(PP[,j]*s.truenull)/colSums(PP[,j]*(1-s.truenull))
		PO.headline = sum(PP[-1,j])/PP[1,j]
		names(PO.headline) <- paste(colnames(s), collapse = " | ")
		# One degree of freedom tests
		PO.pairwise.1df = colSums(PP[,j]*s.pairwise.xor)/colSums(PP[,j]*(1-s.pairwise)) # Note xor in numerator; or in the denominator
		PO.truenull.1df = colSums(PP[,j]*s.truenull.xor)/colSums(PP[,j]*(1-s.truenull)) # Note xor in numerator; or in the denominator
		PO.headline.1df = sum(PP[,j]*(degfree==1))/PP[1,j]
		names(PO.headline.1df) <- paste(colnames(s), collapse = " ^ ")
		# Closed testing procedure p-values under Theorem 2
		p.adj.marginal = pchisq(2*log(PO.marginal/nu/mu[j]/sqrt(xi[j])), 1, low=FALSE); p.adj.marginal[p.adj.marginal>0.02] = 1
		p.adj.pairwise = pchisq(2*log(PO.pairwise/nu/mu[j]/sqrt(xi[j])), 1, low=FALSE); p.adj.pairwise[p.adj.pairwise>0.02] = 1
		p.adj.truenull = pchisq(2*log(PO.truenull/nu/mu[j]/sqrt(xi[j])), 1, low=FALSE); p.adj.truenull[p.adj.truenull>0.02] = 1
		p.adj.headline = pchisq(2*log(PO.headline/nu/mu[j]/sqrt(xi[j])), 1, low=FALSE); p.adj.headline[p.adj.headline>0.02] = 1
		# Unadjusted p-values under Theorem 1
		p.unadj.marginal = pchisq(2*log(PO.marginal/mu[j]/sqrt(xi[j])), 1, low=FALSE); p.unadj.marginal[p.unadj.marginal>0.02] = 1
		names(p.unadj.marginal) <- data@x.names
		doublethink.bma[[j]] =
			new("bmsim_analysisResults",
				analysis = "Bayesian model averaged multivariable Mendelian randomization with Doublethink",
				 data = data@id,
				 m = data@m,
				 names = data@x.names,
				 estimate = colSums(PP[,j]*estimate),
				 stderror = sqrt(colSums(PP[,j]*(stderror^2 + estimate^2)) - colSums(PP[,j]*estimate)^2),
				 signif.neglog10padj = -log10(p.adj.marginal),
				 signif.log10po = log10(PO.marginal),
				 pairwise.signif.neglog10padj = -log10(p.adj.pairwise),
				 pairwise.signif.log10po = log10(PO.pairwise),
				 truenull.signif.neglog10padj = -log10(p.adj.truenull),
				 truenull.signif.log10po = log10(PO.truenull),
				 headline.signif.neglog10padj = -log10(p.adj.headline),
				 headline.signif.log10po = log10(PO.headline),
				 pvalueTests = list(
					"Bonf" = calc.Bonferroni(p.unadj.marginal, nu, params=params),
					"BH" = calc.BH(p.unadj.marginal, nu, params=params),
					"HMP" = calc.HMP(p.unadj.marginal, nu, params=params),
					"Simes" = calc.Simes(p.unadj.marginal, nu, params=params),
					"Hommel" = calc.Hommel(p.unadj.marginal, nu, params=params),
					"Cauchy" = calc.Cauchy(p.unadj.marginal, nu, params=params),
					"Evalue" = calc.evalue(p.unadj.marginal, nu, kappa=e.value.kappa, params=params),
					"Evalue.BF2p" = calc.evalue.BF2p(log10(PO.marginal) - log10(mu[j]), log10(PO.pairwise) - log10(2*mu[j]), log10(PO.truenull) - log10(sum(params==0)*mu[j]), log10(PO.headline) - log10(nu*mu[j]), nu)
				 ),
				 time.secs = as.double(difftime(end_time, start_time, units="secs"))
			)
		# One degree of freedom tests only (breaks the CTP except asymptotically)
		p.adj.pairwise.1df = pchisq(2*log(PO.pairwise.1df/nu/mu[j]/sqrt(xi[j])), 1, low=FALSE); p.adj.pairwise.1df[p.adj.pairwise.1df>0.02] = 1
		p.adj.truenull.1df = pchisq(2*log(PO.truenull.1df/nu/mu[j]/sqrt(xi[j])), 1, low=FALSE); p.adj.truenull.1df[p.adj.truenull.1df>0.02] = 1
		p.adj.headline.1df = pchisq(2*log(PO.headline.1df/nu/mu[j]/sqrt(xi[j])), 1, low=FALSE); p.adj.headline.1df[p.adj.headline.1df>0.02] = 1
		doublethink.bma.1df.preciser[[j]] =
			new("bmsim_analysisResults",
				analysis = "Bayesian model averaged multivariable Mendelian randomization with Doublethink (1df tests preciser)",
				 data = data@id,
				 m = data@m,
				 names = data@x.names,
				 estimate = colSums(PP[,j]*estimate),
				 stderror = sqrt(colSums(PP[,j]*(stderror^2 + estimate^2)) - colSums(PP[,j]*estimate)^2),
				 signif.neglog10padj = -log10(p.adj.marginal),
				 signif.log10po = log10(PO.marginal),
				 pairwise.signif.neglog10padj = -log10(p.adj.pairwise.1df),
				 pairwise.signif.log10po = log10(PO.pairwise.1df),
				 truenull.signif.neglog10padj = -log10(p.adj.truenull.1df),
				 truenull.signif.log10po = log10(PO.truenull.1df),
				 headline.signif.neglog10padj = -log10(p.adj.headline.1df),
				 headline.signif.log10po = log10(PO.headline.1df),
				 pvalueTests = list(
					"Bonf" = calc.Bonferroni(p.unadj.marginal, nu, params=params),
					"BH" = calc.BH(p.unadj.marginal, nu, params=params),
					"HMP" = calc.HMP(p.unadj.marginal, nu, params=params),
					"Simes" = calc.Simes(p.unadj.marginal, nu, params=params),
					"Hommel" = calc.Hommel(p.unadj.marginal, nu, params=params),
					"Cauchy" = calc.Cauchy(p.unadj.marginal, nu, params=params),
					"Evalue" = calc.evalue(p.unadj.marginal, nu, kappa=e.value.kappa, params=params),
					"Evalue.BF2p" = calc.evalue.BF2p(log10(PO.marginal) - log10(mu[j]), log10(PO.pairwise.1df) - log10(2*mu[j]), log10(PO.truenull.1df) - log10(sum(params==0)*mu[j]), log10(PO.headline.1df) - log10(nu*mu[j]), nu)
				 ),
				 time.secs = as.double(difftime(end_time, start_time, units="secs"))
			)		# Closed testing procedure p-values under Theorem 2
		# More precise version of the denominators for when n is not huge (e.g. 25% multiplicative difference for n=145, mu=0.2, h=4)
		denom = (1 + mu[j]*sqrt(xi[j]))^nu - 1
		p.adj.marginal = pchisq(2*log(PO.marginal/denom), 1, low=FALSE); p.adj.marginal[p.adj.marginal>0.02] = 1
		p.adj.pairwise = pchisq(2*log(PO.pairwise/denom), 1, low=FALSE); p.adj.pairwise[p.adj.pairwise>0.02] = 1
		p.adj.truenull = pchisq(2*log(PO.truenull/denom), 1, low=FALSE); p.adj.truenull[p.adj.truenull>0.02] = 1
		p.adj.headline = pchisq(2*log(PO.headline/denom), 1, low=FALSE); p.adj.headline[p.adj.headline>0.02] = 1
		# Unadjusted p-values under Theorem 1
		# No adjustment needed in the case that |V|=1
		p.unadj.marginal = pchisq(2*log(PO.marginal/mu[j]/sqrt(xi[j])), 1, low=FALSE); p.unadj.marginal[p.unadj.marginal>0.02] = 1
		names(p.unadj.marginal) <- data@x.names
		doublethink.bma.preciser[[j]] =
			new("bmsim_analysisResults",
				analysis = "Bayesian model averaged multivariable Mendelian randomization with Doublethink (preciser denominators)",
				 data = data@id,
				 m = data@m,
				 names = data@x.names,
				 estimate = colSums(PP[,j]*estimate),
				 stderror = sqrt(colSums(PP[,j]*(stderror^2 + estimate^2)) - colSums(PP[,j]*estimate)^2),
				 signif.neglog10padj = -log10(p.adj.marginal),
				 signif.log10po = log10(PO.marginal),
				 pairwise.signif.neglog10padj = -log10(p.adj.pairwise),
				 pairwise.signif.log10po = log10(PO.pairwise),
				 truenull.signif.neglog10padj = -log10(p.adj.truenull),
				 truenull.signif.log10po = log10(PO.truenull),
				 headline.signif.neglog10padj = -log10(p.adj.headline),
				 headline.signif.log10po = log10(PO.headline),
				 pvalueTests = list(
					"Bonf" = calc.Bonferroni(p.unadj.marginal, nu, params=params),
					"BH" = calc.BH(p.unadj.marginal, nu, params=params),
					"HMP" = calc.HMP(p.unadj.marginal, nu, params=params),
					"Simes" = calc.Simes(p.unadj.marginal, nu, params=params),
					"Hommel" = calc.Hommel(p.unadj.marginal, nu, params=params),
					"Cauchy" = calc.Cauchy(p.unadj.marginal, nu, params=params),
					"Evalue" = calc.evalue(p.unadj.marginal, nu, kappa=e.value.kappa, params=params),
					"Evalue.BF2p" = calc.evalue.BF2p(log10(PO.marginal) - log10(mu[j]), log10(PO.pairwise) - log10(2*mu[j]), log10(PO.truenull) - log10(sum(params==0)*mu[j]), log10(PO.headline) - log10(nu*mu[j]), nu)
				 ),
				 time.secs = as.double(difftime(end_time, start_time, units="secs"))
			)
	}
	names(doublethink.bma) <- hyper.names
	names(doublethink.bma.preciser) <- hyper.names
	names(doublethink.bma.1df.preciser) <- hyper.names
	# Compute the return objects: model selection results
	doublethink.modelselection = list()
	for(j in 1:nanal) {
		gd = which.max(PP[, j])
		doublethink.modelselection[[j]] = add1drop1(data, s[gd,], nu, analysis.name="Bayesian model selection multivariable Mendelian randomization with Doublethink; leave one out/add one in significance testing", params=params)
	}
	names(doublethink.modelselection) <- hyper.names

	return(list(
		doublethink.internal = list(
			mu = mu,
			h = h,
			s = s,
			loglik = loglik,
			degfree = degfree,
			estimate = estimate,
			stderror = stderror
		),
		doublethink.bma = doublethink.bma,
		doublethink.bma.preciser = doublethink.bma.preciser,
		doublethink.bma.1df.preciser = doublethink.bma.1df.preciser,
		doublethink.modelselection = doublethink.modelselection
   ))
}

# Stepwise regression. By default, BIC criterion (k = log(data@n))
do.stepwise = function(data, k=log(data@n), nu=data@m, params=NULL) {
	dataf = data.frame(y=data@y, data@x)
	model.formula = as.formula(paste(data@y.name, "~ 0 + ."))
	sw = step(lm(model.formula, data=dataf), direction="backward", k=k, trace=0)
	binary.inclusion.vector = !is.na(match(data@x.names, names(coef(sw))))
	add1drop1(data, binary.inclusion.vector, nu, analysis.name="Multivariable Mendelian randomization with stepwise regression; leave one out/add one in significance testing", params=params)
}

# LASSO, elastic net and ridge regression estimates
do.glmnet = function(data, grid.len=11, nu=data@m) {
	elastic.alpha = seq(0, 1, len=grid.len)
	
	# LASSO
	start_time = Sys.time()
	ret.cv.glmnet = cv.glmnet(data@x, data@y, keep=TRUE, intercept=FALSE)
	lasso.estimate = coef(ret.cv.glmnet)[-1]
	
	lasso.ret = new("bmsim_analysisResults",
		analysis = "LASSO",
		data = data@id,
		m = data@m,
		names = data@x.names,
		estimate = lasso.estimate,
		stderror = rep(as.double(NA), data@m),
		signif.neglog10padj = rep(as.double(NA), data@m),
		signif.log10po = rep(as.double(NA), len=data@m),
		time.secs = as.double(difftime(Sys.time(), start_time, units="secs"))
	)
	names(lasso.ret@estimate) <- lasso.ret@names
	names(lasso.ret@stderror) <- lasso.ret@names
	names(lasso.ret@signif.neglog10padj) <- lasso.ret@names
	names(lasso.ret@signif.log10po) <- lasso.ret@names
	validate(lasso.ret)
	
	# Elastic net
	elastic.ret.cv.glmnet = lapply(elastic.alpha, function(ALPHA) cv.glmnet(data@x, data@y, intercept=FALSE, alpha=ALPHA, foldid=ret.cv.glmnet$foldid))
	# Compare the cross-validation error across elastic.alpha
	elastic.cvm = unlist(lapply(elastic.ret.cv.glmnet, function(fit) min(fit$cvm)))
	best.elastic = which.min(elastic.cvm)
	elnet.estimate = coef(elastic.ret.cv.glmnet[[best.elastic]])[-1]

	elnet.ret = new("bmsim_analysisResults",
		analysis = "Elastic net",
		data = data@id,
		m = data@m,
		names = data@x.names,
		estimate = elnet.estimate,
		stderror = rep(as.double(NA), data@m),
		signif.neglog10padj = rep(as.double(NA), data@m),
		signif.log10po = rep(as.double(NA), len=data@m),
		time.secs = as.double(difftime(Sys.time(), start_time, units="secs"))
	)
	names(elnet.ret@estimate) <- elnet.ret@names
	names(elnet.ret@stderror) <- elnet.ret@names
	names(elnet.ret@signif.neglog10padj) <- elnet.ret@names
	names(elnet.ret@signif.log10po) <- elnet.ret@names
	validate(elnet.ret)

	# Ridge regression
	check_time = Sys.time()
	ridge.estimate = coef(elastic.ret.cv.glmnet[[1]])[-1]
	
	ridge.ret = new("bmsim_analysisResults",
		analysis = "Ridge regression",
		data = data@id,
		m = data@m,
		names = data@x.names,
		estimate = ridge.estimate,
		stderror = rep(as.double(NA), data@m),
		signif.neglog10padj = rep(as.double(NA), data@m),
		signif.log10po = rep(as.double(NA), len=data@m),
		time.secs = as.double(difftime(Sys.time(), check_time, units="secs"))
	)
	names(ridge.ret@estimate) <- ridge.ret@names
	names(ridge.ret@stderror) <- ridge.ret@names
	names(ridge.ret@signif.neglog10padj) <- ridge.ret@names
	names(ridge.ret@signif.log10po) <- ridge.ret@names
	validate(ridge.ret)

	return(list('lasso'=lasso.ret, 'elnet'=elnet.ret, 'ridge'=ridge.ret))
}

# Generate a function to perform simulations with sample size up to max.n with fixed independent variables
# Model is for Mendelian randomization with the specified causal effects
gen.simulate = function(data, max.n=1000, seed=NA, independence=FALSE) {
	stopifnot(is(data, "bmsim_data"))
	validate(data)
	x.mu = colMeans(data@x)
	x.sigma = cov(data@x)
	if(independence) x.sigma = diag(diag(x.sigma))
	# Simulate the independent variables just once: do not wish to study randomness in the
	# data generating process of x, but rather condition on x
	if(!is.na(seed)) set.seed(seed)
	x = rmvnorm(max.n, x.mu, x.sigma)
	simulate = function(n, param) {
		n = as.integer(n)
		if(!(data@m<n & n<=max.n)) stop(paste0(data@m, " < n <= ", max.n, " is not TRUE"))
		if(!(length(param)==data@m)) stop(paste0("length(param) == ", data@m, " is not TRUE"))
		ret = new("bmsim_data",
			id = paste0("simulation of ", data@id, " @ ", Sys.time()),
			y = matrix(rnorm(n, as.vector(x[1:n,] %*% matrix(param, ncol=1)), 1), ncol=1),
			x = x[1:n,],
			n = n,
			m = data@m,
			y.name = data@y.name,
			x.names = data@x.names
		)
		colnames(ret@y) <- ret@y.name
		colnames(ret@x) <- ret@x.names
		return(ret)
	}
	return(simulate)
}

# Combine performance for a specific test
rowMeans.robust = function(x) {
	if(is.null(dim(x))) x = matrix(x, nrow=1)
	rowMeans(x)
}
combine.test.performance = function(performance.list, test.name) {
	nsim = length(performance.list)
	ret = new("bmsim_test_performance")
	if(length(slot(performance.list[[1]], test.name)@bias)>0) {
		ret = new("bmsim_test_performance",
			bias = rowMeans.robust(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@bias)),
			typeI = rowMeans.robust(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@typeI)),
			typeII = rowMeans.robust(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@typeII)),
			call = rowMeans.robust(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@call)),
			familywiseI = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@familywiseI), na.rm=TRUE),
			familywiseII = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@familywiseII), na.rm=TRUE),
			pfdr = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@pfdr), na.rm=TRUE),
			fdr = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@fdr)),
			fndr = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@fndr)),
			fpr = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@fpr), na.rm=TRUE),
			fnr = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@fnr), na.rm=TRUE),
			tp = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@tp)),
			fp = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@fp)),
			tn = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@tn)),
			fn = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)@fn))
		)
		names(ret@bias) <- names(slot(performance.list[[1]], test.name)@bias)
		names(ret@typeI) <- names(slot(performance.list[[1]], test.name)@typeI)
		names(ret@typeII) <- names(slot(performance.list[[1]], test.name)@typeII)
		names(ret@call) <- names(slot(performance.list[[1]], test.name)@call)
	}
	return(ret)
}
combine.pvalueTest.performance = function(performance.list, test.name, pvalueTest.name) {
	# test.name is one of pvalueTests.marginal pvalueTests.pairwise pvalueTests.truenull pvalueTests.headline
	nsim = length(performance.list)
	ret = new("bmsim_test_performance")
	if(length(slot(performance.list[[1]], test.name)[[pvalueTest.name]]@bias)>0) {
		ret = new("bmsim_test_performance",
			bias = rowMeans.robust(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@bias)),
			typeI = rowMeans.robust(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@typeI)),
			typeII = rowMeans.robust(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@typeII)),
			call = rowMeans.robust(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@call)),
			familywiseI = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@familywiseI), na.rm=TRUE),
			familywiseII = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@familywiseII), na.rm=TRUE),
			pfdr = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@pfdr), na.rm=TRUE),
			fdr = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@fdr)),
			fndr = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@fndr)),
			fpr = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@fpr), na.rm=TRUE),
			fnr = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@fnr), na.rm=TRUE),
			tp = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@tp)),
			fp = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@fp)),
			tn = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@tn)),
			fn = mean(sapply(1:nsim, function(i) slot(performance.list[[i]], test.name)[[pvalueTest.name]]@fn))
		)
		names(ret@bias) <- names(slot(performance.list[[1]], test.name)[[pvalueTest.name]]@bias)
		names(ret@typeI) <- names(slot(performance.list[[1]], test.name)[[pvalueTest.name]]@typeI)
		names(ret@typeII) <- names(slot(performance.list[[1]], test.name)[[pvalueTest.name]]@typeII)
		names(ret@call) <- names(slot(performance.list[[1]], test.name)[[pvalueTest.name]]@call)
	}
	return(ret)
}

# Combine a list of performance objects applied to comparable analyses
combine.performance = function(performance.list) {
	nsim = length(performance.list)
	for(i in 1:nsim) {
		stopifnot(is(performance.list[[i]], "bmsim_performance"))
		stopifnot(performance.list[[i]]@m == performance.list[[1]]@m)
		stopifnot(all(performance.list[[i]]@names == performance.list[[1]]@names))
	}
	
	ret = new("bmsim_performance",
		analysis = unique(sapply(1:nsim, function(i) performance.list[[i]]@analysis)),
		data = unique(sapply(1:nsim, function(i) performance.list[[i]]@data)),
		m = performance.list[[1]]@m,
		names =  performance.list[[1]]@names
	)

	# Estimate
	ret@estimate.bias = rowMeans(sapply(1:nsim, function(i) performance.list[[i]]@estimate.bias))
	ret@estimate.L1 = mean(sapply(1:nsim, function(i) performance.list[[i]]@estimate.L1))
	ret@estimate.L2 = mean(sapply(1:nsim, function(i) performance.list[[i]]@estimate.L2))
	
	# Std Error
	ret@stderr.coverage = rowMeans(sapply(1:nsim, function(i) performance.list[[i]]@stderr.coverage))
	
	# Hypothesis tests
	ret@freqt.marginal = combine.test.performance(performance.list, "freqt.marginal")
	ret@bayes.marginal = combine.test.performance(performance.list, "bayes.marginal")
	if(length(performance.list[[1]]@freqt.pairwise@bias)>0) ret@freqt.pairwise = combine.test.performance(performance.list, "freqt.pairwise")
	if(length(performance.list[[1]]@bayes.pairwise@bias)>0) ret@bayes.pairwise = combine.test.performance(performance.list, "bayes.pairwise")
	if(length(performance.list[[1]]@freqt.truenull@bias)>0) ret@freqt.headline = combine.test.performance(performance.list, "freqt.truenull")
	if(length(performance.list[[1]]@bayes.truenull@bias)>0) ret@bayes.headline = combine.test.performance(performance.list, "bayes.truenull")
	if(length(performance.list[[1]]@freqt.headline@bias)>0) ret@freqt.headline = combine.test.performance(performance.list, "freqt.headline")
	if(length(performance.list[[1]]@bayes.headline@bias)>0) ret@bayes.headline = combine.test.performance(performance.list, "bayes.headline")

	# pvalue tests
	ret@pvalueTests.marginal = lapply(names(performance.list[[1]]@pvalueTests.marginal), function(pvalueTest.name) combine.pvalueTest.performance(performance.list, "pvalueTests.marginal", pvalueTest.name))
	names(ret@pvalueTests.marginal) <- names(performance.list[[1]]@pvalueTests.marginal)
	ret@pvalueTests.pairwise = lapply(names(performance.list[[1]]@pvalueTests.pairwise), function(pvalueTest.name) combine.pvalueTest.performance(performance.list, "pvalueTests.pairwise", pvalueTest.name))
	names(ret@pvalueTests.pairwise) <- names(performance.list[[1]]@pvalueTests.pairwise)
	ret@pvalueTests.truenull = lapply(names(performance.list[[1]]@pvalueTests.truenull), function(pvalueTest.name) combine.pvalueTest.performance(performance.list, "pvalueTests.truenull", pvalueTest.name))
	names(ret@pvalueTests.truenull) <- names(performance.list[[1]]@pvalueTests.truenull)
	ret@pvalueTests.headline = lapply(names(performance.list[[1]]@pvalueTests.headline), function(pvalueTest.name) combine.pvalueTest.performance(performance.list, "pvalueTests.headline", pvalueTest.name))
	names(ret@pvalueTests.headline) <- names(performance.list[[1]]@pvalueTests.headline)

	# Prediction
	ret@prediction.L1 = mean(sapply(1:nsim, function(i) performance.list[[i]]@prediction.L1))
	ret@prediction.L2 = mean(sapply(1:nsim, function(i) performance.list[[i]]@prediction.L2))

	# Computation time (seconds)
	ret@time.secs = mean(sapply(1:nsim, function(i) performance.list[[i]]@time.secs))

	# Apply names
	names(ret@estimate.bias) <- ret@names
	names(ret@stderr.coverage) <- ret@names
	names(ret@freqt.marginal@bias) <- ret@names
	names(ret@freqt.marginal@typeI) <- ret@names
	names(ret@freqt.marginal@typeII) <- ret@names
	names(ret@bayes.marginal@bias) <- ret@names
	names(ret@bayes.marginal@typeI) <- ret@names
	names(ret@bayes.marginal@typeII) <- ret@names

	return(ret)
}

# Perform a standardized set of analyses of a single dataset
do.analyses = function(data, params, nu=data@m, dblthk.h = c(0.25, 1, 4), dblthk.mu = c(0.05, 0.1, 0.2), mr.bma.nsim=1000, mr.bma.sigma=0.5, mr.bma.prior_prob=0.1) {
	stopifnot(is(data, "bmsim_data"))
	validate(data)
	stopifnot(length(params)==data@m)
	stopifnot(!any(is.na(params)))
	
	results = list()
	# Single model results based on the 'oracle' model
	results[["oracle"]] = add1drop1(data, params!=0, nu, analysis.name="Multivariable Mendelian randomization with oracle model; leave one out/add one in significance testing", params=params)

	# Single model results based on the grand null
	results[["grand null"]] = add1drop1(data, rep(FALSE, data@m), nu, analysis.name="Multivariable Mendelian randomization with no variables; leave one out/add one in significance testing", params=params)

	# Single model results based on the grand alternative
	results[["grand alternative"]] = add1drop1(data, rep(TRUE, data@m), nu, analysis.name="Multivariable Mendelian randomization with all variables; leave one out/add one in significance testing", params=params)

	# Stepwise regression
	results[["backward elimination"]] = do.stepwise(data, nu=nu, params=params)
	
	# Fit LASSO, elastic net and ridge regression
	results.glmnet = do.glmnet(data, nu=nu)
	results[["lasso"]] = results.glmnet[["lasso"]]
	results[["elastic net"]] = results.glmnet[["elnet"]]
	results[["ridge regression"]] = results.glmnet[["ridge"]]
	# LASSO-based model selection
	results[["lasso model selection"]] = add1drop1(data, results.glmnet$lasso@estimate!=0, nu, analysis.name="Multivariable Mendelian randomization with LASSO; leave one out/add one in significance testing", params=params)
	# Elastic net-based model selection
	results[["elastic net model selection"]] = add1drop1(data, results.glmnet$elnet@estimate!=0, nu, analysis.name="Multivariable Mendelian randomization with elastic net; leave one out/add one in significance testing", params=params)

	# Perform Bayesian model-averaged Mendelian randomization using Doublethink
	vec.dblthk.h = rep(dblthk.h, each=length(dblthk.mu))
	vec.dblthk.mu = rep(dblthk.mu, length(dblthk.h))
	results.doublethink = doublethink.x(data, h=vec.dblthk.h, mu=vec.dblthk.mu, nu=nu, params=params)
	#doublethink.names = outer(dblthk.h, dblthk.mu, function(h,mu) paste0("mu = ", mu, "; h = ", h))
	doublethink.names = paste0("mu = ", vec.dblthk.mu, "; h = ", vec.dblthk.h)
	for(doublethink.name in doublethink.names) {
		# Doublethink BMA
		results[[paste0("doublethink bma ", doublethink.name)]] = results.doublethink$doublethink.bma[[doublethink.name]]
		# Doublethink BMA: preciser denominators
		results[[paste0("doublethink bma preciser ", doublethink.name)]] = results.doublethink$doublethink.bma.preciser[[doublethink.name]]
		# Doublethink BMA: 1df preciser denominators
		results[[paste0("doublethink bma 1df preciser ", doublethink.name)]] = results.doublethink$doublethink.bma.1df.preciser[[doublethink.name]]
		# Doublethink MAP-based model selection
		results[[paste0("doublethink model selection ", doublethink.name)]] = results.doublethink$doublethink.modelselection[[doublethink.name]]
	}

	# Perform Bayesian model-averaged Mendelian randomization using MR-BMA (no permutation procedure)
	results.mr.bma = mr.bma.x(data, sigma=mr.bma.sigma, prior_prob=mr.bma.prior_prob, nsim=0, nu=nu)
	results[["mr-bma bma"]] = results.mr.bma$bma
	# results[["mr-bma map"]] = results.mr.bma$map
	# MR-BMA MAP-based model selection
	results[["mr-bma model selection"]] = add1drop1(data, results.mr.bma$map@estimate!=0, nu, analysis.name="Bayesian model selection multivariable Mendelian randomization with MR-BMA; leave one out/add one in significance testing", params=params)

	# Perform Bayesian model-averaged Mendelian randomization using MR-BMA (with permutation procedure)
	# (Do this separately for timing purposes)
	results.mr.bma = mr.bma.x(data, sigma=mr.bma.sigma, prior_prob=mr.bma.prior_prob, nsim=mr.bma.nsim, nu=nu)
	results[["mr-bma bma permute"]] = results.mr.bma$bma
	
	return(results)
}

# Perform a standardized set of analyses of a single dataset
calc.performance = function(analyses, params, freqt.alpha = 0.05, bayes.tau = 19, newdata=NULL) {
	stopifnot(is.list(analyses))
	nanal = length(analyses)
	for(analysis in analyses) {
		stopifnot(is(analysis, "bmsim_analysisResults"))
		validate(analysis)
		stopifnot(analysis@m==length(params))
	}
	if(!is.null(newdata)) {
		stopifnot(newdata@m==length(params))
		validate(newdata)
	}
	
	thresh.neglog10padj = -log10(freqt.alpha)
	thresh.log10po = log10(bayes.tau)
	ret = list()
	for(i in 1:nanal) {
		analysis.name = names(analyses)[i]
		ret[[analysis.name]] = performance(analyses[[analysis.name]], parameter=params, newdata=newdata, thresh.neglog10padj=thresh.neglog10padj, thresh.log10po=thresh.log10po)
	}
	return(ret)
}

# Populate parameters for use in the simulations. Idea is to execute this once and hard-code it into the git repository
# Set seed to NULL to avoid deterministic results
# simulate.parameters(3, 15, "prior", filename=stdout(), data=data, seed=NULL)
simulate.parameters = function(nsim=1000, m=15, model=c("grand_null", "bad_case", "prior")[1], filename=paste0(model, ".parameters.txt"), overwrite=FALSE, data=NULL, tau=19, h=1, mu=0.1, n=145, seed=0) {
	if(!is.null(data)) {
		stopifnot(data@m==m)
		validate(data)
	}
	if(!is.null(seed)) set.seed(seed)
	if(filename!=stdout()) stopifnot(!file.exists(filename) | overwrite==TRUE)
	ret = matrix(as.double(NA), nsim, m)
	model = tolower(model)
	if(model=="grand_null") {
		# All parameters are zero
		ret = matrix(as.double(0.0), nsim, m)
	} else if(model=="bad_case") {
		# Odd numbered parameters are zero
		# Even numbered parameters are non-zero
		# For each non-zero variable, find the highest absolute correlation rho and set beta = rho/(1-rho^2)/sqrt(n)*sqrt(2*log(tau/c))
		# Theory for the two-variable case suggests this generates inflation, albeit not the worst inflation (for which no expression is available)
		stopifnot(!is.null(data))
		c = mu*sqrt(h/(n+h))
		RHO = cor(data@x); diag(RHO) <- NA
		rho.max = RHO[cbind(1:m, apply(abs(RHO), 1, which.max))]
		params = ifelse((1:m) %% 2==0, rho.max/(1-rho.max^2)/sqrt(n)*sqrt(2*log(tau/c)), 0.0)
		ret = matrix(params, nsim, m, byrow=TRUE)
	} else if(model=="prior") {
		# Simulate from the Johnson (2005) prior given n, h, mu and data@x
		stopifnot(!is.null(data))
		unit.info = crossprod(data@x)/data@n
		prior.var = solve(unit.info)/h
		params = rmvnorm(nsim, rep(0.0, m), prior.var)
		params[sample(length(params), rbinom(1, length(params), 1/(1+mu)), replace=FALSE)] = 0.0
		ret = params
	} else error("simulate.parameters(): unrecognized model. Try one of 'grand_null', 'bad_case' or 'prior'")
	write(t(ret), file=filename, ncolumns=m)
}
