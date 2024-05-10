# Load required libraries
library(glmnet)
library(combinat)
library(mvtnorm)
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

# Class for standardized analysis results
setClass("bmsim_analysisResults",
	representation(
		analysis = "character",
		data = "character",
		m = "integer",
		names = "character",
		estimate = "numeric",
		stderror = "numeric",
		signif.neglog10padj = "numeric",
		signif.log10po = "numeric",
		time.secs = "numeric"
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
	pp2log10po = function(pp) log10(pp/(1-pp))
	ret = new("bmsim_analysisResults",
		analysis = "Bayesian model averaged multivariable Mendelian randomization with MR-BMA",
		data = data@id,
		m = data@m,
		names = data@x.names,
		estimate = mr.output@BMAve_Estimate,
		stderror = mr.output@BMAve_StdError,
		signif.neglog10padj = -log10(nu * pval),
		signif.log10po = pp2log10po(mr.output@pp_marginal),
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
add1drop1 = function(data, binary.inclusion.vector, nu=data@m, print.ssq=FALSE, analysis.name="Leave one out/add one in significance testing") {
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
	validate(ret)
	return(ret)
}

# Exhaustive doublethink analysis
# h and mu can be vectors
# By default, unit information prior (h=1) and 10% expected inclusion probability
# By default, nu is the number of variables in the data, but can be set bigger
doublethink.x = function(data, h=1, mu=.1/(1-.1), nu=data@m) {
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
	PO = exp(degfree * logc + (1-xi[1]) * loglik)
	PP = matrix(PO/sum(PO), ncol=1)
	if(nanal>1) {
		for(j in 2:nanal) {
			logc = log(mu[j] * sqrt(xi[j]))
			PO = exp(degfree * logc + (1-xi[j]) * loglik)
			PP = cbind(PP, PO/sum(PO))
		}
	}
	colnames(PP) <- hyper.names
	end_time = Sys.time()
	# Compute the return objects: model-averaged results
	doublethink.bma = list()
	for(j in 1:nanal) {
		PO = colSums(PP[,j]*s)/colSums(PP[,j]*(1-s))
		p.adj = pchisq(2*log(PO/nu/mu[j]/sqrt(xi[j])), 1, low=FALSE)
		p.adj[p.adj>0.02] = 1
		doublethink.bma[[j]] =
			new("bmsim_analysisResults",
				analysis = "Bayesian model averaged multivariable Mendelian randomization with Doublethink",
				 data = data@id,
				 m = data@m,
				 names = data@x.names,
				 estimate = colSums(PP[,j]*estimate),
				 stderror = sqrt(colSums(PP[,j]*(stderror^2 + estimate^2)) - colSums(PP[,j]*estimate)^2),
				 signif.neglog10padj = -log10(p.adj),
				 signif.log10po = log10(PO),
				 time.secs = as.double(difftime(end_time, start_time, units="secs"))
			)
	}
	names(doublethink.bma) <- hyper.names
	# Compute the return objects: model selection results
	doublethink.modelselection = list()
	for(j in 1:nanal) {
		gd = which.max(PP)
		doublethink.modelselection[[j]] = add1drop1(data, s[gd,], nu, analysis.name="Bayesian model selection multivariable Mendelian randomization with Doublethink; leave one out/add one in significance testing")
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
		doublethink.modelselection = doublethink.modelselection
   ))
}

# Stepwise regression. By default, BIC criterion (k = log(data@n))
do.stepwise = function(data, k=log(data@n), nu=data@m) {
	dataf = data.frame(y=data@y, data@x)
	model.formula = as.formula(paste(data@y.name, "~ 0 + ."))
	sw = step(lm(model.formula, data=dataf), direction="backward", k=k, trace=0)
	binary.inclusion.vector = !is.na(match(data@x.names, names(coef(sw))))
	add1drop1(data, binary.inclusion.vector, nu, analysis.name="Multivariable Mendelian randomization with stepwise regression; leave one out/add one in significance testing")
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
		freqt.bias = "numeric",
		freqt.typeI = "numeric",
		freqt.typeII = "numeric",
		freqt.familywiseI = "numeric",
		freqt.familywiseII = "numeric",
		freqt.fdr = "numeric",
		freqt.fndr = "numeric",
		freqt.fpr = "numeric",
		freqt.fnr = "numeric",
		# Bayesian test performance
		bayes.bias = "numeric",
		bayes.typeI = "numeric",
		bayes.typeII = "numeric",
		bayes.familywiseI = "numeric",
		bayes.familywiseII = "numeric",
		bayes.fdr = "numeric",
		bayes.fndr = "numeric",
		bayes.fpr = "numeric",
		bayes.fnr = "numeric",
		# Prediction performance
		prediction.L1 = "numeric",
		prediction.L2 = "numeric",
		# Computation time (seconds)
		time.secs = "numeric"
	)
)

setGeneric("performance", function(obj, ...) standardGeneric("performance"))
setMethod("performance", "bmsim_analysisResults", function(obj, parameter, newdata=NULL, thresh.neglog10padj=NA, thresh.log10po=NA) {
	validate(obj)
	stopifnot(length(parameter)==obj@m)
	if(!is.null(newdata)) stopifnot(newdata@m==obj@m)
	
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
	
	# Std Error
	ret@stderr.coverage = performance.stderr.coverage(obj@estimate, obj@stderror, parameter)
	
	# Hypothesis tests
	param.signif = 1*(parameter!=0)
	
	freqt.signif = 1*(obj@signif.neglog10padj >= thresh.neglog10padj)
	ret@freqt.bias = performance.bias(freqt.signif, param.signif)
	ret@freqt.typeI = performance.typeI(freqt.signif, param.signif)
	ret@freqt.typeII = performance.typeII(freqt.signif, param.signif)
	ret@freqt.familywiseI = performance.familywise(ret@freqt.typeI)
	ret@freqt.familywiseII = performance.familywise(ret@freqt.typeII)
	ret@freqt.fdr = performance.error.rate(ret@freqt.typeI, freqt.signif)
	ret@freqt.fndr = performance.error.rate(ret@freqt.typeII, 1-freqt.signif)
	ret@freqt.fpr = performance.error.rate(ret@freqt.typeI, 1-param.signif)
	ret@freqt.fnr = performance.error.rate(ret@freqt.typeII, param.signif)

	bayes.signif = 1*(obj@signif.log10po >= thresh.log10po)
	ret@bayes.bias = performance.bias(bayes.signif, param.signif)
	ret@bayes.typeI = performance.typeI(bayes.signif, param.signif)
	ret@bayes.typeII = performance.typeII(bayes.signif, param.signif)
	ret@bayes.familywiseI = performance.familywise(ret@bayes.typeI)
	ret@bayes.familywiseII = performance.familywise(ret@bayes.typeII)
	ret@bayes.fdr = performance.error.rate(ret@bayes.typeI, bayes.signif)
	ret@bayes.fndr = performance.error.rate(ret@bayes.typeII, 1-bayes.signif)
	ret@bayes.fpr = performance.error.rate(ret@bayes.typeI, 1-param.signif)
	ret@bayes.fnr = performance.error.rate(ret@bayes.typeII, param.signif)

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
	
	# Apply names
	names(ret@estimate.bias) <- ret@names
	names(ret@stderr.coverage) <- ret@names
	names(ret@freqt.bias) <- ret@names
	names(ret@freqt.typeI) <- ret@names
	names(ret@freqt.typeII) <- ret@names
	names(ret@bayes.bias) <- ret@names
	names(ret@bayes.typeI) <- ret@names
	names(ret@bayes.typeII) <- ret@names

	return(ret)
})

# Generate a function to perform simulations with sample size up to max.n with fixed independent variables
# Model is for Mendelian randomization with the specified causal effects
gen.simulate = function(data, max.n=1000, seed=NA) {
	stopifnot(is(data, "bmsim_data"))
	validate(data)
	x.mu = colMeans(data@x)
	x.sigma = cov(data@x)
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
	ret@freqt.bias = rowMeans(sapply(1:nsim, function(i) performance.list[[i]]@freqt.bias))
	ret@freqt.typeI = rowMeans(sapply(1:nsim, function(i) performance.list[[i]]@freqt.typeI))
	ret@freqt.typeII = rowMeans(sapply(1:nsim, function(i) performance.list[[i]]@freqt.typeII))
	ret@freqt.familywiseI = mean(sapply(1:nsim, function(i) performance.list[[i]]@freqt.familywiseI))
	ret@freqt.familywiseII = mean(sapply(1:nsim, function(i) performance.list[[i]]@freqt.familywiseII))
	ret@freqt.fdr = mean(sapply(1:nsim, function(i) performance.list[[i]]@freqt.fdr))
	ret@freqt.fndr = mean(sapply(1:nsim, function(i) performance.list[[i]]@freqt.fndr))
	ret@freqt.fpr = mean(sapply(1:nsim, function(i) performance.list[[i]]@freqt.fpr))
	ret@freqt.fnr = mean(sapply(1:nsim, function(i) performance.list[[i]]@freqt.fnr))

	ret@bayes.bias = rowMeans(sapply(1:nsim, function(i) performance.list[[i]]@bayes.bias))
	ret@bayes.typeI = rowMeans(sapply(1:nsim, function(i) performance.list[[i]]@bayes.typeI))
	ret@bayes.typeII = rowMeans(sapply(1:nsim, function(i) performance.list[[i]]@bayes.typeII))
	ret@bayes.familywiseI = mean(sapply(1:nsim, function(i) performance.list[[i]]@bayes.familywiseI))
	ret@bayes.familywiseII = mean(sapply(1:nsim, function(i) performance.list[[i]]@bayes.familywiseII))
	ret@bayes.fdr = mean(sapply(1:nsim, function(i) performance.list[[i]]@bayes.fdr))
	ret@bayes.fndr = mean(sapply(1:nsim, function(i) performance.list[[i]]@bayes.fndr))
	ret@bayes.fpr = mean(sapply(1:nsim, function(i) performance.list[[i]]@bayes.fpr))
	ret@bayes.fnr = mean(sapply(1:nsim, function(i) performance.list[[i]]@bayes.fnr))

	# Prediction
	ret@prediction.L1 = mean(sapply(1:nsim, function(i) performance.list[[i]]@prediction.L1))
	ret@prediction.L2 = mean(sapply(1:nsim, function(i) performance.list[[i]]@prediction.L2))

	# Computation time (seconds)
	ret@time.secs = mean(sapply(1:nsim, function(i) performance.list[[i]]@time.secs))

	# Apply names
	names(ret@estimate.bias) <- ret@names
	names(ret@stderr.coverage) <- ret@names
	names(ret@freqt.bias) <- ret@names
	names(ret@freqt.typeI) <- ret@names
	names(ret@freqt.typeII) <- ret@names
	names(ret@bayes.bias) <- ret@names
	names(ret@bayes.typeI) <- ret@names
	names(ret@bayes.typeII) <- ret@names

	return(ret)
}

# Perform a standardized set of analyses of a single dataset
do.analyses = function(data, nu=data@m, dblthk.h = c(0.25, 1, 4), dblthk.mu = c(0.05, 0.1, 0.2), mr.bma.nsim=1000, mr.bma.sigma=0.5, mr.bma.prior_prob=0.1) {
	stopifnot(is(data, "bmsim_data"))
	validate(data)
	
	results = list()
	# Single model results based on various types of model selection
	results[["grand alternative"]] = add1drop1(data, rep(TRUE, data@m), nu, analysis.name="Multivariable Mendelian randomization with all variables; leave one out/add one in significance testing")

	# Stepwise regression
	results[["backward elimination"]] = do.stepwise(data, nu=nu)
	
	# Fit LASSO, elastic net and ridge regression
	results.glmnet = do.glmnet(data, nu=nu)
	results[["lasso"]] = results.glmnet[["lasso"]]
	results[["elastic net"]] = results.glmnet[["elnet"]]
	results[["ridge regression"]] = results.glmnet[["ridge"]]
	# LASSO-based model selection
	results[["lasso model selection"]] = add1drop1(data, results.glmnet$lasso@estimate!=0, nu, analysis.name="Multivariable Mendelian randomization with LASSO; leave one out/add one in significance testing")
	# Elastic net-based model selection
	results[["elastic net model selection"]] = add1drop1(data, results.glmnet$elnet@estimate!=0, nu, analysis.name="Multivariable Mendelian randomization with elastic net; leave one out/add one in significance testing")

	# Perform Bayesian model-averaged Mendelian randomization using Doublethink
	results.doublethink = doublethink.x(data, h=dblthk.h, mu=dblthk.mu, nu=nu)
	doublethink.names = outer(dblthk.h, dblthk.mu, function(h,mu) paste0("mu = ", mu, "; h = ", h))
	for(doublethink.name in doublethink.names) {
		# Doublethink BMA
		results[[paste0("doublethink bma ", doublethink.name)]] = results.doublethink$doublethink.bma[[doublethink.name]]
		# Doublethink MAP-based model selection
		results[[paste0("doublethink model selection ", doublethink.name)]] = results.doublethink$doublethink.modelselection[[doublethink.name]]
	}

	# Perform Bayesian model-averaged Mendelian randomization using MR-BMA (no permutation procedure)
	results.mr.bma = mr.bma.x(data, sigma=mr.bma.sigma, prior_prob=mr.bma.prior_prob, nsim=0, nu=nu)
	results[["mr-bma bma"]] = results.mr.bma$bma
	# results[["mr-bma map"]] = results.mr.bma$map
	# MR-BMA MAP-based model selection
	results[["mr-bma model selection"]] = add1drop1(data, results.mr.bma$map@estimate!=0, nu, analysis.name="Bayesian model selection multivariable Mendelian randomization with MR-BMA; leave one out/add one in significance testing")

	# Perform Bayesian model-averaged Mendelian randomization using MR-BMA (with permutation procedure)
	# (Do this separately for timing purposes)
	results.mr.bma = mr.bma.x(data, sigma=mr.bma.sigma, prior_prob=mr.bma.prior_prob, nsim=mr.bma.nsim, nu=nu)
	results[["mr-bma bma permute"]] = results.mr.bma$bma
	
	return(results)
}

# Perform a standardized set of analyses of a single dataset
calc.performance = function(analyses, params, freqt.alpha = 0.05, bayes.tau = 19) {
	stopifnot(is.list(analyses))
	nanal = length(analyses)
	for(analysis in analyses) {
		stopifnot(is(analysis, "bmsim_analysisResults"))
		validate(analysis)
		stopifnot(analysis@m==length(params))
	}
	thresh.neglog10padj = -log10(freqt.alpha)
	thresh.log10po = log10(bayes.tau)
	ret = list()
	for(i in 1:nanal) {
		analysis.name = names(analyses)[i]
		ret[[analysis.name]] = performance(analyses[[analysis.name]], params, thresh.neglog10padj=thresh.neglog10padj, thresh.log10po=thresh.log10po)
	}
	return(ret)
}
