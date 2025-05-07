# High-performance compute cluster commands (run in interactive mode):
# srun --pty --mem=320GB /bin/bash
# module add R/4.3.2-gfbf-2023a
# OMP_NUM_THREADS=20 R

#############
# Arguments #
#############
# Working directory
workdir = ""
# Local clone of https://github.com/danny-wilson/dblthk_amd_sims
dblthk_amd_sims_repo = "~/GitHub/dblthk_amd_sims/"
# Set the sample size independently of the example data
nsamp = 145
# Number of simulations
nsim = ceiling(50000 * sqrt(nsamp/145))
# Hyper-parameters
mu = 0.1; h = 1; hsim = h; xi = h/(nsamp+h)
logc = log(mu * sqrt(h/(nsamp+h)))
tau = 9
# Filenames
output_filename = paste0(workdir, "/nsamp", nsamp, "_nsim", nsim, "_mu", mu, "_h", h, "_", format(Sys.time(), "%Y-%M-%d"))

###########################
# Libraries and functions #
###########################
require(mvtnorm)
# logarithm of the sum of logged variables (aka logsumexp)
logsum = function(x) { mx = max(x); mx + log(sum(exp(x-mx))) }
# Simulate from the prior. Relies on globals data, L, s, degfree, mu, h, (inputs), prior.var (output)
prior.var = list()
prior.sim.correct = function(nsim) {
	prior.var = list()
	wh.s = sample(L, nsim, replace=TRUE, prob=mu^degfree)
	for(WH in unique(wh.s)) {
		if(WH==1) prior.var[[WH]] = 0.0
		do.calc = length(prior.var)<WH
		if(!do.calc) do.calc = is.null(prior.var[[WH]])
		if(do.calc) prior.var[[WH]] = solve(UFI[s[WH,], s[WH,]])/h
	}
	ret = matrix(0.0, nsim, data@m)
	for(i in 1:nsim) {
		WH = wh.s[i]
		if(WH>1) {
			S = s[WH,]
			ret[i, S] <- rmvnorm(1, sigma=prior.var[[WH]])
		}
	}
	return(ret)
}
# Perform simulations. Relies on globals: nsamp, UFI.inv, L, s, logc, xi, ret.beta.hat (output)
# Later: ret.beta.hat = matrix(0.0, 0, data@m)
do.simulation = function(all_params, fwel.threshold = 1.0, record.marginals=FALSE) {
	nsim = nrow(all_params)
	
	# Simulate beta.hat directly under the grand alternative
	beta.hat = t(all_params + rmvnorm(nsim, sigma=1.0/nsamp * UFI.inv))
	ret.beta.hat = rbind(ret.beta.hat, t(beta.hat))

	# Compute the loglikelihood of every model for every simulation
	system.time((loglik = t(vapply(1:L, function(i) {
		# U = beta.hat.multiplier[[i]] %*% beta.hat
		0.5 * colSums((beta.hat.multiplier[[i]] %*% beta.hat)^2)
	}, FUN.VALUE = rep(0.0, nsim)))))
	#   user  system elapsed
	#288.789  23.517  27.672

	# Convert that into a log posterior odds (under fixed mu, h)
	logPO = degfree * logc + (1.0-xi) * loglik


	####################
	# Hypothesis tests #
	####################

	system.time((logPO.headline = apply(logPO[-1,], 2, logsum)))
	#  user  system elapsed
	#11.043   3.061  14.105

	ret = cbind("logPO.headline" = logPO.headline)
	if(record.marginals) {
		tPP = exp(t(logPO) - apply(logPO, 2, max))
		marg = log(tPP %*% s) - log(tPP %*% !s)
		colnames(marg) <- paste0("logPO.marginal.", 1:ncol(s))
		ret = cbind(ret, marg)
	}

	for(FWEL.THRESHOLD in fwel.threshold) {
		is.nul = ((all_params!=0.0) %*% (data.rsq>=FWEL.THRESHOLD))==0.0

		# Number of true null hypotheses
		ntruenull = rowSums(is.nul)

		# Compute the doublethink test statistic
		system.time((logPO.truenull = t(vapply(1:nsim, function(i) {
			if(ntruenull[i]==0) return(as.double(c(NA, NA)))
			model.truenull = rowSums(s[, is.nul[i,], drop=FALSE])==0
			c(logsum(logPO[!model.truenull, i]), logsum(logPO[model.truenull, i]))
		}, FUN.VALUE=double(2)))))
		# user  system elapsed
		#3.886   0.542   4.428

		# Compute the asymptotically equivalent test statistic (for falsifiability of assumptions)
		system.time((logPO.asymnull = t(vapply(1:nsim, function(i) {
			ntrue = ntruenull[i]
			if(ntrue==0) return(as.double(c(NA, NA)))
			model.null = colSums(is.nul[i,]!=t(s))==data@m
			model.alts = (rowSums(s[, !is.nul[i,], drop=FALSE])==(data@m-ntrue)) & !model.null
			c(logsum(logPO[model.alts, i]), logsum(logPO[model.null, i]))
		}, FUN.VALUE=double(2)))))

		ret = cbind(ret, logPO.truenull, logPO.asymnull, ntruenull)
		colnames(ret)[ncol(ret)-(4:0)] <- paste0(c("logPO.truenull.num.", "logPO.truenull.den.", "logPO.asymnull.num.", "logPO.asymnull.den.", "ntruenull."), FWEL.THRESHOLD)
	}
	return(ret)
}
# Number of elementary hypotheses at different fwel tolerances
fwel.groups = Vectorize(function(FWEL.THRESHOLD) {
	in.gp = data.rsq>=FWEL.THRESHOLD
	gp = 1:data@m
	gp[in.gp[1,]] = gp[1]
	for(j in 2:data@m) {
		# Find the current group numbers of all individuals in a group with j
		gps.to.merge = unique(gp[in.gp[j,]])
		# Merge all those groups with j, not just those individuals in a group with j
		gp[!is.na(match(gp, gps.to.merge))] = gp[j]
	}
	return(gp)
})

######################
# Load code and data #
######################

# Load code from repo then set working directory
setwd(dblthk_amd_sims_repo)
source("biomarker-sim-functions.R")
full.data = load.data()
setwd(workdir)

#######################
# Do pre-computation  #
#######################

# Perform univariable association tests
results.univariable.MR = univariable.MR(full.data)

# Filter example data to the 15 most significant risk factors based on univariable associations
data = filter.correlated.data(full.data, m=as.integer(15), rsq.max.thresh=0.66, results.univariable.MR=results.univariable.MR)
data.rsq = cor(data@x)^2

# FWEL tolerances: pre-specified gradient of FWEL tolerances defining different group sizes
fwel.tol = sort(unique(data.rsq[lower.tri(data.rsq)]), decreasing=TRUE)
gp = fwel.groups(fwel.tol)
ngps = apply(gp, 2, function(GP) length(unique(GP)))
# First occurence of each number of groups
fwel.tol = c(0, fwel.tol[match(2:14, ngps)], 1)
gp = fwel.groups(fwel.tol)
ngps = apply(gp, 2, function(GP) length(unique(GP)))

# Total number of additive models
(L = as.integer(2^data@m))
# [1] 32768
# Model inclusion matrix
s = t(sapply(0:(L-1), function(s) {
	(incl = intToBits(s)[data@m:1]==1)
})); colnames(s) <- data@x.names

# Unit Fisher information matrix assuming known error variance equal to 1
ssq = 1
UFI = crossprod(data@x)/nrow(data@x)/ssq
UFI.inv = solve(UFI)

# See OneNote 'More efficient simulations' 6/6/24, variable M^{(s)}
system.time((beta.hat.multiplier = lapply(1:L, function(i) {
	S = s[i,]
	norm.S = sum(S)
	if(norm.S==0) return(matrix(0.0, nrow=1, ncol=data@m))
	if(norm.S==data@m) return(diag(data@m))
	I.S = diag(norm.S)
	A.S = - UFI.inv[S, !S] %*% solve(UFI.inv[!S, !S])
	# NB: t(chol(X)) %*% chol(X) == X for symmetric X
	M.S = sqrt(nsamp) * ( chol(UFI[S, S]) %*% cbind(I.S, A.S) )
	# E.g. if S = c(T, F, T, F) then od.S = c(1, 3, 2, 4)
	# E.g. S = c(F, F, ..., T, F) and od.S = c(14, 1, 2, ..., 13, 15)
	od.S = order(S, decreasing=TRUE)
	# Re-order the columns of (I.S, A.S) to align with the rows of beta.hat
	# O.S = matrix(0.0, data@m, data@m); O.S[cbind(1:data@m, od.S)] <- 1.0
	# M.S %*% O.S
	# Faster:
	rk.S = order(od.S)
	M.S[, rk.S, drop=FALSE]
})))
#   user  system elapsed 
#  5.192   0.603   5.816 

# For obtaining point estimates (15 x 15 x 32768 array)
system.time((beta.hat.converter = vapply(1:L, function(i) {
	S = s[i,]
	norm.S = sum(S)
	if(norm.S==0) return(matrix(0.0, nrow=data@m, ncol=data@m))
	if(norm.S==data@m) return(diag(data@m))
	I.S = diag(norm.S)
	A.S = - UFI.inv[S, !S] %*% solve(UFI.inv[!S, !S])
	M.S = cbind(I.S, A.S)
	od.S = order(S, decreasing=TRUE)
	rk.S = order(od.S)
	ret = matrix(0.0, data@m, data@m)
	ret[S,] = M.S[, rk.S, drop=FALSE]
	return(ret)
}, FUN.VALUE = matrix(0.0, data@m, data@m))))
#   user  system elapsed 
#  4.689   0.543   5.254 

# Degrees of freedom for each model
degfree = rowSums(s)

#######################
# Perform simulations #
#######################

h.bkp = h
h = hsim
all_params = prior.sim.correct(nsim)
h = h.bkp
stopifnot(nrow(all_params)==nsim)
stopifnot(ncol(all_params)==data@m)

# Negligible speed improvement seen above batch size of 5000
# 5000 is 20% faster than 1000
# 1000 is 100% faster than 100
# Max out at about 100 simulations per second, using ~1 CPU and
# memory: RES 14.5g VIRT 16.9g (for blocksize 10000)
blocksize = 5000
nblock = ceiling(nsim/blocksize)
ret.beta.hat = matrix(0.0, 0, data@m) # Storage for beta.hat
start_time = Sys.time()
for(block in 1:nblock) {
	if(block==1) {
		sim = do.simulation(all_params[1:pmin(nsim, blocksize),], fwel.tol)
	} else {
		sim = rbind(sim, do.simulation(all_params[unique(pmin(nsim, (block-1)*blocksize + (1:blocksize))),], fwel.tol))
	}
	check_time = Sys.time()-start_time
	cat("Done", pmin(nsim, block*blocksize), "simulations in", check_time, units(check_time), "\n")
}
# Done 1e+05 simulations in 14.14855 mins

################
# Save results #
################
colnames(all_params) <- paste0("true.", data@x.names)
stopifnot(!file.exists(output_filename)); saveRDS(cbind(sim, all_params), output_filename)
