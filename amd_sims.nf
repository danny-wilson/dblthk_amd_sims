process simulate_dry_run {
	input:
		val taskid
		tuple val(params), val(n)
	output:
		tuple val(params), val(n)
		path("taskid"), emit: taskid
		path("sim.data.RDS"), emit: sim.data
		path("sim.analyses.RDS"), emit: sim.analyses
		path("sim.perf.RDS"), emit: sim.perf
shell:
'''
	echo !{taskid} !{params} !{n}
'''
}

process simulate {
beforeScript "module add R/4.2.1-foss-2022a"
input:
	val taskid
	tuple val(params), val(n)
output:
	tuple val(params), val(n)
	path("${taskid}.sim.data.RDS"), 	emit: sim_data
	path("${taskid}.sim.anal.RDS"),		emit: sim_anal
	path("${taskid}.sim.perf.RDS"), 	emit: sim_perf
shell:
'''
	#!/usr/bin/env Rscript

	# Get working directory
	wd = getwd()
	setwd("/well/bag/wilson/GitHub/dblthk_amd_sims")
	# Load required source code: this in turn sources summary_mvMR_BF.R
	source("biomarker-sim-functions.R")

	# Read arguments
	taskid = as.integer("!{taskid}")
	params = as.numeric("!{params}")
	n = as.integer("!{n}")
	filename.sim_data = "!{taskid}.sim.data.RDS"
	filename.sim_anal = "!{taskid}.sim.anal.RDS"
	filename.sim_perf = "!{taskid}.sim.perf.RDS"

	# Check arguments
	stopifnot(!is.na(taskid))
	stopifnot(!any(is.na(params)))
	stopifnot(length(params)==1)
	params = rep(params, 15)
	stopifnot(n>15)

	# Load example AMD data
	full.data = load.data()

	# Return to working directory
	setwd(wd)

	# Perform univariable association tests
	results.univariable.MR = univariable.MR(full.data)

	# Filter example data to the 15 most significant risk factors based on univariable associations
	data = filter.data(full.data, as.integer(15), results.univariable.MR)

	# Set seed to fix the independent variables (x)
	set.seed(0)
	# Simulate the independent variables based on the example AMBdata
	simulate = gen.simulate(data, n)
	
	# Set seed again to unfix the dependent variables (y)
	set.seed(taskid)
	# Simulate the dependent variables
	sim.data = simulate(n, params)
	
	# Perform the standard set of analyses
	sim.anal = do.analyses(sim.data, nu=full.data@m, mr.bma.nsim=0)

	# Compute performance metrics for the analyses
	sim.perf = calc.performance(sim.anal, params, freqt.alpha=0.05, bayes.tau=19)

	# Save files
	saveRDS(sim.data, file=filename.sim_data)
	saveRDS(sim.anal, file=filename.sim_anal)
	saveRDS(sim.perf, file=filename.sim_perf)
	
	cat("simulate: Doublethink AMD simulation completed successfully\n")
'''
}

process combine_performance {
	beforeScript "module add R/4.2.1-foss-2022a"
	input:
		path(infiles)
	output:
		path("combined.performance.RDS"), emit: combined_perf
	shell:
'''
	#!/usr/bin/env Rscript
	
	# Get working directory
	wd = getwd()
	setwd("/well/bag/wilson/GitHub/dblthk_amd_sims")
	# Load required source code: this in turn sources summary_mvMR_BF.R
	source("biomarker-sim-functions.R")
	# Return to working directory
	setwd(wd)

	# Read arguments
	infiles.all = "!{infiles}"
	outfile = "combined.performance.RDS"
	
	# Check arguments
	infiles = unlist(strsplit(infiles.all, ' '))
	nfiles = length(infiles)
	stopifnot(file.exists(infiles))
	stopifnot(!file.exists(outfile))
	cat("combine_performance: combining", nfiles, "files\n")
	
	# Load simulation performance metrics
	perf = list()
	nanalyses = list()
	for(i in 1:nfiles) {
		perf[[i]] = readRDS(infiles[i])
		nanalyses[[i]] = length(perf[[i]])
		stopifnot(nanalyses[[i]]==nanalyses[[1]])
		for(j in 1:nanalyses[[1]]) stopifnot(is(perf[[i]][[j]], "bmsim_performance"))
	}
	
	# Combine analysis results
	results = list()
	for(i in 1:nanalyses[[1]]) {
		# Combine the performance metrics
		analysis.name = names(perf[[1]])[i]
		results[[analysis.name]] = combine.performance(lapply(perf, function(PERF) PERF[[analysis.name]]))
	}
	
	# Save results
	saveRDS(results, file="combined.performance.RDS")

	cat("combine_performance: Completed successfully\n")
'''
}

// Default values
params.ntasks = 2
params.params = [0.0] * params.ntasks
params.n = [149] * params.ntasks

// Define the channels
ch_taskid = Channel.of(1..params.ntasks)
ch_params = Channel.fromList(params.params)
ch_n = Channel.fromList(params.n)
ch_args = ch_params.combine(ch_n)

// Print arguments
println 'Arguments'
println 'ntasks:			' + params.ntasks
println 'params:			' + params.params
println 'n:					' + params.n

// Go
workflow {
	// Perform the simulations
	simulate(ch_taskid, ch_args)

	// Combine performance metrics
	combine_performance(simulate.out.sim_perf.collect())
}
